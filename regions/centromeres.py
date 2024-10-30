from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import numpy as np
from pathlib import Path
import logging

from .common import (
    ChromosomeRegion, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator, SequenceFeature
)
from repeats import (
    RepeatGenerator, RepeatType, HORGenerator,
    HomogenizationEngine, HORParams
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CentromereZone(Enum):
    """Типы зон центромерного региона"""
    CORE = "core"  # Центральная область с высокоорганизованными HOR
    PERIPHERAL = "peripheral"  # Периферийные области с вариабельными повторами
    TRANSITION = "transition"  # Переходные зоны

@dataclass
class CENPBParams:
    """Параметры для CENP-B боксов"""
    consensus: str = "TTCGTTGGAAACGGGA"  # Консенсусная последовательность
    spacing: int = 171  # Расстояние между боксами (обычно равно размеру мономера)
    mutation_rate: float = 0.1  # Частота мутаций в боксах
    presence_rate: float = 0.8  # Частота присутствия бокса в допустимых позициях

@dataclass
class CentromereParams:
    """Параметры для генерации центромерного региона"""
    min_core_length: int = 1500000  # 1.5 Mb
    max_core_length: int = 3000000  # 3 Mb
    min_peripheral_length: int = 500000  # 500 kb
    max_peripheral_length: int = 1000000  # 1 Mb
    transition_length: int = 200000  # 200 kb
    
    # Параметры HOR
    core_hor_params: HORParams = None
    peripheral_hor_params: HORParams = None
    
    # Параметры CENP-B боксов
    cenpb_params: CENPBParams = None
    
    def __post_init__(self):
        # Устанавливаем параметры по умолчанию для HOR в core регионе
        if self.core_hor_params is None:
            self.core_hor_params = HORParams(
                unit_size=171,  # Размер альфа-сателлитного мономера
                units_in_hor=12,  # Стандартное количество мономеров в HOR
                hor_copies=100,  # Большое количество копий для core региона
                mutation_rate_between_units=0.02,
                mutation_rate_between_hors=0.01,
                gc_content=0.58
            )
        
        # Параметры HOR для периферийных областей (более вариабельные)
        if self.peripheral_hor_params is None:
            self.peripheral_hor_params = HORParams(
                unit_size=171,
                units_in_hor=8,  # Меньше мономеров в HOR
                hor_copies=50,  # Меньше копий
                mutation_rate_between_units=0.05,  # Выше уровень мутаций
                mutation_rate_between_hors=0.03,
                gc_content=0.56
            )
        
        # Параметры CENP-B боксов
        if self.cenpb_params is None:
            self.cenpb_params = CENPBParams()

class CentromereRegion(ChromosomeRegion):
    """Класс для генерации центромерного региона"""
    
    def __init__(self, params: RegionParams, 
                 centromere_params: Optional[CentromereParams] = None):
        super().__init__(params)
        self.centromere_params = centromere_params or CentromereParams()
        self.repeat_gen = RepeatGenerator()
        self.hor_gen = HORGenerator()
        self.gradient_gen = GradientGenerator()
        self.homogenization = HomogenizationEngine()
    
    def _generate_cenpb_boxes(self, sequence: str, zone: CentromereZone) -> List[Dict]:
        """Генерирует и вставляет CENP-B боксы"""
        boxes = []
        params = self.centromere_params.cenpb_params
        
        # Определяем вероятность вставки бокса в зависимости от зоны
        if zone == CentromereZone.CORE:
            presence_prob = params.presence_rate
        elif zone == CentromereZone.PERIPHERAL:
            presence_prob = params.presence_rate * 0.7  # Меньше боксов на периферии
        else:
            presence_prob = params.presence_rate * 0.4  # Еще меньше в переходных зонах
        
        # Проходим по последовательности с шагом spacing
        for pos in range(0, len(sequence), params.spacing):
            if np.random.random() < presence_prob:
                box_seq = list(params.consensus)
                
                # Вносим мутации
                for i in range(len(box_seq)):
                    if np.random.random() < params.mutation_rate:
                        original_base = box_seq[i]
                        bases = ['A', 'T', 'G', 'C']
                        bases.remove(original_base)
                        box_seq[i] = np.random.choice(bases)
                
                boxes.append({
                    'type': 'CENP-B_box',
                    'start': pos,
                    'end': pos + len(params.consensus),
                    'sequence': ''.join(box_seq),
                    'zone': zone.value
                })
        
        return boxes
    
    def _generate_zone(self, zone: CentromereZone, length: int) -> Tuple[str, List[Dict]]:
        """Генерирует последовательность для определенной зоны центромеры"""
        # Выбираем параметры в зависимости от зоны
        if zone == CentromereZone.CORE:
            hor_params = self.centromere_params.core_hor_params
        else:
            hor_params = self.centromere_params.peripheral_hor_params
        
        self.hor_gen.params = hor_params
        
        # Генерируем последовательность HOR
        sequence, mutations = self.hor_gen.generate_array()
        
        # Обрезаем или дополняем до нужной длины
        if len(sequence) < length:
            copies_needed = (length - len(sequence)) // len(sequence) + 1
            sequence = sequence * copies_needed
        sequence = sequence[:length]
        
        # Добавляем CENP-B боксы
        cenpb_boxes = self._generate_cenpb_boxes(sequence, zone)
        
        # Применяем гомогенизацию к core региону
        if zone == CentromereZone.CORE:
            sequence, conversions = self.homogenization.apply_gene_conversion(
                sequence,
                unit_size=hor_params.unit_size,
                units_in_hor=hor_params.units_in_hor
            )
        
        # Собираем все особенности
        features = []
        features.extend([{**mut, 'zone': zone.value} for mut in mutations])
        features.extend(cenpb_boxes)
        if zone == CentromereZone.CORE:
            features.extend([{**conv, 'zone': zone.value} for conv in conversions])
        
        return sequence, features
    
    def _create_transition_zone(self, zone1_seq: str, zone2_seq: str,
                              length: int) -> Tuple[str, List[Dict]]:
        """Создает переходную зону между двумя регионами"""
        # Создаем градиент для смешивания последовательностей
        gradient_params = GradientParams(
            start_value=0.0,
            end_value=1.0,
            shape='sigmoid',
            steepness=5.0
        )
        gradient = self.gradient_gen.create_gradient(gradient_params, length)
        
        # Генерируем переходную последовательность
        transition_seq = ""
        for i in range(length):
            # Определяем, какую последовательность использовать в данной позиции
            if np.random.random() < gradient[i]:
                base_pool = zone2_seq
            else:
                base_pool = zone1_seq
            
            # Выбираем случайный нуклеотид из соответствующей последовательности
            pos = np.random.randint(0, len(base_pool))
            transition_seq += base_pool[pos]
        
        # Создаем аннотации для переходной зоны
        features = [{
            'type': 'transition',
            'start': 0,
            'end': length,
            'gradient_value': float(g),
            'zone': CentromereZone.TRANSITION.value
        } for g in gradient]
        
        return transition_seq, features
    
    def generate(self) -> Tuple[str, List[Dict]]:
        """Генерирует полный центромерный регион"""
        # Определяем размеры каждой зоны
        core_length = np.random.randint(
            self.centromere_params.min_core_length,
            self.centromere_params.max_core_length
        )
        left_peripheral_length = np.random.randint(
            self.centromere_params.min_peripheral_length,
            self.centromere_params.max_peripheral_length
        )
        right_peripheral_length = np.random.randint(
            self.centromere_params.min_peripheral_length,
            self.centromere_params.max_peripheral_length
        )
        
        # Генерируем основные зоны
        core_seq, core_features = self._generate_zone(CentromereZone.CORE, core_length)
        left_peri_seq, left_features = self._generate_zone(
            CentromereZone.PERIPHERAL,
            left_peripheral_length
        )
        right_peri_seq, right_features = self._generate_zone(
            CentromereZone.PERIPHERAL,
            right_peripheral_length
        )
        
        # Создаем переходные зоны
        left_trans_seq, left_trans_features = self._create_transition_zone(
            left_peri_seq, core_seq,
            self.centromere_params.transition_length
        )
        right_trans_seq, right_trans_features = self._create_transition_zone(
            core_seq, right_peri_seq,
            self.centromere_params.transition_length
        )
        
        # Собираем полную последовательность
        sequence = (left_peri_seq + left_trans_seq + 
                   core_seq + right_trans_seq + right_peri_seq)
        
        # Корректируем позиции в аннотациях и собираем их вместе
        current_pos = 0
        features = []
        
        # Левая периферийная область
        for feature in left_features:
            feature['start'] += current_pos
            feature['end'] += current_pos
            features.append(feature)
        current_pos += len(left_peri_seq)
        
        # Левая переходная зона
        for feature in left_trans_features:
            feature['start'] += current_pos
            feature['end'] += current_pos
            features.append(feature)
        current_pos += len(left_trans_seq)
        
        # Центральная область
        for feature in core_features:
            feature['start'] += current_pos
            feature['end'] += current_pos
            features.append(feature)
        current_pos += len(core_seq)
        
        # Правая переходная зона
        for feature in right_trans_features:
            feature['start'] += current_pos
            feature['end'] += current_pos
            features.append(feature)
        current_pos += len(right_trans_seq)
        
        # Правая периферийная область
        for feature in right_features:
            feature['start'] += current_pos
            feature['end'] += current_pos
            features.append(feature)
        
        return sequence, features