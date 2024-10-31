# yatomat/regions/centromeres.py
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import numpy as np
from pathlib import Path
import logging

from yatomat.regions.common import (
    ChromosomeRegion, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator, SequenceFeature
)
from yatomat.repeats import (
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

    core_hor_params: Optional[HORParams] = None
    peripheral_hor_params: Optional[HORParams] = None
    cenpb_params: Optional[CENPBParams] = None

    def __post_init__(self):
        if self.core_hor_params is None:
            self.core_hor_params = HORParams(
                unit_size=171,
                units_in_hor=12,
                hor_copies=100,
                mutation_rate_between_units=0.02,
                mutation_rate_between_hors=0.01,
                gc_content=0.58
            )

        if self.peripheral_hor_params is None:
            self.peripheral_hor_params = HORParams(
                unit_size=171,
                units_in_hor=8,
                hor_copies=50,
                mutation_rate_between_units=0.05,
                mutation_rate_between_hors=0.03,
                gc_content=0.56
            )

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
        """Generates and inserts CENP-B boxes"""
        if zone != CentromereZone.CORE:
            return []  # Only generate CENP-B boxes for the core region

        boxes = []
        params = self.centromere_params.cenpb_params

        if params is None:
            return boxes

        # Determine the probability of box insertion based on the zone
        presence_prob = params.presence_rate

        # Iterate through the sequence with a step of spacing
        for pos in range(0, len(sequence), params.spacing):
            if np.random.random() < presence_prob:
                # Split the consensus sequence into important regions
                box_seq = list(params.consensus)
                # TTCG and GGGA - important motifs to preserve
                conserved_regions = [(0, 4), (12, 16)]  # positions of TTCG and GGGA

                # Introduce mutations only in non-conserved regions
                for i in range(len(box_seq)):
                    if not any(start <= i < end for start, end in conserved_regions):
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

    def _convert_mutation_to_feature(self, mutation: Dict, position: int, zone: CentromereZone) -> Dict:
        """Преобразует мутацию в полноценный feature с необходимыми полями"""
        feature = {
            'type': mutation['type'],
            'start': position + mutation.get('position', 0),
            'end': position + mutation.get('position', 0) + mutation.get('size', 1),
            'zone': zone.value
        }
        return feature

    def _generate_zone(self, zone: CentromereZone, target_length: int) -> Tuple[str, List[Dict]]:
        """Генерирует последовательность для определенной зоны центромеры"""
        if target_length <= 0:
            return "", []

        if zone == CentromereZone.CORE:
            hor_params = self.centromere_params.core_hor_params
        else:
            hor_params = self.centromere_params.peripheral_hor_params

        if hor_params is None:
            return "", []

        # Рассчитываем базовые параметры
        single_hor_length = hor_params.unit_size * hor_params.units_in_hor
        min_required_copies = (target_length + single_hor_length - 1) // single_hor_length

        # Создаем новые параметры с нужным количеством копий
        new_params = HORParams(
            unit_size=hor_params.unit_size,
            units_in_hor=hor_params.units_in_hor,
            hor_copies=min_required_copies,
            mutation_rate_between_units=hor_params.mutation_rate_between_units,
            mutation_rate_between_hors=hor_params.mutation_rate_between_hors,
            gc_content=hor_params.gc_content
        )

        # Генерируем HOR массив
        self.hor_gen.params = new_params
        sequence, mutations = [], []

        while len(sequence) < target_length:
            seq, muts = self.hor_gen.generate_array()
            current_len = len(sequence)

            # Корректируем позиции мутаций
            for mut in muts:
                if 'position' in mut:
                    mut['position'] += current_len

            sequence.extend(list(seq))
            mutations.extend(muts)

        # Обрезаем до нужной длины
        sequence = ''.join(sequence[:target_length])

        # Создаем features
        features = []

        # Добавляем базовый feature для всей зоны
        features.append({
            'type': f'{zone.value}_region',
            'start': 0,
            'end': target_length,
            'zone': zone.value
        })

        # Добавляем мутации
        for mut in mutations:
            if mut.get('position', 0) < target_length:
                features.append(self._convert_mutation_to_feature(mut, 0, zone))

        # Добавляем CENP-B боксы
        features.extend(self._generate_cenpb_boxes(sequence, zone))

        # Применяем гомогенизацию для core региона
        if zone == CentromereZone.CORE:
            sequence, conversions = self.homogenization.apply_gene_conversion(
                sequence,
                unit_size=hor_params.unit_size,
                units_in_hor=hor_params.units_in_hor
            )

            for conv in conversions:
                if conv['start'] < target_length:
                    features.append({
                        'type': 'gene_conversion',
                        'start': conv['start'],
                        'end': min(conv['start'] + conv['length'], target_length),
                        'donor_start': conv['donor_start'],
                        'zone': zone.value
                    })

        logger.info(f"Generated {zone.value} zone of length {len(sequence)} (target: {target_length})")
        return sequence, features

    def _create_transition_zone(self, zone1_seq: str, zone2_seq: str,
                                length: int, reverse_gradient: bool = False) -> Tuple[str, List[Dict]]:
        """Создает переходную зону между двумя регионами"""
        if length <= 0:
            return "", []

        # Создаем более плавный сигмоидальный градиент
        positions = np.linspace(-6, 6, length)  # Используем более широкий диапазон для плавности
        gradient = 1 / (1 + np.exp(-positions))
        gradient = (gradient - gradient.min()) / (gradient.max() - gradient.min())

        if reverse_gradient:
            gradient = 1 - gradient  # Reverse the gradient

        # Генерируем последовательность
        sequence = []

        # Берем случайные участки из обеих последовательностей
        window_size = min(171, len(zone1_seq), len(zone2_seq))  # Используем размер альфа-сателлита

        for i in range(0, length, window_size):
            current_pos = min(i + window_size, length)
            grad_value = gradient[i:current_pos].mean()

            # Выбираем последовательность на основе градиента
            if np.random.random() < grad_value:
                source_seq = zone2_seq
            else:
                source_seq = zone1_seq

            # Берем случайный участок соответствующей длины
            start_pos = np.random.randint(0, len(source_seq) - window_size)
            sequence.extend(source_seq[start_pos:start_pos + window_size])

        sequence = ''.join(sequence[:length])

        # Создаем features для визуализации градиента
        features = []

        # Add a feature for the entire transition zone
        features.append({
            'type': f'{CentromereZone.TRANSITION.value}_region',
            'start': 0,
            'end': length,
            'zone': CentromereZone.TRANSITION.value
        })

        step = length // 20  # 20 окон для плавного отображения
        for i in range(0, length, step):
            end_pos = min(i + step, length)
            grad_value = float(gradient[i:end_pos].mean())

            features.append({
                'type': 'transition',
                'start': i,
                'end': end_pos,
                'gradient_value': grad_value,
                'zone': CentromereZone.TRANSITION.value
            })

        logger.info(f"Generated transition zone of length {len(sequence)} with {len(features)} gradient features")

        # Проверяем монотонность
        gradient_values = [f['gradient_value'] for f in features if f['type'] == 'transition']
        is_monotonic = all(x <= y for x, y in zip(gradient_values, gradient_values[1:]))
        logger.info(f"Transition zone gradient is {'monotonic' if is_monotonic else 'not monotonic'}")
        logger.info(f"Gradient values: {gradient_values}")

        return sequence, features


    def generate(self) -> Tuple[str, List[Dict]]:
        """Генерирует полный центромерный регион"""
        # Определяем размеры зон
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

        # Генерируем зоны
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
            self.centromere_params.transition_length,
            reverse_gradient=True
        )

        # Собираем полную последовательность
        sequence = (left_peri_seq + left_trans_seq +
                   core_seq + right_trans_seq + right_peri_seq)

        # Корректируем позиции features и собираем их вместе
        features = []
        current_pos = 0

        # Добавляем features с корректными позициями для каждой зоны
        for zone_features, zone_length in [
            (left_features, len(left_peri_seq)),
            (left_trans_features, len(left_trans_seq)),
            (core_features, len(core_seq)),
            (right_trans_features, len(right_trans_seq)),
            (right_features, len(right_peri_seq))
        ]:
            for feature in zone_features:
                feature = feature.copy()  # Создаем копию чтобы не модифицировать оригинал
                feature['start'] += current_pos
                feature['end'] += current_pos
                features.append(feature)
            current_pos += zone_length

        return sequence, features
