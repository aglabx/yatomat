from typing import List, Dict, Tuple, Optional
import numpy as np
from pathlib import Path
import logging
from dataclasses import dataclass

from .common import (
    ChromosomeRegion, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator, SequenceFeature
)
from repeats import RepeatGenerator, RepeatType, HORGenerator, HomogenizationEngine

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class TelomereParams:
    """Параметры для генерации теломерного региона"""
    min_length: int = 3000  # Минимальная длина теломеры (3 kb)
    max_length: int = 15000  # Максимальная длина теломеры (15 kb)
    consensus_repeat: str = 'TTAGGG'  # Консенсусная последовательность
    mutation_rate_5_end: float = 0.01  # Частота мутаций на 5' конце
    mutation_rate_3_end: float = 0.05  # Частота мутаций на 3' конце
    variant_repeat_prob: float = 0.1  # Вероятность вариантных повторов
    variant_types: List[str] = None  # Типы вариантных повторов

    def __post_init__(self):
        if self.variant_types is None:
            # Известные варианты теломерных повторов
            self.variant_types = [
                'TCAGGG',  # C-тип
                'TGAGGG',  # G-тип
                'TTGGGG',  # GG-вариант
                'TTAGAG',  # А-вариант
                'TTAGA',   # Укороченный вариант
            ]

@dataclass
class SubtelomereParams:
    """Параметры для генерации субтеломерного региона"""
    min_length: int = 10000  # Минимальная длина (10 kb)
    max_length: int = 300000  # Максимальная длина (300 kb)
    satellite_density: float = 0.4  # Плотность сателлитных повторов
    transposon_density: float = 0.2  # Плотность транспозонов
    segment_dup_prob: float = 0.3  # Вероятность сегментных дупликаций

class TelomereRegion(ChromosomeRegion):
    """Класс для генерации теломерного региона"""
    
    def __init__(self, params: RegionParams, telomere_params: Optional[TelomereParams] = None):
        super().__init__(params)
        self.telomere_params = telomere_params or TelomereParams()
        self.repeat_gen = RepeatGenerator()
        self.gradient_gen = GradientGenerator()
    
    def _generate_variant_repeat(self) -> str:
        """Генерирует вариантный повтор"""
        return np.random.choice(self.telomere_params.variant_types)
    
    def _apply_mutations(self, sequence: str, mutation_rate: float) -> str:
        """Вносит мутации в последовательность"""
        sequence = list(sequence)
        for i in range(len(sequence)):
            if np.random.random() < mutation_rate:
                # 80% замен, 10% инсерций, 10% делеций
                mutation_type = np.random.choice(['sub', 'ins', 'del'], p=[0.8, 0.1, 0.1])
                
                if mutation_type == 'sub':
                    nucleotides = ['A', 'T', 'G', 'C']
                    nucleotides.remove(sequence[i])
                    sequence[i] = np.random.choice(nucleotides)
                elif mutation_type == 'ins':
                    sequence.insert(i, np.random.choice(['A', 'T', 'G', 'C']))
                elif mutation_type == 'del' and len(sequence) > 1:
                    sequence.pop(i)
                    break
        
        return ''.join(sequence)
    
    def generate(self) -> Tuple[str, List[Dict]]:
        """Генерирует теломерный регион"""
        # Определяем длину теломеры
        telomere_length = np.random.randint(
            self.telomere_params.min_length,
            self.telomere_params.max_length
        )
        
        # Создаем градиент мутаций
        mutation_gradient = self.gradient_gen.create_gradient(
            GradientParams(
                start_value=self.telomere_params.mutation_rate_5_end,
                end_value=self.telomere_params.mutation_rate_3_end,
                shape='exponential',
                steepness=3
            ),
            telomere_length // len(self.telomere_params.consensus_repeat)
        )
        
        # Генерируем последовательность
        sequence = ""
        features = []
        current_pos = 0
        
        while current_pos < telomere_length:
            # Определяем, будет ли это вариантный повтор
            if np.random.random() < self.telomere_params.variant_repeat_prob:
                repeat = self._generate_variant_repeat()
                repeat_type = 'variant'
            else:
                repeat = self.telomere_params.consensus_repeat
                repeat_type = 'consensus'
            
            # Вносим мутации согласно градиенту
            gradient_idx = current_pos // len(self.telomere_params.consensus_repeat)
            if gradient_idx < len(mutation_gradient):
                repeat = self._apply_mutations(repeat, mutation_gradient[gradient_idx])
            
            # Добавляем повтор и его аннотацию
            sequence += repeat
            features.append({
                'type': repeat_type,
                'start': current_pos,
                'end': current_pos + len(repeat),
                'sequence': repeat
            })
            
            current_pos += len(repeat)
        
        self.sequence = sequence[:telomere_length]
        return self.sequence, features

class SubtelomereRegion(ChromosomeRegion):
    """Класс для генерации субтеломерного региона"""
    
    def __init__(self, params: RegionParams, 
                 subtelomere_params: Optional[SubtelomereParams] = None):
        super().__init__(params)
        self.subtelomere_params = subtelomere_params or SubtelomereParams()
        self.repeat_gen = RepeatGenerator()
        self.hor_gen = HORGenerator()
    
    def _generate_satellite_block(self, length: int) -> Tuple[str, Dict]:
        """Генерирует блок сателлитных повторов"""
        # Выбираем тип сателлита
        satellite_type = np.random.choice([
            RepeatType.SATELLITE_1,
            RepeatType.SATELLITE_2,
            RepeatType.SATELLITE_3
        ])
        
        # Генерируем базовый мономер
        monomer = self.repeat_gen.generate_monomer(satellite_type)
        
        # Создаем массив повторов
        copies = length // len(monomer.sequence)
        sequence = monomer.sequence * copies
        
        return sequence, {
            'type': 'satellite',
            'satellite_type': satellite_type.value,
            'monomer_size': len(monomer.sequence),
            'copies': copies
        }
    
    def _generate_transposon(self, min_length: int = 1000, 
                           max_length: int = 5000) -> Tuple[str, Dict]:
        """Генерирует транспозон-подобный элемент"""
        length = np.random.randint(min_length, max_length)
        sequence = ''.join(np.random.choice(['A', 'T', 'G', 'C']) for _ in range(length))
        
        return sequence, {
            'type': 'transposon',
            'length': length
        }
    
    def generate(self) -> Tuple[str, List[Dict]]:
        """Генерирует субтеломерный регион"""
        subtelomere_length = np.random.randint(
            self.subtelomere_params.min_length,
            self.subtelomere_params.max_length
        )
        
        sequence = ""
        features = []
        current_pos = 0
        
        while current_pos < subtelomere_length:
            # Определяем тип следующего элемента
            if np.random.random() < self.subtelomere_params.satellite_density:
                # Генерируем блок сателлитов
                block_size = np.random.randint(1000, 5000)
                block_seq, block_info = self._generate_satellite_block(block_size)
            else:
                # Генерируем транспозон
                block_seq, block_info = self._generate_transposon()
            
            # Добавляем последовательность и аннотацию
            sequence += block_seq
            block_info.update({
                'start': current_pos,
                'end': current_pos + len(block_seq)
            })
            features.append(block_info)
            
            current_pos += len(block_seq)
        
        # Обрезаем до нужной длины
        self.sequence = sequence[:subtelomere_length]
        
        # Корректируем координаты последнего элемента
        if features:
            features[-1]['end'] = min(features[-1]['end'], subtelomere_length)
        
        return self.sequence, features

if __name__ == "__main__":
    # Пример использования
    from common import RegionBuilder
    
    # Создаем теломерный регион
    telomere_params = RegionBuilder(ChromosomeRegionType.TELOMERE)\
        .set_boundaries(0, 10000)\
        .set_gc_content(0.48)\
        .build()
    
    telomere = TelomereRegion(telomere_params)
    sequence, features = telomere.generate()
    
    logger.info(f"Generated telomere sequence of length {len(sequence)}")
    logger.info(f"Number of features: {len(features)}")
    
    # Создаем субтеломерный регион
    subtelomere_params = RegionBuilder(ChromosomeRegionType.SUBTELOMERE)\
        .set_boundaries(10000, 50000)\
        .set_gc_content(0.52)\
        .build()
    
    subtelomere = SubtelomereRegion(subtelomere_params)
    sequence, features = subtelomere.generate()
    
    logger.info(f"Generated subtelomere sequence of length {len(sequence)}")
    logger.info(f"Number of features: {len(features)}")