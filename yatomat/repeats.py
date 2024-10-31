from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import random
import numpy as np
from yatomat.core import SequenceGenerator, MutationEngine, GCProfiler, MutationParams
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RepeatType(Enum):
    """Типы тандемных повторов"""
    ALPHA_SATELLITE = "alpha_satellite"
    BETA_SATELLITE = "beta_satellite"
    SATELLITE_1 = "satellite_1"
    SATELLITE_2 = "satellite_2"
    SATELLITE_3 = "satellite_3"
    SIMPLE_REPEAT = "simple_repeat"
    TELOMERIC = "telomeric"

@dataclass
class RepeatUnit:
    """Базовая единица повтора"""
    sequence: str
    type: RepeatType
    size: int
    gc_content: float
    metadata: Dict = None

@dataclass
class HORParams:
    """Параметры для создания HOR структуры"""
    unit_size: int
    units_in_hor: int
    hor_copies: int
    mutation_rate_between_units: float = 0.05
    mutation_rate_between_hors: float = 0.02
    gc_content: float = 0.5
    gc_std: float = 0.02

class RepeatGenerator:
    """Генератор различных типов повторов"""
    
    def __init__(self):
        self.seq_gen = SequenceGenerator()
        self.mut_engine = MutationEngine()
        self.gc_profiler = GCProfiler()
        
        # Определяем стандартные характеристики известных повторов
        self.repeat_characteristics = {
            RepeatType.SIMPLE_REPEAT: {
                'unit_size': 2,  # Может варьировать от 1 до 6
                'gc_content': 0.5,
                'typical_hor_units': 1
            },
            RepeatType.ALPHA_SATELLITE: {
                'unit_size': 171,
                'gc_content': 0.58,
                'typical_hor_units': 12
            },
            RepeatType.BETA_SATELLITE: {
                'unit_size': 68,
                'gc_content': 0.54,
                'typical_hor_units': 4
            },
            RepeatType.SATELLITE_1: {
                'unit_size': 42,
                'gc_content': 0.48,
                'typical_hor_units': 1
            },
            RepeatType.SATELLITE_2: {
                'unit_size': 5,
                'gc_content': 0.65,
                'typical_hor_units': 1
            },
            RepeatType.SATELLITE_3: {
                'unit_size': 5,
                'gc_content': 0.48,
                'typical_hor_units': 1
            },
            RepeatType.TELOMERIC: {
                'unit_size': 6,
                'gc_content': 0.33,
                'typical_hor_units': 1,
                'consensus': 'TTAGGG'
            }
        }
    
    def generate_monomer(self, repeat_type: RepeatType, 
                        custom_params: Optional[Dict] = None) -> RepeatUnit:
        """
        Генерирует мономер заданного типа повтора
        
        Args:
            repeat_type: тип повтора
            custom_params: пользовательские параметры (переопределяют стандартные)
        """
        # Получаем базовые характеристики
        if repeat_type not in self.repeat_characteristics:
            raise ValueError(f"Неизвестный тип повтора: {repeat_type}")
        
        params = self.repeat_characteristics[repeat_type].copy()
        if custom_params:
            params.update(custom_params)
        
        # Для теломерного повтора используем консенсусную последовательность
        if repeat_type == RepeatType.TELOMERIC:
            sequence = params['consensus']
        elif repeat_type == RepeatType.SIMPLE_REPEAT:
            # Для простых повторов генерируем короткий паттерн
            unit_size = custom_params.get('unit_size', random.randint(1, 6))
            pattern = self.seq_gen.generate_sequence(unit_size, local_gc=params['gc_content'])
            sequence = pattern
        else:
            sequence = self.seq_gen.generate_sequence(
                params['unit_size'],
                local_gc=params['gc_content']
            )
        
        return RepeatUnit(
            sequence=sequence,
            type=repeat_type,
            size=params['unit_size'],
            gc_content=params['gc_content'],
            metadata=params
        )

class HORGenerator:
    """Генератор HOR структур"""
    
    def __init__(self, params: Optional[HORParams] = None):
        self.params = params or HORParams(
            unit_size=171,  # размер альфа-сателлита по умолчанию
            units_in_hor=12,
            hor_copies=10,
            mutation_rate_between_units=0.05,
            mutation_rate_between_hors=0.02,
            gc_content=0.58
        )
        self.mut_engine = MutationEngine()
        self.repeat_gen = RepeatGenerator()
    
    def generate_hor(self, base_repeat: Optional[RepeatUnit] = None) -> Tuple[str, List[Dict]]:
        """
        Генерирует одну HOR единицу
        
        Args:
            base_repeat: базовый повторяющийся элемент (если None, генерируется новый)
        """
        if base_repeat is None:
            base_repeat = self.repeat_gen.generate_monomer(
                RepeatType.ALPHA_SATELLITE,
                {'unit_size': self.params.unit_size, 'gc_content': self.params.gc_content}
            )
        
        mutations_log = []
        hor_sequence = ""
        
        # Создаем варианты мономеров для HOR
        for i in range(self.params.units_in_hor):
            mut_params = MutationParams(
                substitution_rate=self.params.mutation_rate_between_units,
                insertion_rate=self.params.mutation_rate_between_units/5,
                deletion_rate=self.params.mutation_rate_between_units/5
            )
            self.mut_engine.params = mut_params
            
            if i == 0:
                current_monomer = base_repeat.sequence
            else:
                current_monomer, muts = self.mut_engine.mutate_sequence(base_repeat.sequence)
                mutations_log.extend([{**m, 'unit_number': i} for m in muts])
            
            hor_sequence += current_monomer
        
        return hor_sequence, mutations_log

    def generate_array(self) -> Tuple[str, List[Dict]]:
        """Генерирует массив HOR повторов"""
        base_hor, base_mutations = self.generate_hor()
        array_sequence = ""
        mutations_log = base_mutations
        
        # Создаем копии HOR с вариациями
        for i in range(self.params.hor_copies):
            mut_params = MutationParams(
                substitution_rate=self.params.mutation_rate_between_hors,
                insertion_rate=self.params.mutation_rate_between_hors/5,
                deletion_rate=self.params.mutation_rate_between_hors/5
            )
            self.mut_engine.params = mut_params
            
            if i == 0:
                current_hor = base_hor
            else:
                current_hor, muts = self.mut_engine.mutate_sequence(base_hor)
                mutations_log.extend([{**m, 'hor_number': i} for m in muts])
            
            array_sequence += current_hor
        
        return array_sequence, mutations_log

class HomogenizationEngine:
    """Движок для моделирования процессов гомогенизации"""
    
    def __init__(self, conversion_rate: float = 0.1, 
                 min_tract_length: int = 50,
                 max_tract_length: int = 200):
        self.conversion_rate = conversion_rate
        self.min_tract_length = min_tract_length
        self.max_tract_length = max_tract_length
    
    def apply_gene_conversion(self, sequence: str, 
                            unit_size: int,
                            units_in_hor: int) -> Tuple[str, List[Dict]]:
        """
        Применяет gene conversion-подобные события к последовательности
        
        Args:
            sequence: входная последовательность
            unit_size: размер базового повтора
            units_in_hor: количество повторов в HOR
        """
        seq_length = len(sequence)
        hor_size = unit_size * units_in_hor
        conversions = []
        
        # Определяем количество событий конверсии
        n_events = np.random.binomial(seq_length // hor_size, self.conversion_rate)
        
        sequence = list(sequence)
        for _ in range(n_events):
            # Выбираем случайный участок для конверсии
            tract_length = random.randint(self.min_tract_length, 
                                        min(self.max_tract_length, seq_length))
            start = random.randint(0, seq_length - tract_length)
            
            # Выбираем донорный участок
            donor_start = random.randint(0, seq_length - tract_length)
            while abs(donor_start - start) < tract_length:
                donor_start = random.randint(0, seq_length - tract_length)
            
            # Применяем конверсию
            donor_tract = sequence[donor_start:donor_start + tract_length]
            sequence[start:start + tract_length] = donor_tract
            
            conversions.append({
                'type': 'gene_conversion',
                'start': start,
                'length': tract_length,
                'donor_start': donor_start
            })
        
        return ''.join(sequence), conversions

if __name__ == "__main__":
    # Пример использования
    repeat_gen = RepeatGenerator()
    
    # Генерируем теломерный повтор
    telomere = repeat_gen.generate_monomer(RepeatType.TELOMERIC)
    print(f"Telomeric repeat: {telomere.sequence}")
    
    # Генерируем HOR структуру
    hor_params = HORParams(
        unit_size=171,
        units_in_hor=12,
        hor_copies=10,
        mutation_rate_between_units=0.05,
        mutation_rate_between_hors=0.02
    )
    
    hor_gen = HORGenerator(hor_params)
    array_sequence, mutations = hor_gen.generate_array()
    
    print(f"\nGenerated HOR array length: {len(array_sequence)}")
    print(f"Number of mutations: {len(mutations)}")
    
    # Применяем гомогенизацию
    homogenization = HomogenizationEngine()
    homogenized_seq, conversions = homogenization.apply_gene_conversion(
        array_sequence, 
        unit_size=171,
        units_in_hor=12
    )
    
    print(f"\nNumber of gene conversion events: {len(conversions)}")