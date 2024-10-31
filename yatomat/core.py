from typing import List, Tuple, Dict, Optional, Union
import random
import numpy as np
from dataclasses import dataclass
from enum import Enum
import logging

# Настройка логирования
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MutationType(Enum):
    """Типы возможных мутаций"""
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"
    DUPLICATION = "duplication"
    INVERSION = "inversion"

@dataclass
class MutationParams:
    """Параметры для внесения мутаций"""
    substitution_rate: float = 0.01
    insertion_rate: float = 0.005
    deletion_rate: float = 0.005
    duplication_rate: float = 0.002
    inversion_rate: float = 0.001
    min_event_size: int = 1
    max_event_size: int = 10

class SequenceGenerator:
    """Базовый класс для генерации последовательностей"""
    
    def __init__(self, gc_content: float = 0.5, gc_std: float = 0.02):
        """
        Args:
            gc_content: целевой GC-состав (0-1)
            gc_std: стандартное отклонение GC-состава
        """
        if not 0 <= gc_content <= 1:
            raise ValueError("GC-состав должен быть между 0 и 1")
        self.gc_content = gc_content
        self.gc_std = gc_std
        self._validate_params()
    
    def _validate_params(self):
        """Валидация параметров"""
        if self.gc_std < 0:
            raise ValueError("Стандартное отклонение не может быть отрицательным")
        if self.gc_content - 3*self.gc_std < 0 or self.gc_content + 3*self.gc_std > 1:
            logger.warning("GC-состав может выйти за пределы 0-1 с текущими параметрами")
    
    def _get_base_probabilities(self) -> Dict[str, float]:
        """Рассчитывает вероятности для каждого нуклеотида"""
        gc_prob = self.gc_content / 2
        at_prob = (1 - self.gc_content) / 2
        return {
            'G': gc_prob,
            'C': gc_prob,
            'A': at_prob,
            'T': at_prob
        }
    
    def generate_sequence(self, length: int, local_gc: Optional[float] = None) -> str:
        """
        Генерирует последовательность заданной длины
        
        Args:
            length: длина последовательности
            local_gc: локальный GC-состав (если None, используется self.gc_content)
        """
        if local_gc is None:
            local_gc = np.random.normal(self.gc_content, self.gc_std)
            local_gc = np.clip(local_gc, 0, 1)
        
        gc_prob = local_gc / 2
        at_prob = (1 - local_gc) / 2
        
        nucleotides = ['G', 'C', 'A', 'T']
        probabilities = [gc_prob, gc_prob, at_prob, at_prob]
        
        return ''.join(np.random.choice(nucleotides, size=length, p=probabilities))

class MutationEngine:
    """Класс для внесения мутаций в последовательности"""
    
    def __init__(self, params: Optional[MutationParams] = None):
        self.params = params or MutationParams()
        self.seq_gen = SequenceGenerator()
    
    def _apply_substitution(self, sequence: str, position: int) -> str:
        """Замена одного нуклеотида"""
        bases = ['A', 'T', 'G', 'C']
        current_base = sequence[position]
        bases.remove(current_base)
        new_base = random.choice(bases)
        return sequence[:position] + new_base + sequence[position + 1:]
    
    def _apply_insertion(self, sequence: str, position: int, size: int) -> str:
        """Вставка последовательности"""
        insert_seq = self.seq_gen.generate_sequence(size)
        return sequence[:position] + insert_seq + sequence[position:]
    
    def _apply_deletion(self, sequence: str, position: int, size: int) -> str:
        """Удаление участка последовательности"""
        return sequence[:position] + sequence[position + size:]
    
    def _apply_duplication(self, sequence: str, position: int, size: int) -> str:
        """Тандемная дупликация участка"""
        dup_seq = sequence[position:position + size]
        return sequence[:position] + dup_seq + sequence[position:]
    
    def _apply_inversion(self, sequence: str, position: int, size: int) -> str:
        """Инверсия участка"""
        inv_seq = sequence[position:position + size][::-1]
        return sequence[:position] + inv_seq + sequence[position + size:]
    
    def mutate_sequence(self, sequence: str) -> Tuple[str, List[Dict]]:
        """
        Вносит мутации в последовательность согласно заданным параметрам
        
        Returns:
            Tuple[str, List[Dict]]: (мутированная последовательность, список мутаций)
        """
        mutations = []
        seq_length = len(sequence)
        sequence = list(sequence)
        
        # Проходим по каждому типу мутаций
        for mut_type in MutationType:
            rate = getattr(self.params, f"{mut_type.value}_rate")
            
            # Определяем количество мутаций данного типа
            n_mutations = np.random.binomial(seq_length, rate)
            
            for _ in range(n_mutations):
                position = random.randint(0, seq_length - 1)
                size = random.randint(self.params.min_event_size, 
                                    min(self.params.max_event_size, seq_length - position))
                
                # Применяем мутацию
                if mut_type == MutationType.SUBSTITUTION:
                    sequence = self._apply_substitution(''.join(sequence), position)
                elif mut_type == MutationType.INSERTION:
                    sequence = self._apply_insertion(''.join(sequence), position, size)
                elif mut_type == MutationType.DELETION:
                    sequence = self._apply_deletion(''.join(sequence), position, size)
                elif mut_type == MutationType.DUPLICATION:
                    sequence = self._apply_duplication(''.join(sequence), position, size)
                elif mut_type == MutationType.INVERSION:
                    sequence = self._apply_inversion(''.join(sequence), position, size)
                
                mutations.append({
                    'type': mut_type.value,
                    'position': position,
                    'size': size
                })
                
                # Обновляем длину последовательности
                seq_length = len(sequence)
        
        return ''.join(sequence), mutations

class GCProfiler:
    """Класс для анализа и модификации GC-состава"""
    
    @staticmethod
    def calculate_gc_content(sequence: str, window_size: int = 100) -> List[float]:
        """Рассчитывает GC-состав в скользящем окне"""
        if window_size > len(sequence):
            window_size = len(sequence)
        
        gc_profile = []
        for i in range(0, len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            gc_count = window.count('G') + window.count('C')
            gc_profile.append(gc_count / window_size)
        
        return gc_profile
    
    @staticmethod
    def adjust_gc_content(sequence: str, target_gc: float, 
                         window_size: int = 100) -> str:
        """
        Корректирует GC-состав последовательности до целевого значения
        
        Args:
            sequence: исходная последовательность
            target_gc: целевой GC-состав
            window_size: размер окна для анализа
        """
        seq_list = list(sequence)
        current_gc = (sequence.count('G') + sequence.count('C')) / len(sequence)
        
        while abs(current_gc - target_gc) > 0.01:
            # Выбираем случайную позицию
            pos = random.randint(0, len(sequence) - 1)
            current_base = seq_list[pos]
            
            if current_gc < target_gc:
                # Нужно увеличить GC-состав
                if current_base in ['A', 'T']:
                    seq_list[pos] = random.choice(['G', 'C'])
            else:
                # Нужно уменьшить GC-состав
                if current_base in ['G', 'C']:
                    seq_list[pos] = random.choice(['A', 'T'])
            
            current_gc = (seq_list.count('G') + seq_list.count('C')) / len(seq_list)
        
        return ''.join(seq_list)

if __name__ == "__main__":
    # Пример использования
    seq_gen = SequenceGenerator(gc_content=0.45, gc_std=0.05)
    sequence = seq_gen.generate_sequence(1000)
    
    mut_engine = MutationEngine()
    mutated_seq, mutations = mut_engine.mutate_sequence(sequence)
    
    gc_profiler = GCProfiler()
    gc_profile = gc_profiler.calculate_gc_content(mutated_seq)
    
    print(f"Original sequence length: {len(sequence)}")
    print(f"Mutated sequence length: {len(mutated_seq)}")
    print(f"Number of mutations: {len(mutations)}")
    print(f"Average GC content: {sum(gc_profile)/len(gc_profile):.3f}")