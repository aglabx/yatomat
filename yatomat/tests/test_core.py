import matplotlib.pyplot as plt
import numpy as np
from yatomat.core import SequenceGenerator, MutationEngine, GCProfiler, MutationParams
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def visualize_gc_profile(sequence: str, window_size: int = 100, title: str = "GC Profile"):
    """Визуализация GC-профиля последовательности"""
    gc_profiler = GCProfiler()
    gc_profile = gc_profiler.calculate_gc_content(sequence, window_size)
    
    plt.figure(figsize=(12, 6))
    plt.plot(gc_profile)
    plt.axhline(y=np.mean(gc_profile), color='r', linestyle='--', label=f'Mean GC: {np.mean(gc_profile):.3f}')
    plt.title(title)
    plt.xlabel("Position (window)")
    plt.ylabel("GC content")
    plt.legend()
    plt.grid(True)
    plt.show()

def test_sequence_generation():
    """Тест генерации последовательностей с разным GC-составом"""
    logger.info("Testing sequence generation with different GC content...")
    
    # Создаем последовательности с разным GC-составом
    gc_values = [0.3, 0.5, 0.7]
    sequence_length = 10000
    
    for gc in gc_values:
        seq_gen = SequenceGenerator(gc_content=gc, gc_std=0.05)
        sequence = seq_gen.generate_sequence(sequence_length)
        
        # Визуализируем GC-профиль
        visualize_gc_profile(sequence, title=f"GC Profile (target GC = {gc})")
        
        # Проверяем средний GC-состав
        actual_gc = (sequence.count('G') + sequence.count('C')) / len(sequence)
        logger.info(f"Target GC: {gc:.2f}, Actual GC: {actual_gc:.2f}")

def test_mutations():
    """Тест внесения мутаций"""
    logger.info("\nTesting mutation engine...")
    
    # Создаем исходную последовательность
    seq_gen = SequenceGenerator(gc_content=0.5)
    sequence = seq_gen.generate_sequence(5000)
    
    # Создаем движок мутаций с повышенными частотами для демонстрации
    mut_params = MutationParams(
        substitution_rate=0.02,
        insertion_rate=0.01,
        deletion_rate=0.01,
        duplication_rate=0.005,
        inversion_rate=0.005,
        min_event_size=1,
        max_event_size=20
    )
    mut_engine = MutationEngine(params=mut_params)
    
    # Вносим мутации
    mutated_seq, mutations = mut_engine.mutate_sequence(sequence)
    
    # Анализируем результаты
    mutation_types = {}
    for mut in mutations:
        mut_type = mut['type']
        mutation_types[mut_type] = mutation_types.get(mut_type, 0) + 1
    
    logger.info("Mutation statistics:")
    for mut_type, count in mutation_types.items():
        logger.info(f"{mut_type}: {count} events")
    
    # Визуализируем GC-профили до и после мутаций
    visualize_gc_profile(sequence, title="Original Sequence GC Profile")
    visualize_gc_profile(mutated_seq, title="Mutated Sequence GC Profile")

def test_gc_adjustment():
    """Тест корректировки GC-состава"""
    logger.info("\nTesting GC content adjustment...")
    
    # Создаем последовательность с низким GC
    seq_gen = SequenceGenerator(gc_content=0.3)
    sequence = seq_gen.generate_sequence(5000)
    
    # Корректируем до целевого значения
    target_gc = 0.6
    gc_profiler = GCProfiler()
    adjusted_seq = gc_profiler.adjust_gc_content(sequence, target_gc)
    
    # Визуализируем результаты
    visualize_gc_profile(sequence, title="Original Sequence (low GC)")
    visualize_gc_profile(adjusted_seq, title=f"Adjusted Sequence (target GC = {target_gc})")

if __name__ == "__main__":
    # Запускаем все тесты
    test_sequence_generation()
    test_mutations()
    test_gc_adjustment()