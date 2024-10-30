import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from typing import List
import logging

# Добавляем путь к модулям
sys.path.append(str(Path(__file__).parent))

from regions.common import (
    ChromosomeRegionType, RegionParams, GradientParams,
    GradientGenerator, ChromosomeRegion, SequenceFeature,
    RegionBuilder
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TestRegion(ChromosomeRegion):
    """Тестовый класс региона для проверки абстрактного класса"""
    def generate(self):
        sequence = "ATCG" * (self.length // 4)
        if self.length % 4:
            sequence += "ATCG"[:self.length % 4]
        self.sequence = sequence
        return sequence, []

class TestCommonFunctionality(unittest.TestCase):
    """Тесты для базовой функциональности модуля common"""
    
    def setUp(self):
        """Установка начальных данных для тестов"""
        self.region_params = RegionParams(
            region_type=ChromosomeRegionType.CENTROMERE,
            start=1000,
            end=2000,
            gc_content=0.5,
            gc_std=0.02,
            repeat_density=0.7
        )
    
    def test_region_params_validation(self):
        """Тест валидации параметров региона"""
        # Тест некорректных границ
        with self.assertRaises(ValueError):
            RegionParams(
                region_type=ChromosomeRegionType.CENTROMERE,
                start=1000,
                end=500
            )
        
        # Тест некорректного GC-состава
        with self.assertRaises(ValueError):
            RegionParams(
                region_type=ChromosomeRegionType.CENTROMERE,
                start=0,
                end=1000,
                gc_content=1.5
            )
        
        # Тест некорректной плотности повторов
        with self.assertRaises(ValueError):
            RegionParams(
                region_type=ChromosomeRegionType.CENTROMERE,
                start=0,
                end=1000,
                repeat_density=2.0
            )
    
    def test_region_builder(self):
        """Тест построителя параметров региона"""
        builder = RegionBuilder(ChromosomeRegionType.TELOMERE)
        params = (builder
                 .set_boundaries(0, 1000)
                 .set_gc_content(0.45)
                 .set_repeat_density(0.9)
                 .add_metadata('orientation', 'forward')
                 .build())
        
        self.assertEqual(params.start, 0)
        self.assertEqual(params.end, 1000)
        self.assertEqual(params.gc_content, 0.45)
        self.assertEqual(params.repeat_density, 0.9)
        self.assertEqual(params.metadata['orientation'], 'forward')
    
    def test_gradient_generator(self):
        """Тест генератора градиентов"""
        length = 100
        gradient_types = ['linear', 'exponential', 'sigmoid']
        gradients = []
        
        plt.figure(figsize=(12, 8))
        for shape in gradient_types:
            params = GradientParams(
                start_value=0.2,
                end_value=0.8,
                shape=shape,
                steepness=5 if shape in ['exponential', 'sigmoid'] else 1
            )
            gradient = GradientGenerator.create_gradient(params, length)
            gradients.append(gradient)
            
            # Проверяем граничные значения
            self.assertAlmostEqual(gradient[0], 0.2, places=2)
            self.assertAlmostEqual(gradient[-1], 0.8, places=2)
            
            # Проверяем монотонность
            self.assertTrue(np.all(np.diff(gradient) >= 0))
            
            # Визуализация
            plt.plot(gradient, label=shape)
        
        plt.title('Comparison of Different Gradient Types')
        plt.xlabel('Position')
        plt.ylabel('Value')
        plt.legend()
        plt.grid(True)
        plt.show()
    
    def test_sequence_feature(self):
        """Тест работы с особенностями последовательности"""
        feature = SequenceFeature(
            feature_type='repeat',
            sequence='ATCGATCG',
            position=100,
            metadata={'type': 'tandem', 'copies': 2}
        )
        
        feature_dict = feature.to_dict()
        self.assertEqual(feature_dict['start'], 100)
        self.assertEqual(feature_dict['end'], 108)
        self.assertEqual(feature_dict['metadata']['copies'], 2)
    
    def test_chromosome_region(self):
        """Тест базового класса региона"""
        region = TestRegion(self.region_params)
        
        # Тест генерации последовательности
        sequence, _ = region.generate()
        self.assertEqual(len(sequence), region.length)
        
        # Тест добавления аннотаций
        region.add_annotation('repeat', 50, 100, {'type': 'satellite'})
        self.assertEqual(len(region.annotations), 1)
        self.assertEqual(region.annotations[0]['start'], 1050)  # start + 50
        
        # Тест подсчета GC-состава
        gc_content = region.get_gc_content()
        self.assertEqual(gc_content, 0.5)  # для нашей тестовой последовательности ATCG

def visualize_gradients():
    """Визуализация различных типов градиентов"""
    length = 100
    gradient_types = {
        'linear': {'steepness': 1},
        'exponential': {'steepness': 3},
        'sigmoid': {'steepness': 8}
    }
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Gradient Visualization')
    
    # Простые градиенты
    ax = axes[0, 0]
    for shape, params in gradient_types.items():
        gradient_params = GradientParams(
            start_value=0.2,
            end_value=0.8,
            shape=shape,
            steepness=params['steepness']
        )
        gradient = GradientGenerator.create_gradient(gradient_params, length)
        ax.plot(gradient, label=shape)
    
    ax.set_title('Simple Gradients')
    ax.set_xlabel('Position')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True)
    
    # Сигмоидальные градиенты с разной крутизной
    ax = axes[0, 1]
    steepness_values = [2, 5, 10]
    for steepness in steepness_values:
        params = GradientParams(
            start_value=0.2,
            end_value=0.8,
            shape='sigmoid',
            steepness=steepness
        )
        gradient = GradientGenerator.create_gradient(params, length)
        ax.plot(gradient, label=f'Steepness {steepness}')
    
    ax.set_title('Sigmoid Gradients with Different Steepness')
    ax.set_xlabel('Position')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True)
    
    # Комбинированный градиент (GC-состав)
    ax = axes[1, 0]
    gc_params = GradientParams(
        start_value=0.3,
        end_value=0.7,
        shape='sigmoid',
        steepness=5
    )
    gradient = GradientGenerator.create_gradient(gc_params, length)
    noise = np.random.normal(0, 0.02, length)
    noisy_gradient = np.clip(gradient + noise, 0, 1)
    
    ax.plot(gradient, label='Smooth')
    ax.plot(noisy_gradient, label='With noise', alpha=0.7)
    ax.set_title('GC Content Gradient with Noise')
    ax.set_xlabel('Position')
    ax.set_ylabel('GC Content')
    ax.legend()
    ax.grid(True)
    
    # Двунаправленный градиент
    ax = axes[1, 1]
    params1 = GradientParams(start_value=0.5, end_value=0.8, shape='sigmoid', steepness=5)
    params2 = GradientParams(start_value=0.5, end_value=0.2, shape='sigmoid', steepness=5)
    
    gradient1 = GradientGenerator.create_gradient(params1, length//2)
    gradient2 = GradientGenerator.create_gradient(params2, length//2)
    combined = np.concatenate([gradient1, gradient2])
    
    ax.plot(combined)
    ax.set_title('Bidirectional Gradient')
    ax.set_xlabel('Position')
    ax.set_ylabel('Value')
    ax.grid(True)
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # Запускаем визуализацию
    visualize_gradients()
    
    # Запускаем тесты
    unittest.main(argv=['first-arg-is-ignored'], exit=False)