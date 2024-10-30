from typing import List, Dict, Optional, Union, Tuple, Any
from dataclasses import dataclass
from enum import Enum, auto
import numpy as np
from abc import ABC, abstractmethod
import logging
from pathlib import Path

# Убедимся, что директория regions существует
Path("regions").mkdir(exist_ok=True)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ChromosomeRegionType(Enum):
    """Типы хромосомных регионов"""
    TELOMERE = auto()
    SUBTELOMERE = auto()
    CENTROMERE = auto()
    PERICENTROMERE = auto()
    TRANSITION = auto()
    TAD = auto()

@dataclass
class RegionParams:
    """Базовые параметры для всех регионов"""
    region_type: ChromosomeRegionType
    start: int
    end: int
    gc_content: float = 0.5
    gc_std: float = 0.02
    repeat_density: float = 0.0  # доля повторов в регионе
    metadata: Dict = None

    def __post_init__(self):
        if self.end <= self.start:
            raise ValueError("end должен быть больше start")
        if not 0 <= self.gc_content <= 1:
            raise ValueError("gc_content должен быть между 0 и 1")
        if not 0 <= self.repeat_density <= 1:
            raise ValueError("repeat_density должен быть между 0 и 1")
        self.metadata = self.metadata or {}

@dataclass
class GradientParams:
    """Параметры для создания градиентов"""
    start_value: float
    end_value: float
    shape: str = 'linear'  # linear, exponential, sigmoid
    steepness: float = 1.0  # для exponential и sigmoid
    
    def __post_init__(self):
        if self.shape not in ['linear', 'exponential', 'sigmoid']:
            raise ValueError(f"Неподдерживаемый тип градиента: {self.shape}")

class GradientGenerator:
    """Генератор градиентов для различных параметров"""
    
    @staticmethod
    def create_gradient(params: GradientParams, length: int) -> np.ndarray:
        """
        Создает градиент заданной формы
        
        Args:
            params: параметры градиента
            length: длина градиента
        
        Returns:
            np.ndarray: массив значений градиента
        """
        x = np.linspace(0, 1, length)
        
        if params.shape == 'linear':
            values = x

        elif params.shape == 'exponential':
            values = (np.exp(params.steepness * x) - 1) / (np.exp(params.steepness) - 1)

        elif params.shape == 'sigmoid':
            values = 1 / (1 + np.exp(-params.steepness * (x - 0.5)))
            # Нормализуем сигмоид к диапазону [0, 1]
            values = (values - values[0]) / (values[-1] - values[0])
        
        # Масштабируем значения к целевому диапазону
        values = values * (params.end_value - params.start_value) + params.start_value
        return values

class ChromosomeRegion(ABC):
    """Абстрактный базовый класс для всех хромосомных регионов"""
    
    def __init__(self, params: RegionParams):
        """
        Args:
            params: параметры региона
        """
        self.params = params
        self.sequence: str = ""
        self.annotations: List[Dict] = []
    
    @abstractmethod
    def generate(self) -> Tuple[str, List[Dict]]:
        """
        Генерирует последовательность региона и его аннотации
        
        Returns:
            Tuple[str, List[Dict]]: (последовательность, аннотации)
        """
        pass
    
    @property
    def length(self) -> int:
        """Длина региона"""
        return self.params.end - self.params.start
    
    def add_annotation(self, feature_type: str, start: int, end: int, 
                      metadata: Optional[Dict] = None):
        """
        Добавляет аннотацию для участка последовательности
        
        Args:
            feature_type: тип особенности (повтор, домен и т.д.)
            start: начальная позиция
            end: конечная позиция
            metadata: дополнительная информация
        """
        annotation = {
            'type': feature_type,
            'start': self.params.start + start,
            'end': self.params.start + end,
            'metadata': metadata or {}
        }
        self.annotations.append(annotation)
    
    def get_gc_content(self, start: int = None, end: int = None) -> float:
        """
        Вычисляет GC-состав для участка последовательности
        
        Args:
            start: начальная позиция (если None, то с начала)
            end: конечная позиция (если None, то до конца)
        """
        if not self.sequence:
            return 0.0
        
        seq = self.sequence[start:end]
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq) if seq else 0.0

class SequenceFeature:
    """Класс для работы с особенностями последовательности"""
    
    def __init__(self, feature_type: str, sequence: str, 
                 position: int, metadata: Optional[Dict] = None):
        """
        Args:
            feature_type: тип особенности
            sequence: последовательность
            position: позиция в регионе
            metadata: дополнительная информация
        """
        self.type = feature_type
        self.sequence = sequence
        self.position = position
        self.metadata = metadata or {}
        self.length = len(sequence)
    
    def to_dict(self) -> Dict:
        """Преобразует особенность в словарь для аннотации"""
        return {
            'type': self.type,
            'start': self.position,
            'end': self.position + self.length,
            'metadata': self.metadata
        }

class RegionBuilder:
    """Вспомогательный класс для построения регионов"""
    
    def __init__(self, region_type: ChromosomeRegionType):
        self.params = RegionParams(
            region_type=region_type,
            start=0,
            end=1  # Устанавливаем минимальную валидную длину
        )
    
    def set_boundaries(self, start: int, end: int) -> 'RegionBuilder':
        """Установка границ региона с валидацией"""
        if end <= start:
            raise ValueError("end должен быть больше start")
        self.params.start = start
        self.params.end = end
        return self
    
    def set_gc_content(self, gc_content: float, gc_std: float = 0.02) -> 'RegionBuilder':
        self.params.gc_content = gc_content
        self.params.gc_std = gc_std
        return self
    
    def set_repeat_density(self, density: float) -> 'RegionBuilder':
        self.params.repeat_density = density
        return self
    
    def add_metadata(self, key: str, value: Any) -> 'RegionBuilder':
        if self.params.metadata is None:
            self.params.metadata = {}
        self.params.metadata[key] = value
        return self
    
    def build(self) -> RegionParams:
        return self.params

def create_transition_gradient(region1: ChromosomeRegion, 
                             region2: ChromosomeRegion,
                             transition_length: int,
                             params: Optional[GradientParams] = None) -> np.ndarray:
    """
    Создает градиент для переходной зоны между регионами
    
    Args:
        region1: первый регион
        region2: второй регион
        transition_length: длина переходной зоны
        params: параметры градиента (если None, используется линейный)
    """
    if params is None:
        params = GradientParams(
            start_value=region1.params.gc_content,
            end_value=region2.params.gc_content,
            shape='linear'
        )
    
    gradient_gen = GradientGenerator()
    return gradient_gen.create_gradient(params, transition_length)

# Примеры использования
if __name__ == "__main__":
    # Создание параметров региона через builder
    params = (RegionBuilder(ChromosomeRegionType.CENTROMERE)
             .set_boundaries(1000, 5000)
             .set_gc_content(0.58)
             .set_repeat_density(0.8)
             .add_metadata('orientation', 'forward')
             .build())
    
    logger.info(f"Created region params: {params}")
    
    # Пример создания градиента
    gradient_params = GradientParams(start_value=0.3, end_value=0.7, shape='sigmoid', steepness=5)
    gradient = GradientGenerator.create_gradient(gradient_params, 100)
    
    logger.info(f"Created gradient with shape {gradient_params.shape}")
    logger.info(f"Gradient start: {gradient[0]:.2f}, end: {gradient[-1]:.2f}")