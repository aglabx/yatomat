#yatomat/regions/chromosome_regions.py
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
import numpy as np
from pathlib import Path
import logging

from yatomat.regions.common import (
    ChromosomeRegion, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator, SequenceFeature
)
from yatomat.regions.telomers import (
    TelomereRegion, SubtelomereRegion, TelomereParams, SubtelomereParams
)
from yatomat.regions.centromeres import (
    CentromereRegion, CentromereParams, CentromereZone
)
from yatomat.regions.pericentromeres import (
    PericentromereRegion, PericentromereParams, PericentromereZone
)
from yatomat.regions.boundaries import (
    BoundaryRegion, BoundaryParams, ChromatinState
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class ChromosomeArmParams:
    """Parameters for a chromosome arm"""
    length: int  # Total arm length
    gc_content: float = 0.45  # Average GC content
    gc_std: float = 0.02  # Standard deviation of GC content
    repeat_density: float = 0.3  # Overall repeat density
    tad_organization: float = 0.8  # Strength of TAD organization (0-1)
    chromatin_state: ChromatinState = ChromatinState.EUCHROMATIN

@dataclass
class ChromosomeParams:
    """Parameters for entire chromosome generation"""
    total_length: int  # Total chromosome length
    p_arm_params: ChromosomeArmParams  # Parameters for p arm
    q_arm_params: ChromosomeArmParams  # Parameters for q arm

    # Parameters for specific regions
    telomere_params: TelomereParams = field(
        default_factory=lambda: TelomereParams()
    )
    subtelomere_params: SubtelomereParams = field(
        default_factory=lambda: SubtelomereParams()
    )
    centromere_params: CentromereParams = field(
        default_factory=lambda: CentromereParams()
    )
    pericentromere_params: PericentromereParams = field(
        default_factory=lambda: PericentromereParams()
    )
    boundary_params: BoundaryParams = field(
        default_factory=lambda: BoundaryParams()
    )

    def __post_init__(self):
        """Ensure all parameters are initialized"""
        if self.telomere_params is None:
            self.telomere_params = TelomereParams()
        if self.subtelomere_params is None:
            self.subtelomere_params = SubtelomereParams()
        if self.centromere_params is None:
            self.centromere_params = CentromereParams()
        if self.pericentromere_params is None:
            self.pericentromere_params = PericentromereParams()
        if self.boundary_params is None:
            self.boundary_params = BoundaryParams()

        # Validate total length against arm lengths
        total_arms = self.p_arm_params.length + self.q_arm_params.length
        if total_arms > self.total_length:
            raise ValueError(
                f"Combined arm lengths ({total_arms}) exceed total length ({self.total_length})"
            )

class ChromosomeAssembly:
    """Class for assembling complete chromosome sequences"""

    def __init__(self, params: ChromosomeParams):
        # Ensure all parameters are initialized
        self.params = params

        # Убедимся что все параметры инициализированы
        if self.params.telomere_params is None:
            self.params.telomere_params = TelomereParams()
        if self.params.subtelomere_params is None:
            self.params.subtelomere_params = SubtelomereParams()
        if self.params.centromere_params is None:
            self.params.centromere_params = CentromereParams()
        if self.params.pericentromere_params is None:
            self.params.pericentromere_params = PericentromereParams()

        # Расчет базовых метрик с проверкой
        self.centromere_size = int(min(
            self.params.centromere_params.max_core_length,
            self.params.total_length * 0.05
        ))

        # Initialize generators
        self.gradient_gen = GradientGenerator()

        # Calculate basic metrics
        self.centromere_size = int(min(
            params.centromere_params.max_core_length,
            params.total_length * 0.05  # Max 5% of total length
        ))

        # Adjust arm lengths to accommodate centromere
        remaining_length = params.total_length - self.centromere_size
        p_arm_fraction = params.p_arm_params.length / (
            params.p_arm_params.length + params.q_arm_params.length
        )
        self.p_arm_length = int(remaining_length * p_arm_fraction)
        self.q_arm_length = int(remaining_length * (1 - p_arm_fraction))

        # Initialize region generators
        self._init_region_generators()

        # Track features and boundaries
        self.features: List[Dict] = []
        self.current_position: int = 0

        logger.info(f"Initialized chromosome assembly:")
        logger.info(f"Total length: {params.total_length}")
        logger.info(f"P arm length: {self.p_arm_length}")
        logger.info(f"Centromere size: {self.centromere_size}")
        logger.info(f"Q arm length: {self.q_arm_length}")

    def _init_region_generators(self):
        """Initialize all region generators"""
        # Initialize with minimal valid sizes - they will be updated during generation
        self.telomere_gen = TelomereRegion(
            RegionParams(
                region_type=ChromosomeRegionType.TELOMERE,
                start=0, end=1000,  # Minimal initial size
                gc_content=self.params.p_arm_params.gc_content
            ),
            self.params.telomere_params or TelomereParams()
        )

        self.subtelomere_gen = SubtelomereRegion(
            RegionParams(
                region_type=ChromosomeRegionType.SUBTELOMERE,
                start=0, end=1000,
                gc_content=self.params.p_arm_params.gc_content
            ),
            self.params.subtelomere_params or SubtelomereParams()
        )

        self.centromere_gen = CentromereRegion(
            RegionParams(
                region_type=ChromosomeRegionType.CENTROMERE,
                start=0, end=1000,
                gc_content=0.58  # Typical centromeric GC content
            ),
            self.params.centromere_params or CentromereParams()
        )

        self.pericentromere_gen = PericentromereRegion(
            RegionParams(
                region_type=ChromosomeRegionType.PERICENTROMERE,
                start=0, end=1000,
                gc_content=0.56  # Typical pericentromeric GC content
            ),
            self.params.pericentromere_params or PericentromereParams()
        )

        self.boundary_gen = BoundaryRegion(
            RegionParams(
                region_type=ChromosomeRegionType.TRANSITION,
                start=0, end=1000,
                gc_content=self.params.p_arm_params.gc_content
            ),
            self.params.boundary_params or BoundaryParams()
        )

    def _update_region_params(self, generator: ChromosomeRegion,
                            start: int, end: int, gc_content: float):
        """Update region parameters for a generator"""
        generator.params.start = start
        generator.params.end = end
        generator.params.gc_content = gc_content

    def _add_features(self, features: List[Dict], offset: int = 0):
        """Add features with position adjustment"""
        if offset > 0:
            for feature in features:
                feature['start'] += offset
                feature['end'] += offset
        self.features.extend(features)

    def _generate_arm(self, arm_params: ChromosomeArmParams,
                     is_p_arm: bool = True) -> Tuple[str, List[Dict]]:
        """Generate a chromosome arm with appropriate organization"""
        sequence = ""
        features = []
        current_pos = self.current_position

        # Generate telomere
        self._update_region_params(
            self.telomere_gen,
            current_pos,
            current_pos + self.params.telomere_params.max_length,
            arm_params.gc_content
        )
        if is_p_arm:
            self.telomere_gen.params.orientation = "forward"
            tel_seq, tel_features = self.telomere_gen.generate()
            sequence += tel_seq
            features.extend(tel_features)
            current_pos += len(tel_seq)

        # Generate subtelomere
        self._update_region_params(
            self.subtelomere_gen,
            current_pos,
            current_pos + self.params.subtelomere_params.max_length,
            arm_params.gc_content
        )
        subtel_seq, subtel_features = self.subtelomere_gen.generate()
        sequence += subtel_seq
        features.extend(subtel_features)
        current_pos += len(subtel_seq)

        # Generate pericentromere if this is adjacent to centromere
        remaining_length = (self.p_arm_length if is_p_arm else self.q_arm_length) - len(sequence)
        if is_p_arm:  # Add pericentromere at end for p arm
            peri_length = min(
                remaining_length // 3,
                self.params.pericentromere_params.min_total_length
            )
            self._update_region_params(
                self.pericentromere_gen,
                current_pos,
                current_pos + peri_length,
                0.56  # Pericentromeric GC content
            )
            peri_seq, peri_features = self.pericentromere_gen.generate()
            sequence += peri_seq
            features.extend(peri_features)
            current_pos += len(peri_seq)
            remaining_length -= peri_length

        # Fill remaining space with appropriate chromatin structure
        if remaining_length > 0:
            # Create gradient for smooth transition
            gc_gradient = self.gradient_gen.create_gradient(
                GradientParams(
                    start_value=arm_params.gc_content,
                    end_value=arm_params.gc_content + arm_params.gc_std,
                    shape='sigmoid'
                ),
                remaining_length
            )

            # Generate sequence in chunks
            chunk_size = 50000  # 50 kb chunks
            current_state = arm_params.chromatin_state

            for i in range(0, remaining_length, chunk_size):
                chunk_length = min(chunk_size, remaining_length - i)
                local_gc = float(np.mean(gc_gradient[i:i + chunk_length]))

                # Possibly switch chromatin state
                if np.random.random() < 0.2:  # 20% chance to switch
                    new_state = (ChromatinState.HETEROCHROMATIN
                               if current_state == ChromatinState.EUCHROMATIN
                               else ChromatinState.EUCHROMATIN)

                    # Add transition region feature
                    features.append({
                        'type': 'chromatin_transition',
                        'start': current_pos,
                        'end': current_pos + chunk_length,
                        'from_state': current_state.value,
                        'to_state': new_state.value
                    })
                    current_state = new_state

                # Generate chunk with appropriate properties
                chunk_seq = self.telomere_gen.repeat_gen.seq_gen.generate_sequence(
                    chunk_length,
                    local_gc=local_gc
                )
                sequence += chunk_seq

                # Add feature annotations
                features.append({
                    'type': 'chromatin_region',
                    'start': current_pos,
                    'end': current_pos + chunk_length,
                    'chromatin_state': current_state.value,
                    'gc_content': local_gc
                })

                # Add TAD feature
                features.append({
                    'type': 'TAD',
                    'start': current_pos,
                    'end': current_pos + chunk_length,
                    'chromatin_state': current_state.value
                })

                current_pos += chunk_length

        if not is_p_arm:
            self.telomere_gen.telomere_params.orientation = "reverse"
            tel_seq, tel_features = self.telomere_gen.generate()
            sequence += tel_seq
            features.extend(tel_features)
            current_pos += len(tel_seq)

        return sequence, features

    def generate(self) -> Tuple[str, List[Dict]]:
        """Generate complete chromosome sequence with correct global coordinates"""
        sequence = ""
        self.features = []
        self.current_position = 0

        # 1. Generate P arm
        p_arm_seq, p_arm_features = self._generate_arm(
            self.params.p_arm_params,
            is_p_arm=True,
        )
        sequence += p_arm_seq
        self._add_features(p_arm_features)
        self.current_position += len(p_arm_seq)

        # 2. Generate centromere
        self._update_region_params(
            self.centromere_gen,
            self.current_position,
            self.current_position + self.centromere_size,
            0.58  # Typical centromeric GC content
        )
        cent_seq, cent_features = self.centromere_gen.generate()
        sequence += cent_seq
        for f in cent_features:
            f['start'] += self.current_position
            f['end'] += self.current_position
        self._add_features(cent_features)
        self.current_position += len(cent_seq)

        # 3. Generate Q arm
        q_arm_seq, q_arm_features = self._generate_arm(
            self.params.q_arm_params,
            is_p_arm=False
        )
        for f in q_arm_features:
            f['start'] += self.current_position
            f['end'] += self.current_position
        sequence += q_arm_seq
        self._add_features(q_arm_features)
        self.current_position += len(q_arm_seq)

        # # Final length check and adjustment if needed
        # if len(sequence) != self.params.total_length:
        #     if len(sequence) < self.params.total_length:
        #         additional_length = self.params.total_length - len(sequence)
        #         additional_seq = self.telomere_gen.repeat_gen.seq_gen.generate_sequence(
        #             additional_length,
        #             local_gc=self.params.q_arm_params.gc_content
        #         )

        #         self.features.append({
        #             'type': 'chromatin_region',
        #             'start': self.current_position,
        #             'end': self.current_position + additional_length,
        #             'chromatin_state': ChromatinState.EUCHROMATIN.value,
        #             'gc_content': self.params.q_arm_params.gc_content
        #         })

        #         sequence += additional_seq
        #     else:
        #         sequence = sequence[:self.params.total_length]
        #         self.features = [f for f in self.features if f['start'] < self.params.total_length]
        #         for f in self.features:
        #             if f['end'] > self.params.total_length:
        #                 f['end'] = self.params.total_length

        return sequence, self.features

if __name__ == "__main__":
    # Example usage
    p_arm_params = ChromosomeArmParams(
        length=30_000_000,  # 30 Mb
        gc_content=0.45,
        gc_std=0.02,
        repeat_density=0.3,
        tad_organization=0.8,
        chromatin_state=ChromatinState.EUCHROMATIN
    )

    q_arm_params = ChromosomeArmParams(
        length=50_000_000,  # 50 Mb
        gc_content=0.48,
        gc_std=0.02,
        repeat_density=0.35,
        tad_organization=0.75,
        chromatin_state=ChromatinState.EUCHROMATIN
    )

    chromosome_params = ChromosomeParams(
        total_length=80_000_000,  # 80 Mb total
        p_arm_params=p_arm_params,
        q_arm_params=q_arm_params
    )

    assembly = ChromosomeAssembly(chromosome_params)
    sequence, features = assembly.generate()

    logger.info(f"Generated chromosome of length: {len(sequence)}")
    logger.info(f"Number of features: {len(features)}")

    # Count feature types
    feature_types = {}
    for feature in features:
        feat_type = feature['type']
        feature_types[feat_type] = feature_types.get(feat_type, 0) + 1

    logger.info("\nFeature distribution:")
    for feat_type, count in feature_types.items():
        logger.info(f"{feat_type}: {count}")
