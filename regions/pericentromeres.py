# yatomat/regions/pericentromeres.py
from typing import List, Dict, Optional, Tuple, Union, cast
from dataclasses import dataclass, field
from enum import Enum
import numpy as np
from numpy.typing import ArrayLike
from pathlib import Path
import logging

from regions.common import (
    ChromosomeRegion, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator, SequenceFeature
)
from repeats import (
    RepeatGenerator, RepeatType, HORGenerator,
    MutationEngine, HORParams, MutationParams
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PericentromereZone(Enum):
    """Types of pericentromeric zones"""
    PROXIMAL = "proximal"  # Closest to centromere
    INTERMEDIATE = "intermediate"  # Middle region
    DISTAL = "distal"  # Furthest from centromere
    TRANSITION = "transition"  # Transition zones

@dataclass
class SatelliteDistributionParams:
    """Parameters for satellite distribution in pericentromeric regions"""
    alpha_satellite_density: float = 0.3  # Proportion of alpha satellites
    beta_satellite_density: float = 0.2   # Proportion of beta satellites
    gamma_satellite_density: float = 0.15  # Proportion of other satellites
    satellite_mutation_rate: float = 0.1   # Base mutation rate for satellites
    min_block_size: int = 1000            # Minimum size of satellite blocks
    max_block_size: int = 10000           # Maximum size of satellite blocks

@dataclass
class MobileElementParams:
    """Parameters for mobile element integration"""
    density: float = 0.2                  # Overall density of mobile elements
    min_size: int = 500                   # Minimum size of mobile elements
    max_size: int = 5000                  # Maximum size of mobile elements
    types: List[str] = field(default_factory=lambda: ['LINE', 'SINE', 'LTR', 'DNA'])

@dataclass
class PericentromereParams:
    """Parameters for pericentromeric region generation"""
    min_total_length: int = 2_000_000     # 2 Mb minimum
    max_total_length: int = 5_000_000     # 5 Mb maximum
    transition_length: int = 100_000      # 100 kb transitions

    # Parameters for different components
    satellite_params: SatelliteDistributionParams = field(
        default_factory=SatelliteDistributionParams
    )
    mobile_element_params: MobileElementParams = field(
        default_factory=MobileElementParams
    )

    # GC content parameters
    proximal_gc: float = 0.56             # GC content near centromere
    distal_gc: float = 0.52              # GC content far from centromere
    gc_std: float = 0.03                 # Standard deviation of GC content

class PericentromereRegion(ChromosomeRegion):
    """Class for generating pericentromeric regions"""

    def __init__(self, params: RegionParams,
                 pericentromere_params: Optional[PericentromereParams] = None):
        super().__init__(params)
        self.pericentromere_params = pericentromere_params or PericentromereParams()
        self.repeat_gen = RepeatGenerator()
        self.gradient_gen = GradientGenerator()
        self.hor_gen = HORGenerator()

    def _generate_mobile_element(self, size: int) -> Tuple[str, Dict[str, Union[str, int, float]]]:
        """Generate a mobile element sequence and its annotation"""
        params = self.pericentromere_params.mobile_element_params
        element_type = np.random.choice(params.types)

        # Generate sequence with specific nucleotide bias based on element type
        if element_type in ['LINE', 'LTR']:
            gc_content = 0.42  # AT-rich
        else:
            gc_content = 0.53  # Slightly GC-rich

        sequence = self.repeat_gen.seq_gen.generate_sequence(size, local_gc=gc_content)

        return sequence, {
            'type': 'mobile_element',
            'element_type': element_type,
            'length': size,
            'gc_content': gc_content
        }

    def _generate_satellite_block(self, size: int, zone: PericentromereZone) -> Tuple[str, Dict[str, Union[str, int, float]]]:
        """Generate a satellite repeat block"""
        params = self.pericentromere_params.satellite_params

        # Adjust satellite probabilities based on zone
        if zone == PericentromereZone.PROXIMAL:
            alpha_prob = params.alpha_satellite_density * 1.5
            beta_prob = params.beta_satellite_density
            gamma_prob = params.gamma_satellite_density * 0.5
        elif zone == PericentromereZone.DISTAL:
            alpha_prob = params.alpha_satellite_density * 0.5
            beta_prob = params.beta_satellite_density
            gamma_prob = params.gamma_satellite_density * 1.5
        else:
            alpha_prob = params.alpha_satellite_density
            beta_prob = params.beta_satellite_density
            gamma_prob = params.gamma_satellite_density

        # Normalize probabilities
        total = alpha_prob + beta_prob + gamma_prob
        probs = np.array([alpha_prob/total, beta_prob/total, gamma_prob/total])

        # Choose satellite type
        satellite_options = np.array([
            RepeatType.ALPHA_SATELLITE,
            RepeatType.BETA_SATELLITE,
            RepeatType.SATELLITE_1
        ])
        satellite_type = cast(RepeatType, np.random.choice(satellite_options, p=probs))

        # Generate base monomer
        monomer = self.repeat_gen.generate_monomer(satellite_type)

        # Calculate number of copies needed
        copies = size // len(monomer.sequence)
        sequence = monomer.sequence * copies

        # Apply mutations with rate depending on distance from centromere
        if zone == PericentromereZone.PROXIMAL:
            mutation_rate = params.satellite_mutation_rate * 0.5
        elif zone == PericentromereZone.DISTAL:
            mutation_rate = params.satellite_mutation_rate * 2.0
        else:
            mutation_rate = params.satellite_mutation_rate

        # Apply mutations
        mut_params = MutationParams(
            substitution_rate=mutation_rate,
            insertion_rate=mutation_rate / 5,
            deletion_rate=mutation_rate / 5
        )
        mut_engine = MutationEngine(mut_params)
        sequence, mutations = mut_engine.mutate_sequence(sequence)

        return sequence[:size], {
            'type': 'satellite',
            'satellite_type': satellite_type.value,
            'monomer_size': len(monomer.sequence),
            'copies': copies,
            'mutations': len(mutations),
            'zone': zone.value
        }

    def _create_zone_gradient(self, start_value: float, end_value: float,
                            length: int) -> np.ndarray:
        """Create a gradient for transitioning parameters within a zone"""
        params = GradientParams(
            start_value=start_value,
            end_value=end_value,
            shape='sigmoid',
            steepness=3.0
        )
        return self.gradient_gen.create_gradient(params, length)

    def _generate_zone(self, zone: PericentromereZone,
                      length: int) -> Tuple[str, List[Dict[str, Union[str, int, float]]]]:
        """Generate a specific zone of the pericentromeric region"""
        sequence = ""
        features: List[Dict[str, Union[str, int, float]]] = []
        current_pos = 0

        # Create gradient for GC content
        if zone == PericentromereZone.PROXIMAL:
            gc_gradient = np.full(length, self.pericentromere_params.proximal_gc)
        elif zone == PericentromereZone.DISTAL:
            gc_gradient = np.full(length, self.pericentromere_params.distal_gc)
        else:
            gc_gradient = self._create_zone_gradient(
                self.pericentromere_params.proximal_gc,
                self.pericentromere_params.distal_gc,
                length
            )

        while current_pos < length:
            # Determine block type and size
            if np.random.random() < self.pericentromere_params.mobile_element_params.density:
                # Generate mobile element
                size = np.random.randint(
                    self.pericentromere_params.mobile_element_params.min_size,
                    self.pericentromere_params.mobile_element_params.max_size
                )
                block_seq, block_info = self._generate_mobile_element(size)
            else:
                # Generate satellite block
                size = np.random.randint(
                    self.pericentromere_params.satellite_params.min_block_size,
                    self.pericentromere_params.satellite_params.max_block_size
                )
                block_seq, block_info = self._generate_satellite_block(size, zone)

            # Ensure we don't exceed target length
            if current_pos + len(block_seq) > length:
                block_seq = block_seq[:length - current_pos]

            # Add sequence and feature
            sequence += block_seq
            block_info.update({
                'start': current_pos,
                'end': current_pos + len(block_seq),
                'zone': str(zone.value)
            })
            features.append(block_info)
            current_pos += len(block_seq)

        return sequence, features

    def generate(self) -> Tuple[str, List[Dict[str, Union[str, int, float]]]]:
        """Generate complete pericentromeric region"""
        # Determine total length
        total_length = np.random.randint(
            self.pericentromere_params.min_total_length,
            self.pericentromere_params.max_total_length
        )

        # Calculate zone lengths
        main_length = total_length - 2 * self.pericentromere_params.transition_length
        zone_lengths = {
            PericentromereZone.PROXIMAL: int(main_length * 0.4),
            PericentromereZone.INTERMEDIATE: int(main_length * 0.3),
            PericentromereZone.DISTAL: int(main_length * 0.3),
            PericentromereZone.TRANSITION: self.pericentromere_params.transition_length
        }

        # Generate each zone
        sequence = ""
        features: List[Dict[str, Union[str, int, float]]] = []
        current_pos = 0

        for zone in [PericentromereZone.PROXIMAL, PericentromereZone.INTERMEDIATE,
                    PericentromereZone.DISTAL]:
            # Generate main zone
            zone_seq, zone_features = self._generate_zone(zone, zone_lengths[zone])

            # Add transition zone if not the last zone
            if zone != PericentromereZone.DISTAL:
                trans_seq, trans_features = self._generate_zone(
                    PericentromereZone.TRANSITION,
                    zone_lengths[PericentromereZone.TRANSITION]
                )
                zone_seq += trans_seq
                zone_features.extend(trans_features)

            # Update positions and add features
            for feature in zone_features:
                feature['start'] = current_pos + int(feature['start'])
                feature['end'] = current_pos + int(feature['end'])

            sequence += zone_seq
            features.extend(zone_features)
            current_pos += len(zone_seq)

        return sequence, features

if __name__ == "__main__":
    # Example usage
    params = RegionParams(
        region_type=ChromosomeRegionType.PERICENTROMERE,
        start=0,
        end=5_000_000,
        gc_content=0.54
    )

    pericentromere = PericentromereRegion(params, PericentromereParams())
    sequence, features = pericentromere.generate()

    logger.info(f"Generated pericentromeric sequence of length {len(sequence)}")
    logger.info(f"Number of features: {len(features)}")

    # Analyze feature distribution
    feature_types: Dict[str, int] = {}
    for feature in features:
        feat_type = str(feature['type'])
        feature_types[feat_type] = feature_types.get(feat_type, 0) + 1

    logger.info("Feature distribution:")
    for feat_type, count in feature_types.items():
        logger.info(f"{feat_type}: {count}")
