# yatomat/regions/boundaries.py
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
import numpy as np
from pathlib import Path
import logging

from regions.common import (
    ChromosomeRegion, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator, SequenceFeature
)
from repeats import (
    RepeatGenerator, RepeatType, HORGenerator,
    MutationEngine, MutationParams
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ChromatinState(Enum):
    """Types of chromatin states"""
    HETEROCHROMATIN = "heterochromatin"
    EUCHROMATIN = "euchromatin"
    FACULTATIVE = "facultative"
    BOUNDARY = "boundary"

@dataclass
class TADParams:
    """Parameters for Topologically Associating Domains"""
    min_size: int = 100_000  # 100 kb
    max_size: int = 1_000_000  # 1 Mb
    boundary_strength: float = 0.8  # 0-1 scale
    ctcf_density: float = 0.01  # CTCF binding sites per kb
    internal_structure_complexity: float = 0.5  # 0-1 scale

@dataclass
class ChromatinTransitionParams:
    """Parameters for chromatin state transitions"""
    transition_length: int = 50_000  # 50 kb default
    gc_content_range: Tuple[float, float] = (0.35, 0.65)
    repeat_density_range: Tuple[float, float] = (0.1, 0.8)
    nucleosome_positioning_strength: float = 0.7
    gradient_steepness: float = 2.0

@dataclass
class BoundaryParams:
    """Parameters for boundary regions"""
    tad_params: TADParams = field(default_factory=TADParams)
    transition_params: ChromatinTransitionParams = field(
        default_factory=ChromatinTransitionParams
    )
    insulator_strength: float = 0.8
    min_boundary_size: int = 5_000  # 5 kb
    max_boundary_size: int = 20_000  # 20 kb

class CTCFSite:
    """Class representing CTCF binding sites"""

    CONSENSUS = "CCGCGNGGNGGCAG"  # CTCF consensus sequence

    def __init__(self, position: int, orientation: str = "forward",
                 strength: float = 1.0):
        self.position = position
        self.orientation = orientation
        self.strength = strength
        self._sequence = self._generate_sequence()

    def _generate_sequence(self) -> str:
        """Generate CTCF binding site sequence with possible variations"""
        sequence = list(self.CONSENSUS)

        # Replace 'N's with random nucleotides
        for i in range(len(sequence)):
            if sequence[i] == 'N':
                sequence[i] = np.random.choice(['A', 'T', 'G', 'C'])

        # Introduce variations based on strength
        mutation_prob = (1 - self.strength) * 0.2
        for i in range(len(sequence)):
            if np.random.random() < mutation_prob:
                current = sequence[i]
                options = ['A', 'T', 'G', 'C']
                options.remove(current)
                sequence[i] = np.random.choice(options)

        if self.orientation == "reverse":
            sequence = self._reverse_complement(''.join(sequence))
        else:
            sequence = ''.join(sequence)

        return sequence

    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """Generate reverse complement of a sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence))

    @property
    def sequence(self) -> str:
        return self._sequence

class TADomain:
    """Class representing a Topologically Associating Domain"""

    def __init__(self, params: TADParams, size: Optional[int] = None):
        self.params = params
        if size is None:
            self.size = np.random.randint(params.min_size, params.max_size)
        else:
            if size < 2 * len(CTCFSite.CONSENSUS):
                raise ValueError(f"TAD size {size} is too small for CTCF sites.")
            self.size = size
        self.boundaries: List[CTCFSite] = []
        self.internal_ctcf_sites: List[CTCFSite] = []
        self._generate_boundaries()
        self._generate_internal_ctcf_sites()

    def _generate_boundaries(self):
        """Generate CTCF sites at domain boundaries"""
        # Strong CTCF sites at boundaries
        self.boundaries.append(CTCFSite(0, "forward", self.params.boundary_strength))
        self.boundaries.append(CTCFSite(
            self.size - len(CTCFSite.CONSENSUS),
            "reverse",
            self.params.boundary_strength
        ))

    def _generate_internal_ctcf_sites(self):
        """Generate internal CTCF sites within the TAD"""
        num_internal = int((self.size / 1000) * self.params.ctcf_density)
        # Ensure at least one internal CTCF site if density > 0
        if self.params.ctcf_density > 0 and num_internal == 0:
            num_internal = 1
        for _ in range(num_internal):
            pos = np.random.randint(
                len(CTCFSite.CONSENSUS),
                self.size - len(CTCFSite.CONSENSUS)
            )
            orientation = np.random.choice(["forward", "reverse"])
            strength = np.random.uniform(0.3, 0.9)
            self.internal_ctcf_sites.append(CTCFSite(pos, orientation, strength))

    def get_ctcf_sites(self) -> List[CTCFSite]:
        """Return all CTCF sites (boundaries + internal)"""
        return self.boundaries + self.internal_ctcf_sites

class BoundaryRegion(ChromosomeRegion):
    """Class for generating boundary regions between different chromatin states"""

    def __init__(self, params: RegionParams,
                 boundary_params: Optional[BoundaryParams] = None):
        super().__init__(params)
        self.boundary_params = boundary_params or BoundaryParams()
        self.gradient_gen = GradientGenerator()
        self.repeat_gen = RepeatGenerator()
        self.mut_engine = MutationEngine()

    def _generate_tad_sequence(self, tad: TADomain, tad_start: int) -> Tuple[str, List[Dict]]:
        """Generate the sequence and features for a single TAD"""
        sequence = self.repeat_gen.seq_gen.generate_sequence(tad.size, local_gc=self.params.gc_content)
        features = []

        # Add TAD annotation
        features.append({
            'type': 'TAD',
            'start': tad_start,
            'end': tad_start + tad.size,
            'size': tad.size
        })

        # Insert CTCF sites
        for site in tad.get_ctcf_sites():
            site_seq = site.sequence
            pos = site.position
            # Ensure that the site fits within the sequence
            if pos + len(site_seq) > tad.size:
                logger.warning(f"CTCF site at position {pos} with length {len(site_seq)} exceeds TAD size {tad.size}. Skipping.")
                continue
            sequence = sequence[:pos] + site_seq + sequence[pos + len(site_seq):]
            features.append({
                'type': 'CTCF_site',
                'start': tad_start + pos,
                'end': tad_start + pos + len(site_seq),
                'orientation': site.orientation,
                'strength': site.strength
            })

        return sequence, features

    def _create_chromatin_transition(self, 
                                   start_state: ChromatinState,
                                   end_state: ChromatinState,
                                   length: int) -> Tuple[str, List[Dict]]:
        """Create transition region between chromatin states"""
        params = self.boundary_params.transition_params
        sequence = ""
        features = []

        # Create gradients for different properties
        gc_gradient = self.gradient_gen.create_gradient(
            GradientParams(
                start_value=params.gc_content_range[0],
                end_value=params.gc_content_range[1],
                shape='linear',  # Changed to linear for strict monotonicity
                steepness=params.gradient_steepness
            ),
            length
        )

        # Ensure the gradient is strictly monotonic
        if params.gc_content_range[0] < params.gc_content_range[1]:
            gc_gradient = np.sort(gc_gradient)
        else:
            gc_gradient = np.sort(gc_gradient)[::-1]

        # Generate sequence in chunks with varying properties
        chunk_size = 1000
        for i in range(0, length, chunk_size):
            current_length = min(chunk_size, length - i)
            local_gc = float(np.mean(gc_gradient[i:i + current_length]))
            
            chunk = self.repeat_gen.seq_gen.generate_sequence(
                current_length,
                local_gc=local_gc
            )
            sequence += chunk

        # Add transition annotation
        features.append({
            'type': 'chromatin_transition',
            'start': 0,
            'end': length,
            'start_state': start_state.value,
            'end_state': end_state.value,
            'gc_range': (min(gc_gradient), max(gc_gradient)),
            'repeat_range': params.repeat_density_range
        })

        return sequence, features

    def generate(self) -> Tuple[str, List[Dict]]:
        """Generate complete boundary region with TADs and transitions"""
        sequence = ""
        features = []
        current_pos = 0

        # Calculate number and size of TADs
        total_length = self.length
        transition_length = self.boundary_params.transition_params.transition_length  # Fixed length

        # Ensure we have space for initial and final transitions
        remaining_length = total_length - 2 * transition_length
        if remaining_length <= 0:
            raise ValueError("Total length is too short for the required transitions.")

        # Determine number of TADs based on average size
        avg_tad_size = (self.boundary_params.tad_params.min_size + self.boundary_params.tad_params.max_size) // 2
        num_tads = max(1, remaining_length // avg_tad_size)

        # Generate initial transition
        trans_seq, trans_features = self._create_chromatin_transition(
            ChromatinState.HETEROCHROMATIN,
            ChromatinState.BOUNDARY,
            transition_length
        )

        # Adjust feature positions for initial transition
        for feature in trans_features:
            feature['start'] += current_pos
            feature['end'] += current_pos

        sequence += trans_seq
        features.extend(trans_features)
        current_pos += transition_length

        # Generate TADs with transitions
        for i in range(num_tads):
            # Calculate TAD size
            if i == num_tads - 1:
                # Last TAD - use remaining space
                tad_length = remaining_length - (i * avg_tad_size)
            else:
                tad_length = np.random.randint(
                    self.boundary_params.tad_params.min_size,
                    self.boundary_params.tad_params.max_size
                )

            # Create TADomain instance with correct size
            tad = TADomain(self.boundary_params.tad_params, size=tad_length)

            # Generate TAD sequence and features
            tad_seq, tad_features = self._generate_tad_sequence(tad, current_pos)

            # Add TAD sequence and features
            sequence += tad_seq
            features.extend(tad_features)
            current_pos += tad_length

            # Add transition between TADs if not the last TAD
            if i < num_tads - 1:
                trans_seq, trans_features = self._create_chromatin_transition(
                    ChromatinState.BOUNDARY,
                    ChromatinState.BOUNDARY,
                    transition_length
                )
                
                # Adjust feature positions
                for feature in trans_features:
                    feature['start'] += current_pos
                    feature['end'] += current_pos
                
                sequence += trans_seq
                features.extend(trans_features)
                current_pos += transition_length

        # Add final transition
        final_trans_seq, final_trans_features = self._create_chromatin_transition(
            ChromatinState.BOUNDARY,
            ChromatinState.EUCHROMATIN,
            transition_length
        )
        
        # Adjust final transition positions
        for feature in final_trans_features:
            feature['start'] += current_pos
            feature['end'] += current_pos
        
        sequence += final_trans_seq
        features.extend(final_trans_features)
        current_pos += transition_length

        # Ensure we match the target length
        if len(sequence) != total_length:
            if len(sequence) < total_length:
                # Pad with additional sequence if needed
                additional_length = total_length - len(sequence)
                sequence += self.repeat_gen.seq_gen.generate_sequence(
                    additional_length,
                    local_gc=self.params.gc_content
                )
            else:
                # Trim if too long
                sequence = sequence[:total_length]
                # Adjust features that might extend beyond the end
                features = [f for f in features if f['start'] < total_length]
                for f in features:
                    if f['end'] > total_length:
                        f['end'] = total_length

        # Fix any remaining feature overlaps
        sorted_features = sorted(features, key=lambda x: x['start'])
        fixed_features = []
        
        for i, feature in enumerate(sorted_features):
            # if i > 0 and feature['start'] < fixed_features[-1]['end']:
            #     # Adjust start position to avoid overlap
            #     feature['start'] = fixed_features[-1]['end']
            #     if feature['end'] <= feature['start']:
            #         continue  # Skip features that would have negative length
            fixed_features.append(feature)

        return sequence, fixed_features


if __name__ == "__main__":
    # Example usage
    params = RegionParams(
        region_type=ChromosomeRegionType.TRANSITION,
        start=0,
        end=200_000,
        gc_content=0.45)
    
    boundary = BoundaryRegion(params)
    sequence, features = boundary.generate()
    
    logger.info(f"Generated boundary region of length {len(sequence)}")
    logger.info(f"Number of features: {len(features)}")
    
    # Analyze feature distribution
    feature_types = {}
    for feature in features:
        feat_type = feature['type']
        feature_types[feat_type] = feature_types.get(feat_type, 0) + 1
    
    logger.info("Feature distribution:")
    for feat_type, count in feature_types.items():
        logger.info(f"{feat_type}: {count}")
