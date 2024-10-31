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
        return ''.join(complement[base] for base in reversed(sequence))

    @property
    def sequence(self) -> str:
        return self._sequence

class TADomain:
    """Class representing a Topologically Associating Domain"""

    def __init__(self, params: TADParams):
        self.params = params
        self.size = np.random.randint(params.min_size, params.max_size)
        self.boundaries: List[CTCFSite] = []
        self.internal_features: List[Dict] = []
        self._generate_boundaries()

    def _generate_boundaries(self):
        """Generate CTCF sites at domain boundaries"""
        # Strong CTCF sites at boundaries
        self.boundaries.append(CTCFSite(0, "forward", self.params.boundary_strength))
        self.boundaries.append(CTCFSite(
            self.size - len(CTCFSite.CONSENSUS),
            "reverse",
            self.params.boundary_strength
        ))

        # Internal CTCF sites
        num_internal = int(self.size * self.params.ctcf_density / 1000)
        for _ in range(num_internal):
            pos = np.random.randint(
                len(CTCFSite.CONSENSUS),
                self.size - len(CTCFSite.CONSENSUS)
            )
            orientation = np.random.choice(["forward", "reverse"])
            strength = np.random.uniform(0.3, 0.9)
            self.boundaries.append(CTCFSite(pos, orientation, strength))

class BoundaryRegion(ChromosomeRegion):
    """Class for generating boundary regions between different chromatin states"""

    def __init__(self, params: RegionParams,
                 boundary_params: Optional[BoundaryParams] = None):
        super().__init__(params)
        self.boundary_params = boundary_params or BoundaryParams()
        self.gradient_gen = GradientGenerator()
        self.repeat_gen = RepeatGenerator()
        self.mut_engine = MutationEngine()

    def _apply_ctcf_sites(self, sequence: str, tad_start: int, tad_end: int) -> Tuple[str, List[Dict]]:
        """Generate and place CTCF binding sites for a TAD"""
        sites = []
        length = len(sequence)
        
        # Start boundary - forward orientation
        start_site = CTCFSite(0, "forward", self.boundary_params.tad_params.boundary_strength)
        sequence = sequence[:0] + start_site.sequence + sequence[len(start_site.sequence):]
        sites.append({
            'type': 'CTCF_site',
            'start': tad_start + 0,
            'end': tad_start + len(start_site.sequence),
            'orientation': 'forward',
            'strength': start_site.strength
        })

        # End boundary - reverse orientation
        end_pos = length - len(CTCFSite.CONSENSUS)
        end_site = CTCFSite(end_pos, "reverse", self.boundary_params.tad_params.boundary_strength)
        sequence = sequence[:end_pos] + end_site.sequence + sequence[end_pos + len(end_site.sequence):]
        sites.append({
            'type': 'CTCF_site',
            'start': tad_start + end_pos,
            'end': tad_start + end_pos + len(end_site.sequence),
            'orientation': 'reverse',
            'strength': end_site.strength
        })

        # Internal CTCF sites
        num_internal = int(length * self.boundary_params.tad_params.ctcf_density / 1000)
        for _ in range(num_internal):
            pos = np.random.randint(
                len(CTCFSite.CONSENSUS),
                length - len(CTCFSite.CONSENSUS)
            )
            orientation = np.random.choice(["forward", "reverse"])
            strength = np.random.uniform(0.3, 0.9)

            site = CTCFSite(pos, orientation, strength)
            sequence = sequence[:pos] + site.sequence + sequence[pos + len(site.sequence):]

            sites.append({
                'type': 'CTCF_site',
                'start': tad_start + pos,
                'end': tad_start + pos + len(site.sequence),
                'orientation': orientation,
                'strength': strength
            })

        return sequence, sites

    def _generate_tad(self, start_pos: int, length: int) -> Tuple[str, List[Dict]]:
        """Generate a TAD with appropriate features"""
        sequence = self.repeat_gen.seq_gen.generate_sequence(length, local_gc=self.params.gc_content)
        features = []

        # Add TAD annotation
        features.append({
            'type': 'TAD',
            'start': start_pos,
            'end': start_pos + length,
            'size': length
        })

        # Add CTCF sites
        sequence, ctcf_features = self._apply_ctcf_sites(sequence, start_pos, start_pos + length)
        features.extend(ctcf_features)

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
                shape='sigmoid',
                steepness=params.gradient_steepness
            ),
            length
        )

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
        min_tad_size = self.boundary_params.tad_params.min_size
        max_tad_size = min(self.boundary_params.tad_params.max_size, total_length // 3)
        transition_length = self.boundary_params.transition_params.transition_length  # Fixed line


        # Ensure we have space for at least one TAD and two transitions
        remaining_length = total_length - 2 * transition_length
        target_tad_size = (min_tad_size + max_tad_size) // 2
        num_tads = max(1, remaining_length // target_tad_size)

        # Generate initial transition
        trans_seq, trans_features = self._create_chromatin_transition(
            ChromatinState.HETEROCHROMATIN,
            ChromatinState.BOUNDARY,
            transition_length
        )

        # Generate TADs
        for i in range(num_tads):
            # Calculate TAD size
            if i == num_tads - 1:
                # Last TAD - use remaining space
                tad_length = remaining_length - current_pos
            else:
                tad_length = np.random.randint(min_tad_size, max_tad_size)

            # Generate TAD
            tad_seq, tad_features = self._generate_tad(current_pos, tad_length)

            # Add initial transition features after adjusting positions
            if i == 0:
                for feature in trans_features:
                    feature['start'] += current_pos
                    feature['end'] += current_pos
                sequence += trans_seq
                features.extend(trans_features)
                current_pos += transition_length
            
            # Add transition between TADs
            elif i > 0:
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

            # Add TAD sequence and features
            sequence += tad_seq
            features.extend(tad_features)
            current_pos += tad_length

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
            if i > 0 and feature['start'] < fixed_features[-1]['end']:
                # Adjust start position to avoid overlap
                feature['start'] = fixed_features[-1]['end']
                if feature['end'] <= feature['start']:
                    continue  # Skip features that would have negative length
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
