import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import logging
from collections import Counter, defaultdict
from typing import List, Dict, Tuple

from regions.chromosome_regions import (
    ChromosomeAssembly, ChromosomeParams, ChromosomeArmParams,
    ChromatinState
)
from regions.common import (
    RegionBuilder, ChromosomeRegionType
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TestChromosomeAssembly(unittest.TestCase):
    """Tests for chromosome assembly and generation"""
    
    def setUp(self):
        """Initialize test parameters"""
        # Create arm parameters
        self.p_arm_params = ChromosomeArmParams(
            length=30_000_000,  # 30 Mb
            gc_content=0.45,
            gc_std=0.02,
            repeat_density=0.3,
            tad_organization=0.8,
            chromatin_state=ChromatinState.EUCHROMATIN
        )
        
        self.q_arm_params = ChromosomeArmParams(
            length=50_000_000,  # 50 Mb
            gc_content=0.48,
            gc_std=0.02,
            repeat_density=0.35,
            tad_organization=0.75,
            chromatin_state=ChromatinState.EUCHROMATIN
        )
        
        # Create chromosome parameters
        self.chromosome_params = ChromosomeParams(
            total_length=80_000_000,  # 80 Mb total
            p_arm_params=self.p_arm_params,
            q_arm_params=self.q_arm_params
        )
        
        # Create assembly
        self.assembly = ChromosomeAssembly(self.chromosome_params)

    def test_initialization(self):
        """Test proper initialization of chromosome assembly"""
        # Test that all parameters are properly initialized
        self.assertIsNotNone(self.assembly.params.telomere_params)
        self.assertIsNotNone(self.assembly.params.centromere_params)
        self.assertIsNotNone(self.assembly.params.pericentromere_params)
        
        # Test length calculations
        total_structural = (
            self.assembly.p_arm_length + 
            self.assembly.centromere_size + 
            self.assembly.q_arm_length
        )
        self.assertEqual(
            total_structural,
            self.chromosome_params.total_length,
            "Total structural length doesn't match target length"
        )
        
        # Test arm proportions
        arm_ratio = self.assembly.p_arm_length / self.assembly.q_arm_length
        expected_ratio = self.p_arm_params.length / self.q_arm_params.length
        self.assertAlmostEqual(
            arm_ratio,
            expected_ratio,
            places=2,
            msg="Arm ratio not preserved"
        )

    def test_sequence_generation(self):
        """Test basic sequence generation properties"""
        sequence, features = self.assembly.generate()
        
        # Test sequence length
        self.assertEqual(
            len(sequence),
            self.chromosome_params.total_length,
            "Generated sequence length doesn't match target"
        )
        
        # Test sequence composition
        valid_nucleotides = set('ATGC')
        self.assertTrue(
            all(n in valid_nucleotides for n in sequence),
            "Invalid nucleotides in sequence"
        )
        
        # Test basic GC content
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence)
        expected_gc = (
            self.p_arm_params.gc_content * self.assembly.p_arm_length +
            self.q_arm_params.gc_content * self.assembly.q_arm_length
        ) / (self.assembly.p_arm_length + self.assembly.q_arm_length)
        
        self.assertAlmostEqual(
            gc_content,
            expected_gc,
            places=2,
            msg="Overall GC content significantly different from expected"
        )

    def test_feature_integrity(self):
        """Test integrity and organization of chromosome features"""
        sequence, features = self.assembly.generate()
        
        # Test that we have all essential features
        feature_types = set(f['type'] for f in features)
        essential_types = {
            'telomere', 'centromere', 'TAD', 'chromatin_region'
        }
        for feature_type in essential_types:
            self.assertIn(
                feature_type,
                feature_types,
                f"Missing essential feature type: {feature_type}"
            )
        
        # Test feature positions
        for feature in features:
            self.assertGreaterEqual(feature['start'], 0)
            self.assertLess(feature['start'], len(sequence))
            self.assertGreater(feature['end'], feature['start'])
            self.assertLessEqual(feature['end'], len(sequence))
        
        # Test feature ordering and overlap
        sorted_features = sorted(features, key=lambda x: x['start'])
        for i in range(len(sorted_features) - 1):
            current = sorted_features[i]
            next_feat = sorted_features[i + 1]
            
            # Allow small overlaps only for transition features
            if (current['type'] not in {'chromatin_transition', 'TAD'} and 
                next_feat['type'] not in {'chromatin_transition', 'TAD'}):
                self.assertLessEqual(
                    current['end'],
                    next_feat['start'],
                    "Invalid feature overlap"
                )

    def test_region_organization(self):
        """Test proper organization of chromosome regions"""
        sequence, features = self.assembly.generate()
        
        # Get sorted features
        sorted_features = sorted(features, key=lambda x: x['start'])
        
        # Test telomere placement
        telomeres = [f for f in features if f['type'] == 'telomere']
        self.assertEqual(len(telomeres), 2, "Should have exactly 2 telomeres")
        
        # Find centromere
        centromeres = [f for f in features if f['type'] == 'centromere']
        self.assertEqual(len(centromeres), 1, "Should have exactly 1 centromere")
        centromere = centromeres[0]
        
        # Test centromere position
        centromere_center = (centromere['start'] + centromere['end']) / 2
        expected_center = len(sequence) * self.assembly.p_arm_length / (
            self.assembly.p_arm_length + self.assembly.q_arm_length
        )
        
        self.assertAlmostEqual(
            centromere_center,
            expected_center,
            delta=len(sequence) * 0.1,  # Allow 10% deviation
            msg="Centromere not properly positioned"
        )

    def test_chromatin_states(self):
        """Test chromatin state distribution and transitions"""
        sequence, features = self.assembly.generate()
        
        # Get chromatin state features
        chromatin_features = [
            f for f in features
            if 'chromatin_state' in f
        ]
        
        # Test that we have both eu- and heterochromatin
        states = set(f['chromatin_state'] for f in chromatin_features)
        self.assertGreaterEqual(
            len(states),
            2,
            "Should have multiple chromatin states"
        )
        
        # Test state transitions
        sorted_features = sorted(
            [f for f in features if 'chromatin_state' in f],
            key=lambda x: x['start']
        )
        
        for i in range(len(sorted_features) - 1):
            current = sorted_features[i]
            next_feat = sorted_features[i + 1]
            
            if (current['chromatin_state'] != next_feat['chromatin_state']):
                # Should have a transition feature between different states
                transition_features = [
                    f for f in features
                    if f['type'] == 'chromatin_transition'
                    and f['start'] >= current['end']
                    and f['end'] <= next_feat['start']
                ]
                
                self.assertGreaterEqual(
                    len(transition_features),
                    1,
                    "Missing transition between different chromatin states"
                )

    def test_gc_content_distribution(self):
        """Test GC content distribution along the chromosome"""
        sequence, features = self.assembly.generate()
        
        # Calculate GC content in windows
        window_size = 10000  # 10 kb windows
        gc_profile = []
        positions = []
        
        for i in range(0, len(sequence), window_size):
            window = sequence[i:i + window_size]
            gc_count = window.count('G') + window.count('C')
            gc_content = gc_count / len(window)
            gc_profile.append(gc_content)
            positions.append(i + window_size/2)
        
        # Visualize GC content distribution
        plt.figure(figsize=(15, 6))
        plt.plot(positions, gc_profile, 'b-', alpha=0.5)
        
        # Add feature annotations
        feature_colors = {
            'telomere': 'red',
            'centromere': 'green',
            'TAD': 'gray',
            'chromatin_region': 'blue'
        }
        
        y_min, y_max = min(gc_profile), max(gc_profile)
        for feature in features:
            if feature['type'] in feature_colors:
                plt.axvspan(
                    feature['start'],
                    feature['end'],
                    color=feature_colors[feature['type']],
                    alpha=0.2
                )
        
        plt.title('GC Content Distribution Along Chromosome')
        plt.xlabel('Position')
        plt.ylabel('GC Content')
        plt.grid(True, alpha=0.3)
        plt.show()

if __name__ == '__main__':
    # Run tests
    unittest.main(argv=['first-arg-is-ignored'], exit=False)