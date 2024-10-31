import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import logging
from collections import Counter, defaultdict
from typing import List, Dict, Tuple

from regions.boundaries import (
    BoundaryRegion, BoundaryParams, ChromatinState,
    TADParams, ChromatinTransitionParams, CTCFSite
)
from regions.common import (
    RegionBuilder, ChromosomeRegionType, GradientGenerator
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TestBoundaryRegion(unittest.TestCase):
    """Tests for boundary region generation"""
    
    def setUp(self):
        """Initialize test parameters"""
        # Basic region parameters
        self.region_params = RegionBuilder(ChromosomeRegionType.TRANSITION)\
            .set_boundaries(0, 200_000)\
            .set_gc_content(0.45)\
            .build()

        # TAD parameters
        self.tad_params = TADParams(
            min_size=100_000,
            max_size=1_000_000,
            boundary_strength=0.8,
            ctcf_density=0.01,
            internal_structure_complexity=0.5
        )

        # Chromatin transition parameters
        self.transition_params = ChromatinTransitionParams(
            transition_length=50_000,
            gc_content_range=(0.35, 0.65),
            repeat_density_range=(0.1, 0.8),
            nucleosome_positioning_strength=0.7,
            gradient_steepness=2.0
        )

        # Complete boundary parameters
        self.boundary_params = BoundaryParams(
            tad_params=self.tad_params,
            transition_params=self.transition_params,
            insulator_strength=0.8,
            min_boundary_size=5_000,
            max_boundary_size=20_000
        )

        # Create boundary region instance
        self.region = BoundaryRegion(self.region_params, self.boundary_params)

    def test_sequence_generation_basic(self):
        """Test basic sequence generation properties"""
        sequence, features = self.region.generate()

        # Test sequence generation
        self.assertIsNotNone(sequence)
        self.assertGreater(len(sequence), 0)
        self.assertEqual(len(sequence), self.region_params.end - self.region_params.start)

        # Test that sequence contains only valid nucleotides
        valid_nucleotides = set('ATGCN')  # Include 'N' if necessary
        self.assertTrue(all(n in valid_nucleotides for n in sequence))

        # Test feature generation
        self.assertIsNotNone(features)
        self.assertGreater(len(features), 0)

        logger.info(f"Generated sequence length: {len(sequence)}")
        logger.info(f"Number of features: {len(features)}")

    def test_tad_properties(self):
        """Test properties of Topologically Associating Domains"""
        sequence, features = self.region.generate()

        # Get TAD features
        tad_features = [f for f in features if f['type'] == 'TAD']
        self.assertGreater(len(tad_features), 0, "No TAD features found")

        for tad in tad_features:
            # Test TAD size constraints
            size = tad['end'] - tad['start']
            self.assertGreaterEqual(
                size,
                self.tad_params.min_size,
                "TAD too small"
            )
            self.assertLessEqual(
                size,
                self.tad_params.max_size,
                "TAD too large"
            )

            # Test CTCF site presence
            ctcf_sites = [
                f for f in features
                if f['type'] == 'CTCF_site'
                and f['start'] >= tad['start']
                and f['end'] <= tad['end']
            ]

            self.assertGreater(len(ctcf_sites), 0, "No CTCF sites in TAD")

            # Check CTCF orientation at boundaries
            boundary_sites = [
                s for s in ctcf_sites
                if abs(s['start'] - tad['start']) < 100
                or abs(s['end'] - tad['end']) < 100
            ]
            self.assertGreaterEqual(len(boundary_sites), 2, "Missing boundary CTCF sites")

    def test_ctcf_site_properties(self):
        """Test properties of CTCF binding sites"""
        sequence, features = self.region.generate()

        ctcf_sites = [f for f in features if f['type'] == 'CTCF_site']
        self.assertGreater(len(ctcf_sites), 0, "No CTCF sites found")

        for site in ctcf_sites:
            # Test site structure
            self.assertIn('orientation', site)
            self.assertIn('strength', site)
            self.assertGreaterEqual(site['strength'], 0)
            self.assertLessEqual(site['strength'], 1)

            # Check orientation
            self.assertIn(site['orientation'], ['forward', 'reverse'])

            # Check spacing between sites
            nearby_sites = [
                s for s in ctcf_sites
                if s != site and abs(s['start'] - site['start']) < 1000
            ]
            self.assertLess(
                len(nearby_sites), 3,
                "Too many CTCF sites in close proximity"
            )

    def test_chromatin_transitions(self):
        """Test chromatin state transition properties"""
        sequence, features = self.region.generate()

        transition_features = [
            f for f in features
            if f['type'] == 'chromatin_transition'
        ]
        self.assertGreater(len(transition_features), 0)

        for transition in transition_features:
            # Test transition length
            length = transition['end'] - transition['start']
            self.assertEqual(
                length,
                self.transition_params.transition_length,
                "Incorrect transition length"
            )

            # Test GC content range
            gc_range = transition['gc_range']
            self.assertGreaterEqual(gc_range[0], self.transition_params.gc_content_range[0])
            self.assertLessEqual(gc_range[1], self.transition_params.gc_content_range[1])

            # Test repeat density range
            repeat_range = transition['repeat_range']
            self.assertGreaterEqual(
                repeat_range[0],
                self.transition_params.repeat_density_range[0]
            )
            self.assertLessEqual(
                repeat_range[1],
                self.transition_params.repeat_density_range[1]
            )

            # Test state transitions
            self.assertIn('start_state', transition)
            self.assertIn('end_state', transition)
            self.assertIsInstance(
                ChromatinState(transition['start_state']),
                ChromatinState
            )
            self.assertIsInstance(
                ChromatinState(transition['end_state']),
                ChromatinState
            )

    def test_gradient_properties(self):
        """Test gradient-based properties across the region"""
        sequence, features = self.region.generate()

        # Calculate GC content in windows
        window_size = 1000
        gc_profile = []
        for i in range(0, len(sequence), window_size):
            window = sequence[i:i + window_size]
            gc_count = window.count('G') + window.count('C')
            gc_profile.append(gc_count / len(window))

        # Test GC content gradient in transition regions
        transition_features = [
            f for f in features
            if f['type'] == 'chromatin_transition'
        ]

        for transition in transition_features:
            start_idx = transition['start'] // window_size
            end_idx = transition['end'] // window_size
            transition_gc = gc_profile[start_idx:end_idx]

            # Test monotonicity of gradient
            differences = np.diff(transition_gc)
            is_monotonic = all(diff >= 0 for diff in differences) or all(diff <= 0 for diff in differences)
            self.assertTrue(
                is_monotonic,
                "GC content gradient is not monotonic in transition region"
            )

    def test_feature_integrity(self):
        """Test integrity and consistency of generated features"""
        sequence, features = self.region.generate()

        # Sort features by start position
        sorted_features = sorted(features, key=lambda x: x['start'])

        # Check feature positions and overlaps
        for i in range(len(sorted_features) - 1):
            current = sorted_features[i]
            next_feat = sorted_features[i + 1]

            # Basic position checks
            self.assertGreaterEqual(current['start'], 0)
            self.assertGreater(current['end'], current['start'])
            self.assertLessEqual(current['end'], len(sequence))

            # Check for invalid overlaps
            if current['type'] != 'chromatin_transition' and next_feat['type'] != 'chromatin_transition':
                self.assertLessEqual(
                    current['end'],
                    next_feat['start'],
                    "Invalid feature overlap"
                )

    def visualize_boundary_structure(self):
        """Visualize the structure of the boundary region"""
        sequence, features = self.region.generate()
        
        plt.figure(figsize=(15, 10))
        
        # Create color map for different feature types
        color_map = {
            'TAD': 'blue',
            'CTCF_site': 'red',
            'chromatin_transition': 'green',
            'repeat': 'orange'
        }
        
        # Plot features
        for feature in features:
            start = feature['start']
            end = feature['end']
            feature_type = feature['type']
            
            if feature_type == 'CTCF_site':
                # Plot CTCF sites as vertical lines
                plt.axvline(
                    x=start,
                    color=color_map[feature_type],
                    alpha=0.5,
                    linewidth=2
                )
            else:
                # Plot other features as spans
                plt.axvspan(
                    start,
                    end,
                    alpha=0.3,
                    color=color_map.get(feature_type, 'gray'),
                    label=feature_type
                )
        
        # Remove duplicate labels
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
        
        plt.title('Boundary Region Structure')
        plt.xlabel('Position (bp)')
        plt.ylabel('Feature Type')
        plt.grid(True, alpha=0.3)
        plt.show()

    def test_visualize_gradients(self):
        """Test and visualize gradients in the boundary region"""
        sequence, features = self.region.generate()

        # Extract transition features
        transition_features = [
            f for f in features
            if f['type'] == 'chromatin_transition'
        ]

        if not transition_features:
            return

        plt.figure(figsize=(15, 10))

        # Plot GC content gradient
        window_size = 1000
        positions = []
        gc_values = []

        for i in range(0, len(sequence), window_size):
            window = sequence[i:i + window_size]
            gc_count = window.count('G') + window.count('C')
            gc_content = gc_count / len(window)
            positions.append(i + window_size/2)
            gc_values.append(gc_content)

        plt.subplot(211)
        plt.plot(positions, gc_values, label='GC content')
        plt.title('GC Content Gradient')
        plt.xlabel('Position (bp)')
        plt.ylabel('GC Content')
        plt.grid(True, alpha=0.3)
        plt.legend()

        # Plot feature density gradient
        plt.subplot(212)
        feature_density = np.zeros(len(positions))
        
        for feature in features:
            if feature['type'] in ['repeat', 'CTCF_site']:
                start_idx = feature['start'] // window_size
                end_idx = feature['end'] // window_size
                if start_idx < len(feature_density):
                    feature_density[start_idx:end_idx+1] += 1

        plt.plot(positions, feature_density, label='Feature density')
        plt.title('Feature Density Gradient')
        plt.xlabel('Position (bp)')
        plt.ylabel('Feature Count')
        plt.grid(True, alpha=0.3)
        plt.legend()

        plt.tight_layout()
        plt.show()

if __name__ == '__main__':
    # Create instance for visualization
    params = RegionBuilder(ChromosomeRegionType.TRANSITION)\
        .set_boundaries(0, 200_000)\
        .set_gc_content(0.45)\
        .build()
    
    boundary = BoundaryRegion(params)
    
    # Run visualization
    test = TestBoundaryRegion()
    test.setUp()
    test.visualize_boundary_structure()
    test.test_visualize_gradients()
    
    # Run tests
    unittest.main(argv=['first-arg-is-ignored'], exit=False)