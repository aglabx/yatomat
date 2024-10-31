import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import logging
from collections import Counter, defaultdict
from typing import List, Dict, Tuple, Optional
from dataclasses import asdict

sys.path.append(str(Path(__file__).parent.parent))

from yatomat.regions.centromeres import (
    CentromereRegion, CentromereParams, CentromereZone, HORParams, CENPBParams
)
from yatomat.regions.common import RegionBuilder, ChromosomeRegionType, GradientGenerator


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TestCentromereRegion(unittest.TestCase):
    """Tests for centromere region generation"""

    def setUp(self):
        """Initialize test parameters"""
        # Basic region parameters
        self.region_params = RegionBuilder(ChromosomeRegionType.CENTROMERE)\
            .set_boundaries(0, 3_000_000)\
            .set_gc_content(0.58)\
            .build()

        # Core HOR parameters
        self.core_hor = HORParams(
            unit_size=171,       # Alpha satellite monomer size
            units_in_hor=12,     # Number of monomers in HOR
            hor_copies=100,      # Number of HOR copies
            mutation_rate_between_units=0.02,
            mutation_rate_between_hors=0.01,
            gc_content=0.58
        )

        # Peripheral HOR parameters
        self.peripheral_hor = HORParams(
            unit_size=171,
            units_in_hor=8,
            hor_copies=50,
            mutation_rate_between_units=0.05,
            mutation_rate_between_hors=0.03,
            gc_content=0.56
        )

        # CENP-B box parameters
        self.cenpb = CENPBParams()

        # Complete centromere parameters
        self.centromere_params = CentromereParams(
            min_core_length=1_500_000,
            max_core_length=3_000_000,
            min_peripheral_length=500_000,
            max_peripheral_length=1_000_000,
            transition_length=200_000,
            core_hor_params=self.core_hor,
            peripheral_hor_params=self.peripheral_hor,
            cenpb_params=self.cenpb
        )

        # Create centromere region instance
        self.region = CentromereRegion(self.region_params, self.centromere_params)

    def verify_feature(self, feature: Dict) -> bool:
        """Verify if a feature has all required fields"""
        required = {'type', 'start', 'end', 'zone'}
        return all(field in feature for field in required)

    def extract_zone_features(self, features: List[Dict]) -> Dict[CentromereZone, List[Dict]]:
        """Group features by zone"""
        zone_features = defaultdict(list)
        for feature in features:
            if self.verify_feature(feature):
                zone = CentromereZone(feature['zone'])
                zone_features[zone].append(feature)
        return dict(zone_features)

    def test_sequence_generation_basic(self):
        """Test basic sequence generation properties"""
        sequence, features = self.region.generate()

        # Test that sequence was generated
        self.assertIsNotNone(sequence)
        self.assertGreater(len(sequence), 0)

        # Test that sequence contains only valid nucleotides
        valid_nucleotides = set('ATGC')
        self.assertTrue(all(n in valid_nucleotides for n in sequence))

        # Test that features were generated
        self.assertIsNotNone(features)
        self.assertGreater(len(features), 0)

    def test_zone_lengths(self):
        """Test the lengths of different centromeric zones"""
        sequence, features = self.region.generate()

        # Extract region features
        region_features = [f for f in features if f['type'].endswith('_region')]
        zone_lengths = defaultdict(int)

        # Группируем feature по их зонам
        for feature in region_features:
            zone = CentromereZone(feature['zone'])
            length = feature['end'] - feature['start']
            zone_lengths[zone] += length

        logger.info(f"Zone lengths: {zone_lengths}")

        # Test core region length
        self.assertGreaterEqual(
            zone_lengths[CentromereZone.CORE],
            self.centromere_params.min_core_length,
            "Core region too short"
        )
        self.assertLessEqual(
            zone_lengths[CentromereZone.CORE],
            self.centromere_params.max_core_length,
            "Core region too long"
        )

        # Test individual peripheral region lengths
        peripheral_features = [
            f for f in region_features
            if CentromereZone(f['zone']) == CentromereZone.PERIPHERAL
        ]
        for feature in peripheral_features:
            length = feature['end'] - feature['start']
            self.assertGreaterEqual(
                length,
                self.centromere_params.min_peripheral_length,
                "Individual peripheral region too short"
            )
            self.assertLessEqual(
                length,
                self.centromere_params.max_peripheral_length,
                "Individual peripheral region too long"
            )

        # Test transition region lengths
        transition_features = [
            f for f in region_features
            if CentromereZone(f['zone']) == CentromereZone.TRANSITION
        ]
        for feature in transition_features:
            length = feature['end'] - feature['start']
            self.assertEqual(
                length,
                self.centromere_params.transition_length,
                "Incorrect transition region length"
            )

        # Test total sequence length matches sum of regions
        total_length = sum(zone_lengths.values())
        self.assertEqual(len(sequence), total_length)

        logger.info(f"Sequence length: {len(sequence)}, Total from regions: {total_length}")


    def test_cenp_b_boxes(self):
        """Test CENP-B box properties and distribution"""
        sequence, features = self.region.generate()

        # Get all CENP-B boxes
        cenpb_boxes = [
            f for f in features
            if f['type'] == 'CENP-B_box' and 'sequence' in f
        ]

        self.assertGreater(len(cenpb_boxes), 0, "No CENP-B boxes found")

        # Test box properties
        for box in cenpb_boxes:
            # Check sequence length
            self.assertEqual(len(box['sequence']), 16)

            # Check conserved motifs
            sequence = box['sequence']
            self.assertTrue(
                sequence.startswith('TTCG') or 'GGGA' in sequence,
                f"Invalid CENP-B box sequence: {sequence}"
            )

        # Test spacing
        box_positions = sorted(box['start'] for box in cenpb_boxes)
        spacings = np.diff(box_positions)
        median_spacing = np.median(spacings)

        self.assertAlmostEqual(
            median_spacing, 171, delta=20,
            msg="Incorrect CENP-B box spacing"
        )

    def test_cenp_b_boxes_in_core_only(self):
        """Test that CENP-B boxes are only present in the core region"""
        sequence, features = self.region.generate()

        # Get all CENP-B boxes
        cenpb_boxes = [
            f for f in features
            if f['type'] == 'CENP-B_box' and 'sequence' in f
        ]

        self.assertGreater(len(cenpb_boxes), 0, "No CENP-B boxes found")

        # Check that all CENP-B boxes are in the core region
        for box in cenpb_boxes:
            zone = box['zone']
            self.assertEqual(zone, CentromereZone.CORE.value, f"CENP-B box found in non-core zone: {zone}")

    def test_mutations(self):
        """Test mutation patterns in different zones"""
        _, features = self.region.generate()
        zone_features = self.extract_zone_features(features)

        # Calculate mutation rates for each zone
        mutation_rates = {}
        for zone, feats in zone_features.items():
            mutations = [f for f in feats if f['type'] == 'substitution']
            mutation_rates[zone] = len(mutations) / len(feats) if feats else 0

        # Core should have lower mutation rate than peripheral
        self.assertLess(
            mutation_rates.get(CentromereZone.CORE, 0),
            mutation_rates.get(CentromereZone.PERIPHERAL, 0),
            "Core mutation rate should be lower than peripheral"
        )

    def test_transitions(self):
        """Test properties of transition zones"""
        _, features = self.region.generate()

        # Get transition features with gradient values
        transition_features = sorted(
            [f for f in features if (
                f['zone'] == CentromereZone.TRANSITION.value and
                'gradient_value' in f and
                'start' in f
            )],
            key=lambda x: x['start']
        )

        self.assertGreater(len(transition_features), 0, "No transition features found")

        # Assuming we have two transition zones, split the features
        half = len(transition_features) // 2
        left_transition = transition_features[:half]
        right_transition = transition_features[half:]

        # Test left transition gradient
        left_gradient_values = [f['gradient_value'] for f in left_transition]
        self.assertLess(left_gradient_values[0], 0.1, "Left transition initial gradient value too high")
        self.assertGreater(left_gradient_values[-1], 0.9, "Left transition final gradient value too low")
        # Test monotonicity
        self.assertTrue(
            all(a <= b for a, b in zip(left_gradient_values, left_gradient_values[1:])),
            "Left transition gradient is not monotonically increasing"
        )

        # Test right transition gradient
        right_gradient_values = [f['gradient_value'] for f in right_transition]
        self.assertGreater(right_gradient_values[0], 0.9, "Right transition initial gradient value too low")
        self.assertLess(right_gradient_values[-1], 0.1, "Right transition final gradient value too high")
        # Test monotonicity
        self.assertTrue(
            all(a >= b for a, b in zip(right_gradient_values, right_gradient_values[1:])),
            "Right transition gradient is not monotonically decreasing"
        )

    def test_visualize_centromere(self):
        """Test visualization of the centromere structure"""
        sequence, features = self.region.generate()
        visualize_centromere(sequence, features)
        # Since this test is for visualization, we do not have assertions here.
        # It will display the plot when the test is run.

def visualize_centromere(sequence: str, features: List[Dict]):
    """Visualize centromere structure"""
    plt.figure(figsize=(15, 8))
    colors = {
        CentromereZone.CORE.value: 'blue',
        CentromereZone.PERIPHERAL.value: 'green',
        CentromereZone.TRANSITION.value: 'orange'
    }

    # Plot zones
    zone_labels = set()
    for feature in features:
        if feature['type'].endswith('_region'):
            zone = feature['zone']
            color = colors.get(zone, 'gray')
            label = zone if zone not in zone_labels else None
            zone_labels.add(zone)
            plt.axvspan(
                feature['start'],
                feature['end'],
                alpha=0.3,
                color=color,
                label=label
            )

    # Plot CENP-B boxes
    boxes = [f for f in features if f['type'] == 'CENP-B_box']
    if boxes:
        positions = [box['start'] for box in boxes]
        plt.vlines(positions, 0, 1, color='red', alpha=0.5, label='CENP-B box')

    plt.title('Centromere Structure')
    plt.xlabel('Position (bp)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
