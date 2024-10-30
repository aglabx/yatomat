import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import logging
from collections import Counter, defaultdict
from typing import List, Dict, Tuple

# Add path to modules
sys.path.append(str(Path(__file__).parent.parent))

from regions.pericentromeres import (
    PericentromereRegion, PericentromereParams, PericentromereZone,
    SatelliteDistributionParams, MobileElementParams
)
from regions.common import RegionBuilder, ChromosomeRegionType, GradientGenerator
from repeats import RepeatType

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TestPericentromereRegion(unittest.TestCase):
    """Tests for pericentromeric region generation"""

    def setUp(self):
        """Initialize test parameters"""
        # Basic region parameters
        self.region_params = RegionBuilder(ChromosomeRegionType.PERICENTROMERE)\
            .set_boundaries(0, 5_000_000)\
            .set_gc_content(0.54)\
            .build()

        # Satellite distribution parameters
        self.satellite_params = SatelliteDistributionParams(
            alpha_satellite_density=0.3,
            beta_satellite_density=0.2,
            gamma_satellite_density=0.15,
            satellite_mutation_rate=0.1,
            min_block_size=1000,
            max_block_size=10000
        )

        # Mobile element parameters
        self.mobile_params = MobileElementParams(
            density=0.2,
            min_size=500,
            max_size=5000,
            types=['LINE', 'SINE', 'LTR', 'DNA']
        )

        # Complete pericentromere parameters
        self.pericentromere_params = PericentromereParams(
            min_total_length=2_000_000,
            max_total_length=5_000_000,
            transition_length=100_000,
            satellite_params=self.satellite_params,
            mobile_element_params=self.mobile_params,
            proximal_gc=0.56,
            distal_gc=0.52,
            gc_std=0.03
        )

        # Create pericentromere region instance
        self.region = PericentromereRegion(
            self.region_params,
            self.pericentromere_params
        )

    def test_sequence_generation_basic(self):
        """Test basic sequence generation properties"""
        sequence, features = self.region.generate()

        # Test sequence generation
        self.assertIsNotNone(sequence)
        self.assertGreater(len(sequence), 0)
        self.assertLessEqual(
            len(sequence),
            self.pericentromere_params.max_total_length
        )
        self.assertGreaterEqual(
            len(sequence),
            self.pericentromere_params.min_total_length
        )

        # Test that sequence contains only valid nucleotides
        valid_nucleotides = set('ATGC')
        self.assertTrue(all(n in valid_nucleotides for n in sequence))

        # Test feature generation
        self.assertIsNotNone(features)
        self.assertGreater(len(features), 0)

        logger.info(f"Generated sequence length: {len(sequence)}")
        logger.info(f"Number of features: {len(features)}")

    def test_zone_structure(self):
        """Test the structure and organization of different zones"""
        sequence, features = self.region.generate()

        # Group features by zone
        zone_features = defaultdict(list)
        for feature in features:
            zone = feature['zone']
            zone_features[zone].append(feature)

        # Test presence of all zones
        expected_zones = {
            PericentromereZone.PROXIMAL.value,
            PericentromereZone.INTERMEDIATE.value,
            PericentromereZone.DISTAL.value,
            PericentromereZone.TRANSITION.value
        }
        self.assertEqual(
            set(zone_features.keys()),
            expected_zones,
            "Not all expected zones are present"
        )

        # Test zone lengths
        zone_lengths = {}
        for zone, feats in zone_features.items():
            zone_lengths[zone] = sum(f['end'] - f['start'] for f in feats)

        logger.info("Zone lengths:")
        for zone, length in zone_lengths.items():
            logger.info(f"{zone}: {length:,} bp")

        # Test transition zone length
        transition_length = sum(
            f['end'] - f['start']
            for f in features
            if f['zone'] == PericentromereZone.TRANSITION.value
        ) / 2  # Divided by 2 because we have two transition zones

        self.assertAlmostEqual(
            transition_length,
            self.pericentromere_params.transition_length,
            delta=1000,  # Allow 1kb variation
            msg="Incorrect transition zone length"
        )

    def test_satellite_distribution(self):
        """Test satellite repeat distribution patterns"""
        sequence, features = self.region.generate()

        # Group satellites by zone and type
        zone_satellites = defaultdict(lambda: defaultdict(int))
        for feature in features:
            if feature['type'] == 'satellite':
                zone = feature['zone']
                sat_type = feature['satellite_type']
                length = feature['end'] - feature['start']
                zone_satellites[zone][sat_type] += length

        # Test satellite distribution in different zones
        for zone in zone_satellites:
            total_length = sum(zone_satellites[zone].values())
            if total_length > 0:
                distribution = {
                    sat_type: length / total_length
                    for sat_type, length in zone_satellites[zone].items()
                }
                logger.info(f"\nSatellite distribution in {zone}:")
                for sat_type, prop in distribution.items():
                    logger.info(f"{sat_type}: {prop:.3f}")

                # Test if proximal zone has more alpha satellites
                if zone == PericentromereZone.PROXIMAL.value:
                    alpha_prop = distribution.get(RepeatType.ALPHA_SATELLITE.value, 0)
                    other_props = [
                        p for t, p in distribution.items()
                        if t != RepeatType.ALPHA_SATELLITE.value
                    ]
                    if other_props:  # Only test if we have other satellite types
                        self.assertGreater(
                            alpha_prop,
                            max(other_props),
                            "Alpha satellites should dominate in proximal zone"
                        )

    def test_mobile_elements(self):
        """Test mobile element properties and distribution"""
        sequence, features = self.region.generate()

        # Collect mobile element features
        mobile_elements = [
            f for f in features
            if f['type'] == 'mobile_element'
        ]

        self.assertGreater(len(mobile_elements), 0, "No mobile elements found")

        # Test mobile element properties
        for element in mobile_elements:
            # Check size constraints
            size = element['end'] - element['start']
            self.assertGreaterEqual(
                size,
                self.mobile_params.min_size,
                "Mobile element too small"
            )
            self.assertLessEqual(
                size,
                self.mobile_params.max_size,
                "Mobile element too large"
            )

            # Check element type
            self.assertIn(
                element['element_type'],
                self.mobile_params.types,
                "Unknown mobile element type"
            )

        # Calculate overall mobile element density
        total_length = len(sequence)
        mobile_length = sum(
            element['end'] - element['start']
            for element in mobile_elements
        )
        density = mobile_length / total_length

        self.assertAlmostEqual(
            density,
            self.mobile_params.density,
            delta=0.1,
            msg="Mobile element density outside acceptable range"
        )

        logger.info(f"Mobile element density: {density:.3f}")

    def test_gc_content_gradient(self):
        """Test GC content variation across the region"""
        sequence, features = self.region.generate()

        # Calculate GC content in windows
        window_size = 10000
        gc_profile = []

        for i in range(0, len(sequence), window_size):
            window = sequence[i:i + window_size]
            gc_count = window.count('G') + window.count('C')
            gc_profile.append(gc_count / len(window))

        # Test GC content ranges
        proximal_gc = np.mean(gc_profile[:10])  # First 10 windows
        distal_gc = np.mean(gc_profile[-10:])   # Last 10 windows

        self.assertAlmostEqual(
            proximal_gc,
            self.pericentromere_params.proximal_gc,
            delta=0.05,
            msg="Incorrect proximal GC content"
        )

        self.assertAlmostEqual(
            distal_gc,
            self.pericentromere_params.distal_gc,
            delta=0.05,
            msg="Incorrect distal GC content"
        )

        # Visualize GC content profile
        self.visualize_gc_profile(gc_profile)

    @staticmethod
    def visualize_gc_profile(gc_profile: List[float]):
        """Helper method to visualize GC content distribution"""
        plt.figure(figsize=(12, 6))
        plt.plot(gc_profile)
        plt.axhline(y=np.mean(gc_profile), color='r', linestyle='--',
                   label=f'Mean GC: {np.mean(gc_profile):.3f}')
        plt.title('GC Content Profile Across Pericentromeric Region')
        plt.xlabel('Window Number (10kb windows)')
        plt.ylabel('GC Content')
        plt.legend()
        plt.grid(True)
        plt.show()

def visualize_region_structure(sequence: str, features: List[Dict]):
    """Visualize the structure of the pericentromeric region"""
    plt.figure(figsize=(15, 8))

    # Create color map for different feature types and zones
    feature_colors = {
        'satellite': 'blue',
        'mobile_element': 'red'
    }

    zone_colors = {
        PericentromereZone.PROXIMAL.value: 'lightblue',
        PericentromereZone.INTERMEDIATE.value: 'lightgreen',
        PericentromereZone.DISTAL.value: 'lightyellow',
        PericentromereZone.TRANSITION.value: 'lightgray'
    }

    # Plot features
    for feature in features:
        start = feature['start']
        end = feature['end']
        feature_type = feature['type']
        zone = feature['zone']

        # Plot zone background
        plt.axvspan(
            start, end,
            alpha=0.2,
            color=zone_colors.get(zone, 'white')
        )

        # Plot feature
        plt.axvspan(
            start, end,
            alpha=0.5,
            color=feature_colors.get(feature_type, 'gray'),
            label=f"{feature_type} ({zone})"
        )

    # Remove duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='center left',
              bbox_to_anchor=(1, 0.5))

    plt.title('Pericentromeric Region Structure')
    plt.xlabel('Position (bp)')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # Create and visualize a sample pericentromeric region
    params = RegionBuilder(ChromosomeRegionType.PERICENTROMERE)\
        .set_boundaries(0, 5_000_000)\
        .set_gc_content(0.54)\
        .build()

    pericentromere = PericentromereRegion(params)
    sequence, features = pericentromere.generate()

    # Visualize the region structure
    visualize_region_structure(sequence, features)

    # Run tests
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
