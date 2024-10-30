import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import logging
from collections import Counter, defaultdict
from typing import List, Dict, Tuple
from scipy.stats import gaussian_kde

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
    
    def test_satellite_distribution_detailed(self):
        """Test detailed satellite repeat distribution patterns"""
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
        
        # Test block size distribution with zone-specific constraints
        block_sizes = defaultdict(list)
        for feature in features:
            if feature['type'] == 'satellite':
                zone = feature['zone']
                size = feature['end'] - feature['start']
                block_sizes[zone].append(size)
        
        for zone, sizes in block_sizes.items():
            if sizes:
                mean_size = np.mean(sizes)
                std_size = np.std(sizes)
                logger.info(f"\nSatellite block size stats for {zone}:")
                logger.info(f"Mean size: {mean_size:.0f}")
                logger.info(f"Std size: {std_size:.0f}")
                
                # Define zone-specific size constraints
                if zone == PericentromereZone.PROXIMAL.value:
                    min_size = self.satellite_params.min_block_size
                    max_size = self.satellite_params.max_block_size
                elif zone == PericentromereZone.DISTAL.value:
                    # Allow smaller blocks in distal regions
                    min_size = self.satellite_params.min_block_size // 2
                    max_size = self.satellite_params.max_block_size * 2
                else:
                    # Intermediate and transition zones use default constraints
                    min_size = self.satellite_params.min_block_size
                    max_size = self.satellite_params.max_block_size
                
                # Test size constraints with some tolerance for edge cases
                size_violations = [s for s in sizes if s < min_size or s > max_size]
                violation_rate = len(size_violations) / len(sizes)
                
                self.assertLess(
                    violation_rate,
                    0.1,  # Allow up to 10% violations
                    f"Too many block size violations in {zone} zone"
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
    
    def test_spatial_organization(self):
        """Test spatial organization and clustering of features"""
        sequence, features = self.region.generate()
        
        # Calculate nearest neighbor distances for different feature types
        def calculate_distances(features_subset):
            positions = [(f['start'] + f['end']) // 2 for f in features_subset]
            distances = []
            for i in range(len(positions) - 1):
                distances.append(positions[i + 1] - positions[i])
            return distances if distances else [0]
        
        # Analyze distributions for different feature types
        feature_distances = defaultdict(list)
        for feature_type in ['satellite', 'mobile_element']:
            type_features = sorted(
                [f for f in features if f['type'] == feature_type],
                key=lambda x: x['start']
            )
            feature_distances[feature_type] = calculate_distances(type_features)
        
        for feat_type, distances in feature_distances.items():
            if distances:
                mean_dist = np.mean(distances)
                std_dist = np.std(distances)
                logger.info(f"\nSpatial distribution for {feat_type}:")
                logger.info(f"Mean distance: {mean_dist:.0f}")
                logger.info(f"Std distance: {std_dist:.0f}")
                
                # Test for clustering (non-random distribution)
                # Calculate coefficient of variation (CV)
                cv = std_dist / mean_dist if mean_dist > 0 else 0
                logger.info(f"Coefficient of variation: {cv:.3f}")
                
                # CV > 1 indicates clustering
                self.assertGreater(
                    cv, 0.5,
                    f"{feat_type} features should show some clustering"
                )
        
        # Visualize spatial distribution
        self.visualize_spatial_distribution(features)
    
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
    def visualize_gc_profile(gc_profile: List[float], zones: Optional[List[str]] = None):
        """Enhanced helper method to visualize GC content distribution with zones"""
        plt.figure(figsize=(12, 6))
        plt.figure(figsize=(15, 8))
        
        # Plot GC profile
        x = np.arange(len(gc_profile))
        plt.plot(x, gc_profile, 'b-', alpha=0.7, label='GC content')
        
        # Add rolling average
        window = 5
        rolling_mean = np.convolve(gc_profile, np.ones(window)/window, mode='valid')
        plt.plot(x[window-1:], rolling_mean, 'r-', linewidth=2, label='Rolling average')
        
        # Add mean line
        mean_gc = np.mean(gc_profile)
        plt.axhline(y=mean_gc, color='k', linestyle='--',
                   label=f'Mean GC: {mean_gc:.3f}')
        
        # Add zone backgrounds if provided
        if zones:
            zone_boundaries = np.linspace(0, len(gc_profile), len(zones) + 1)
            zone_colors = {'proximal': 'lightblue',
                         'intermediate': 'lightgreen',
                         'distal': 'lightyellow',
                         'transition': 'lightgray'}
            
            for i, zone in enumerate(zones):
                plt.axvspan(zone_boundaries[i], zone_boundaries[i+1],
                          alpha=0.2, color=zone_colors.get(zone, 'white'),
                          label=f'{zone} zone')
        
        plt.title('GC Content Profile Across Pericentromeric Region', pad=20)
        plt.xlabel('Window Number (10kb windows)', labelpad=10)
        plt.ylabel('GC Content', labelpad=10)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        # Add confidence intervals
        std_gc = np.std(gc_profile)
        plt.fill_between(x, gc_profile - std_gc, gc_profile + std_gc,
                        color='blue', alpha=0.1)
        
        plt.tight_layout()
        plt.show()

    # Add these methods to the TestPericentromereRegion class

    def test_mutation_patterns(self):
        """Test detailed mutation patterns in different zones"""
        sequence, features = self.region.generate()
        
        # Collect mutation data by zone
        zone_mutations = defaultdict(list)
        for feature in features:
            if feature['type'] == 'satellite' and 'mutations' in feature:
                zone = feature['zone']
                mutations = feature.get('mutations', 0)  # Get mutation count directly
                sequence_length = feature['end'] - feature['start']
                mutation_rate = mutations / sequence_length if sequence_length > 0 else 0
                zone_mutations[zone].append({
                    'rate': mutation_rate,
                    'count': mutations,  # Store count directly
                    'length': sequence_length
                })
        
        # Analyze mutation patterns for each zone
        for zone, mutations in zone_mutations.items():
            if mutations:
                rates = [m['rate'] for m in mutations]
                mean_rate = np.mean(rates)
                std_rate = np.std(rates)
                
                logger.info(f"\nMutation analysis for {zone}:")
                logger.info(f"Mean mutation rate: {mean_rate:.4f}")
                logger.info(f"Std deviation: {std_rate:.4f}")
                
                # Test mutation rate constraints
                if zone == PericentromereZone.PROXIMAL.value:
                    self.assertLess(
                        mean_rate,
                        self.satellite_params.satellite_mutation_rate,
                        "Proximal zone mutation rate should be lower than base rate"
                    )
                elif zone == PericentromereZone.DISTAL.value:
                    self.assertGreater(
                        mean_rate,
                        self.satellite_params.satellite_mutation_rate,
                        "Distal zone mutation rate should be higher than base rate"
                    )

    def test_feature_boundaries(self):
        """Test feature boundary integrity and transitions"""
        sequence, features = self.region.generate()
        
        # Sort features by start position
        sorted_features = sorted(features, key=lambda x: x['start'])
        
        # Check for excessive overlaps
        for i in range(len(sorted_features) - 1):
            current = sorted_features[i]
            next_feat = sorted_features[i + 1]
            
            # Small overlaps are allowed for transition features
            overlap_size = current['end'] - next_feat['start']
            max_allowed_overlap = 500  # 500bp maximum allowed overlap
            
            if overlap_size > 0:  # If there is an overlap
                # Only allow small overlaps in transition zones
                if current['zone'] == PericentromereZone.TRANSITION.value or \
                next_feat['zone'] == PericentromereZone.TRANSITION.value:
                    self.assertLessEqual(
                        overlap_size,
                        max_allowed_overlap,
                        f"Excessive overlap ({overlap_size} bp) between transition features"
                    )
                else:
                    self.assertLessEqual(
                        overlap_size,
                        0,
                        "Non-transition features should not overlap"
                    )
            
            # Test for large gaps
            if overlap_size < 0:  # If there is a gap
                gap_size = abs(overlap_size)
                self.assertLess(
                    gap_size,
                    1000,  # Maximum allowed gap size (1kb)
                    f"Large gap ({gap_size} bp) found between features"
                )
            
            # Test zone transitions
            if current['zone'] != next_feat['zone']:
                # Only certain zone transitions are allowed
                allowed_transitions = {
                    PericentromereZone.PROXIMAL.value: {PericentromereZone.TRANSITION.value},
                    PericentromereZone.TRANSITION.value: {
                        PericentromereZone.INTERMEDIATE.value,
                        PericentromereZone.PROXIMAL.value,
                        PericentromereZone.DISTAL.value
                    },
                    PericentromereZone.INTERMEDIATE.value: {PericentromereZone.TRANSITION.value},
                    PericentromereZone.DISTAL.value: {PericentromereZone.TRANSITION.value}
                }
                
                self.assertIn(
                    next_feat['zone'],
                    allowed_transitions.get(current['zone'], set()),
                    f"Invalid zone transition: {current['zone']} to {next_feat['zone']}"
                )

    def test_internal_consistency(self):
        """Test internal consistency of features and their properties"""
        sequence, features = self.region.generate()
        
        for feature in features:
            # Test basic feature structure
            required_fields = {'type', 'start', 'end', 'zone'}
            self.assertTrue(
                all(field in feature for field in required_fields),
                f"Missing required fields in feature: {feature}"
            )
            
            # Test position consistency
            self.assertGreaterEqual(feature['start'], 0)
            self.assertGreater(feature['end'], feature['start'])
            self.assertLessEqual(feature['end'], len(sequence))
            
            # Test type-specific properties
            if feature['type'] == 'satellite':
                self.assertIn('satellite_type', feature)
                self.assertIn('monomer_size', feature)
                if 'mutations' in feature:
                    self.assertIsInstance(feature['mutations'], int)
                    
            elif feature['type'] == 'mobile_element':
                self.assertIn('element_type', feature)
                self.assertIn(
                    feature['element_type'],
                    self.mobile_params.types
                )
                
            # Test zone-specific constraints
            if feature['zone'] == PericentromereZone.PROXIMAL.value:
                if feature['type'] == 'satellite':
                    self.assertLessEqual(
                        feature['end'] - feature['start'],
                        self.satellite_params.max_block_size,
                        "Proximal satellite block too large"
                    )

    def visualize_spatial_distribution(self, features: List[Dict]):
        """Visualize spatial distribution of features"""
        plt.figure(figsize=(15, 8))
        
        # Get total sequence length
        total_length = max(f['end'] for f in features)
        
        # Create separate axes for different feature types
        feature_types = {'satellite', 'mobile_element'}
        positions = defaultdict(list)
        
        # Collect positions for each feature type
        for feature in features:
            if feature['type'] in feature_types:
                pos = (feature['start'] + feature['end']) / (2 * total_length)  # normalized position
                positions[feature['type']].append(pos)
        
        # Plot density for each feature type
        colors = {'satellite': 'blue', 'mobile_element': 'red'}
        for feat_type, pos_list in positions.items():
            if pos_list:
                density = gaussian_kde(pos_list)
                xs = np.linspace(0, 1, 200)
                plt.plot(xs, density(xs), label=f'{feat_type}', color=colors.get(feat_type, 'gray'))
                
                # Add scatter points
                plt.scatter(pos_list, [0.01] * len(pos_list), 
                        color=colors.get(feat_type, 'gray'),
                        alpha=0.3, marker='|')
        
        # Add zone backgrounds
        zone_ranges = {
            'proximal': (0, 0.3),
            'intermediate': (0.3, 0.7),
            'distal': (0.7, 1.0)
        }
        
        zone_colors = {
            'proximal': 'lightblue',
            'intermediate': 'lightgreen',
            'distal': 'lightyellow'
        }
        
        for zone, (start, end) in zone_ranges.items():
            plt.axvspan(start, end, alpha=0.2, 
                    color=zone_colors[zone], 
                    label=f'{zone} zone')
        
        plt.title('Spatial Distribution of Features')
        plt.xlabel('Normalized Position')
        plt.ylabel('Density')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    # Add these visualization methods to the test class

    def visualize_mutation_distribution(self, features: List[Dict]):
        """Visualize mutation distribution across zones"""
        plt.figure(figsize=(15, 8))
        
        # Collect mutation data
        zone_mutations = defaultdict(list)
        zone_positions = defaultdict(list)
        
        for feature in features:
            if feature['type'] == 'satellite' and 'mutations' in feature:
                zone = feature['zone']
                mutations = feature.get('mutations', 0)
                position = (feature['start'] + feature['end']) / 2
                sequence_length = feature['end'] - feature['start']
                mutation_rate = mutations / sequence_length if sequence_length > 0 else 0
                
                zone_mutations[zone].append(mutation_rate)
                zone_positions[zone].append(position)
        
        # Plot mutation rates for each zone
        zone_colors = {
            PericentromereZone.PROXIMAL.value: 'blue',
            PericentromereZone.INTERMEDIATE.value: 'green',
            PericentromereZone.DISTAL.value: 'red',
            PericentromereZone.TRANSITION.value: 'gray'
        }
        
        for zone in zone_mutations:
            positions = np.array(zone_positions[zone])
            rates = np.array(zone_mutations[zone])
            
            # Plot scatter points
            plt.scatter(positions, rates, c=zone_colors.get(zone, 'black'),
                    alpha=0.5, label=f'{zone}')
            
            # Add trend line if enough points
            if len(positions) > 1:
                z = np.polyfit(positions, rates, 1)
                p = np.poly1d(z)
                plt.plot(positions, p(positions), 
                        c=zone_colors.get(zone, 'black'),
                        linestyle='--', alpha=0.8)
        
        plt.title('Mutation Rate Distribution Across Zones')
        plt.xlabel('Position in Region')
        plt.ylabel('Mutation Rate')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    def visualize_feature_density(self, features: List[Dict], window_size: int = 50000):
        """Visualize feature density distribution"""
        plt.figure(figsize=(15, 10))
        
        # Get sequence length
        total_length = max(f['end'] for f in features)
        windows = np.arange(0, total_length, window_size)
        
        # Calculate density for different feature types
        feature_types = ['satellite', 'mobile_element']
        densities = {ft: np.zeros(len(windows)) for ft in feature_types}
        
        for feature in features:
            if feature['type'] in feature_types:
                # Find windows that this feature spans
                start_window = feature['start'] // window_size
                end_window = min(feature['end'] // window_size + 1, len(windows))
                feature_length = feature['end'] - feature['start']
                
                for w in range(start_window, end_window):
                    window_start = w * window_size
                    window_end = (w + 1) * window_size
                    
                    # Calculate overlap with window
                    overlap_start = max(feature['start'], window_start)
                    overlap_end = min(feature['end'], window_end)
                    overlap = max(0, overlap_end - overlap_start)
                    
                    densities[feature['type']][w] += overlap / window_size
        
        # Plot densities
        for ft in feature_types:
            plt.plot(windows + window_size/2, densities[ft], 
                    label=ft, alpha=0.7)
        
        # Add zone backgrounds
        zone_boundaries = {
            PericentromereZone.PROXIMAL.value: [0, total_length * 0.3],
            PericentromereZone.INTERMEDIATE.value: [total_length * 0.3, total_length * 0.7],
            PericentromereZone.DISTAL.value: [total_length * 0.7, total_length]
        }
        
        zone_colors = {
            PericentromereZone.PROXIMAL.value: 'lightblue',
            PericentromereZone.INTERMEDIATE.value: 'lightgreen',
            PericentromereZone.DISTAL.value: 'lightyellow'
        }
        
        for zone, (start, end) in zone_boundaries.items():
            plt.axvspan(start, end, alpha=0.2, 
                    color=zone_colors.get(zone, 'white'),
                    label=f'{zone} zone')
        
        plt.title('Feature Density Distribution')
        plt.xlabel('Position in Region')
        plt.ylabel('Feature Density')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

def visualize_region_structure(sequence: str, features: List[Dict], detailed: bool = True):
    """Visualize the structure of the pericentromeric region"""
    if detailed:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12), height_ratios=[3, 1])
    else:
        fig, ax1 = plt.subplots(figsize=(15, 8))
        ax2 = None
    
    # Create color map for different feature types and zones
    feature_colors = {
        'satellite': {
            RepeatType.ALPHA_SATELLITE.value: 'blue',
            RepeatType.BETA_SATELLITE.value: 'green',
            RepeatType.SATELLITE_1.value: 'purple'
        },
        'mobile_element': {
            'LINE': 'red',
            'SINE': 'orange',
            'LTR': 'pink',
            'DNA': 'brown'
        }
    }
    
    zone_colors = {
        PericentromereZone.PROXIMAL.value: 'lightblue',
        PericentromereZone.INTERMEDIATE.value: 'lightgreen',
        PericentromereZone.DISTAL.value: 'lightyellow',
        PericentromereZone.TRANSITION.value: 'lightgray'
    }
    
    # Plot features in main axis
    total_length = max(f['end'] for f in features)
    
    # Plot zone backgrounds first
    for zone in zone_colors:
        zone_features = [f for f in features if f['zone'] == zone]
        if zone_features:
            for feature in zone_features:
                ax1.axvspan(
                    feature['start'] / total_length,
                    feature['end'] / total_length,
                    alpha=0.2,
                    color=zone_colors[zone]
                )
    
    # Plot features
    for feature in features:
        start = feature['start'] / total_length
        end = feature['end'] / total_length
        feature_type = feature['type']
        
        if feature_type == 'satellite':
            sat_type = feature.get('satellite_type', 'unknown')
            color = feature_colors['satellite'].get(sat_type, 'gray')
            label = f"{sat_type}"
        elif feature_type == 'mobile_element':
            element_type = feature.get('element_type', 'unknown')
            color = feature_colors['mobile_element'].get(element_type, 'gray')
            label = f"{element_type}"
        else:
            color = 'gray'
            label = feature_type
        
        ax1.axvspan(start, end, alpha=0.5, color=color, label=label)
        
        # Add feature details text if zoomed in enough
        if end - start > 0.02:
            ax1.text((start + end) / 2, 0.5, label,
                    horizontalalignment='center',
                    verticalalignment='center',
                    rotation=90)
    
    if detailed and ax2 is not None:
        # Create density plot in lower axis
        positions = []
        types = []
        for feature in features:
            pos = (feature['start'] + feature['end']) / (2 * total_length)
            positions.append(pos)
            types.append(feature['type'])
        
        # Calculate density for different feature types
        for feat_type in ['satellite', 'mobile_element']:
            type_positions = [pos for pos, t in zip(positions, types) if t == feat_type]
            if type_positions:
                density = gaussian_kde(type_positions)
                xs = np.linspace(0, 1, 200)
                ax2.plot(xs, density(xs), label=f'{feat_type} density')
    
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