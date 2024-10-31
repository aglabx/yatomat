import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import logging
from collections import Counter
from typing import List, Dict, Tuple

# Add path to modules
sys.path.append(str(Path(__file__).parent.parent))

from yatomat.regions.telomers import (
    TelomereRegion, SubtelomereRegion, TelomereParams, SubtelomereParams
)
from yatomat.regions.common import RegionBuilder, ChromosomeRegionType, GradientGenerator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TestTelomereRegion(unittest.TestCase):
    """Tests for the TelomereRegion class"""
    
    def setUp(self):
        """Set up test parameters"""
        self.params = RegionBuilder(ChromosomeRegionType.TELOMERE)\
            .set_boundaries(0, 10000)\
            .set_gc_content(0.48)\
            .build()
        
        self.telomere_params = TelomereParams(
            min_length=3000,
            max_length=15000,
            mutation_rate_5_end=0.01,
            mutation_rate_3_end=0.05,
            variant_repeat_prob=0.1
        )
        
        self.region = TelomereRegion(self.params, self.telomere_params)
    
    def test_length_constraints(self):
        """Test if generated telomere length is within specified bounds"""
        sequence, _ = self.region.generate()
        self.assertGreaterEqual(len(sequence), self.telomere_params.min_length)
        self.assertLessEqual(len(sequence), self.telomere_params.max_length)
    
    def test_mutation_gradient(self):
        """Test if mutations follow the expected gradient pattern"""
        sequence, features = self.region.generate()
        
        # Split sequence into windows and calculate mutation rates
        window_size = 500
        windows = [sequence[i:i+window_size] 
                  for i in range(0, len(sequence), window_size)]
        
        consensus = self.telomere_params.consensus_repeat
        mutation_rates = []
        
        for window in windows:
            # Count perfect repeats vs mutated repeats
            repeats = [window[i:i+len(consensus)] 
                      for i in range(0, len(window), len(consensus))]
            mutated = sum(1 for repeat in repeats if repeat != consensus)
            mutation_rates.append(mutated / len(repeats) if repeats else 0)
        
        # Check if mutation rate increases from 5' to 3' end
        self.assertLess(mutation_rates[0], mutation_rates[-1])
        
        # Visualize mutation gradient
        plt.figure(figsize=(10, 6))
        plt.plot(mutation_rates)
        plt.title('Telomere Mutation Rate Gradient')
        plt.xlabel('Window Number')
        plt.ylabel('Mutation Rate')
        plt.grid(True)
        plt.show()
    
    def test_variant_repeat_distribution(self):
        """Test the distribution of canonical vs variant repeats"""
        _, features = self.region.generate()
        
        repeat_types = [feature['type'] for feature in features]
        type_counts = Counter(repeat_types)
        
        # Calculate the proportion of variant repeats
        total_repeats = len(repeat_types)
        variant_proportion = type_counts.get('variant', 0) / total_repeats
        
        # Check if proportion is close to specified probability
        self.assertAlmostEqual(
            variant_proportion, 
            self.telomere_params.variant_repeat_prob,
            delta=0.05
        )
    
    def test_feature_annotations(self):
        """Test if feature annotations are correct and complete"""
        sequence, features = self.region.generate()
        sequence_length = len(sequence)
        
        # Check that features cover the entire sequence
        covered_positions = set()
        for feature in features:
            # Check feature structure
            self.assertIn('type', feature)
            self.assertIn('start', feature)
            self.assertIn('end', feature)
            self.assertIn('sequence', feature)
            
            # Check positions
            start = feature['start']
            end = min(feature['end'], sequence_length)  # Truncate to sequence length
            
            self.assertGreaterEqual(start, 0)
            self.assertLessEqual(start, sequence_length)  # Start must be within sequence
            self.assertGreater(end, start)
            
            # Add positions to covered set (only within sequence bounds)
            if start < sequence_length:
                covered_positions.update(range(start, end))
        
        # Verify complete coverage
        self.assertEqual(len(covered_positions), len(sequence))

class TestSubtelomereRegion(unittest.TestCase):
    """Tests for the SubtelomereRegion class"""
    
    def setUp(self):
        """Set up test parameters"""
        self.params = RegionBuilder(ChromosomeRegionType.SUBTELOMERE)\
            .set_boundaries(10000, 50000)\
            .set_gc_content(0.52)\
            .build()
        
        self.subtelomere_params = SubtelomereParams(
            min_length=10000,
            max_length=300000,
            satellite_density=0.4,
            transposon_density=0.2
        )
        
        self.region = SubtelomereRegion(self.params, self.subtelomere_params)
    
    def test_length_constraints(self):
        """Test if generated subtelomere length is within specified bounds"""
        sequence, _ = self.region.generate()
        self.assertGreaterEqual(len(sequence), self.subtelomere_params.min_length)
        self.assertLessEqual(len(sequence), self.subtelomere_params.max_length)
    
    def test_repeat_density(self):
        """Test if repeat densities match specified parameters"""
        sequence, features = self.region.generate()
        
        # Calculate actual densities using sequence length
        total_length = len(sequence)
        satellite_length = sum(
            min(feature['end'], total_length) - feature['start']
            for feature in features
            if feature['type'] == 'satellite'
            and feature['start'] < total_length
        )
        transposon_length = sum(
            min(feature['end'], total_length) - feature['start']
            for feature in features
            if feature['type'] == 'transposon'
            and feature['start'] < total_length
        )
        
        # Calculate densities for sequence covered by features
        satellite_density = satellite_length / total_length
        transposon_density = transposon_length / total_length
        
        # Check densities with some tolerance
        self.assertAlmostEqual(
            satellite_density,
            self.subtelomere_params.satellite_density,
            delta=0.1
        )
        self.assertAlmostEqual(
            transposon_density,
            self.subtelomere_params.transposon_density,
            delta=0.1
        )
    
    def test_feature_distribution(self):
        """Test the spatial distribution of features"""
        _, features = self.region.generate()
        
        # Sort features by start position
        sorted_features = sorted(features, key=lambda x: x['start'])
        
        # Check for overlaps
        for i in range(len(sorted_features) - 1):
            current = sorted_features[i]
            next_feature = sorted_features[i + 1]
            self.assertLessEqual(current['end'], next_feature['start'])
    
    def test_satellite_properties(self):
        """Test properties of satellite repeat blocks"""
        _, features = self.region.generate()
        
        satellite_features = [f for f in features if f['type'] == 'satellite']
        for feature in satellite_features:
            # Check satellite type
            self.assertIn(
                feature.get('satellite_type'),
                ['satellite_1', 'satellite_2', 'satellite_3']
            )
            
            # Check monomer size and copies
            self.assertGreater(feature.get('monomer_size', 0), 0)
            self.assertGreater(feature.get('copies', 0), 0)

def visualize_region_structure(sequence: str, features: List[Dict],
                             title: str = "Region Structure"):
    """Helper function to visualize the structure of a region"""
    plt.figure(figsize=(15, 6))
    
    # Create color map for different feature types
    color_map = {
        'consensus': 'blue',
        'variant': 'red',
        'satellite': 'green',
        'transposon': 'orange'
    }
    
    # Plot features
    for feature in features:
        start = feature['start']
        end = feature['end']
        feature_type = feature['type']
        
        plt.axvspan(
            start, end,
            alpha=0.3,
            color=color_map.get(feature_type, 'gray'),
            label=feature_type
        )
    
    # Remove duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Feature Type")
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    # Run visualization for a sample telomere and subtelomere
    telomere_region_params = RegionBuilder(ChromosomeRegionType.TELOMERE)\
        .set_boundaries(0, 10000)\
        .set_gc_content(0.48)\
        .build()
    
    telomere_params = TelomereParams(
        min_length=3000,
        max_length=15000,
        mutation_rate_5_end=0.01,
        mutation_rate_3_end=0.05,
        variant_repeat_prob=0.1
    )
    
    telomere = TelomereRegion(telomere_region_params, telomere_params)
    sequence, features = telomere.generate()
    visualize_region_structure(sequence, features, "Telomere Structure")
    
    subtelomere_region_params = RegionBuilder(ChromosomeRegionType.SUBTELOMERE)\
        .set_boundaries(10000, 50000)\
        .set_gc_content(0.52)\
        .build()
        
    subtelomere_params = SubtelomereParams(
        min_length=10000,
        max_length=300000,
        satellite_density=0.4,
        transposon_density=0.2
    )
    
    subtelomere = SubtelomereRegion(subtelomere_region_params, subtelomere_params)
    sequence, features = subtelomere.generate()
    visualize_region_structure(sequence, features, "Subtelomere Structure")
    
    # Run tests
    unittest.main(argv=['first-arg-is-ignored'], exit=False)