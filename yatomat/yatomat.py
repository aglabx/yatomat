#!/usr/bin/env python3
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
import argparse
import logging
from pathlib import Path
import json
import sys

from yatomat.regions.chromosome_regions import (
    ChromosomeParams, ChromosomeArmParams, ChromosomeAssembly,
    ChromatinState, RegionParams, ChromosomeRegionType,
    GradientParams, GradientGenerator
)
from yatomat.regions.boundaries import BoundaryParams
from yatomat.regions.centromeres import CentromereParams
from yatomat.regions.telomers import TelomereParams, SubtelomereParams
from yatomat.regions.pericentromeres import PericentromereParams

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class OutputConfig:
    """Configuration for output files"""
    output_dir: Path
    prefix: str
    compress: bool = False

    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    @property
    def fasta_path(self) -> Path:
        suffix = '.fa.gz' if self.compress else '.fa'
        return self.output_dir / f"{self.prefix}{suffix}"

    @property
    def gff_path(self) -> Path:
        suffix = '.gff3.gz' if self.compress else '.gff3'
        return self.output_dir / f"{self.prefix}{suffix}"


def write_fasta(sequence: str, name: str, path: Path, compress: bool = False):
    """Write sequence to FASTA format"""
    # Split sequence into lines of 60 characters
    sequence_lines = [sequence[i:i+60] for i in range(0, len(sequence), 60)]

    if compress:
        import gzip
        with gzip.open(path, 'wt') as f:
            f.write(f">{name}\n")
            f.write('\n'.join(sequence_lines))
            f.write('\n')
    else:
        with open(path, 'w') as f:
            f.write(f">{name}\n")
            f.write('\n'.join(sequence_lines))
            f.write('\n')

    logger.info(f"Wrote FASTA file to {path}")

def format_gff_attributes(attributes: Dict) -> str:
    """Format GFF attributes dictionary to string"""
    formatted = []
    for key, value in attributes.items():
        if isinstance(value, (list, tuple)):
            value = ','.join(map(str, value))
        elif isinstance(value, bool):
            value = str(value).lower()
        elif isinstance(value, (int, float)):
            value = str(value)
        formatted.append(f"{key}={value}")
    return ';'.join(formatted)

def write_gff(features: List[Dict], sequence_name: str, path: Path, compress: bool = False):
    """Write features to GFF3 format"""
    # Filter features to include only specified types
    filtered_features = features
    # allowed_feature_types = {'telomere', 'centromere', 'pericentromere', 'repeats'}
    # filtered_features = [f for f in features if f['type'] in allowed_feature_types]



    gff_lines = [
        "##gff-version 3",
        f"##sequence-region {sequence_name} 1 {max(f['end'] for f in filtered_features)}"
    ]

    for feature in filtered_features:
        # Basic GFF fields
        seqid = sequence_name
        source = "YATOMAT"
        type = feature['type']
        start = feature['start'] + 1  # Convert to 1-based coordinates
        end = feature['end']
        score = '.'
        strand = feature.get('strand', '.')
        phase = '.'

        # Prepare attributes
        attributes = {k: v for k, v in feature.items()
                     if k not in {'type', 'start', 'end', 'strand'}}
        attr_string = format_gff_attributes(attributes)

        # Combine fields
        line = f"{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attr_string}"
        gff_lines.append(line)

    if compress:
        import gzip
        with gzip.open(path, 'wt') as f:
            f.write('\n'.join(gff_lines))
            f.write('\n')
    else:
        with open(path, 'w') as f:
            f.write('\n'.join(gff_lines))
            f.write('\n')

    logger.info(f"Wrote GFF file to {path}")

def load_config(config_path: Path) -> Dict:
    """Load configuration from JSON file"""
    with open(config_path) as f:
        config = json.load(f)
    return config

def create_chromosome_params(config: Dict) -> ChromosomeParams:
    """Create ChromosomeParams from configuration"""
    p_arm_params = ChromosomeArmParams(
        length=config['p_arm']['length'],
        gc_content=config['p_arm'].get('gc_content', 0.45),
        gc_std=config['p_arm'].get('gc_std', 0.02),
        repeat_density=config['p_arm'].get('repeat_density', 0.3),
        tad_organization=config['p_arm'].get('tad_organization', 0.8),
        chromatin_state=ChromatinState[config['p_arm'].get('chromatin_state', 'EUCHROMATIN')]
    )

    q_arm_params = ChromosomeArmParams(
        length=config['q_arm']['length'],
        gc_content=config['q_arm'].get('gc_content', 0.48),
        gc_std=config['q_arm'].get('gc_std', 0.02),
        repeat_density=config['q_arm'].get('repeat_density', 0.35),
        tad_organization=config['q_arm'].get('tad_organization', 0.75),
        chromatin_state=ChromatinState[config['q_arm'].get('chromatin_state', 'EUCHROMATIN')]
    )

    telomere_params = TelomereParams(**config.get('telomere', {}))
    subtelomere_params = SubtelomereParams(**config.get('subtelomere', {}))
    centromere_params = CentromereParams(**config.get('centromere', {}))
    pericentromere_params = PericentromereParams(**config.get('pericentromere', {}))
    boundary_params = BoundaryParams(**config.get('boundary', {}))

    return ChromosomeParams(
        total_length=config['total_length'],
        p_arm_params=p_arm_params,
        q_arm_params=q_arm_params,
        telomere_params=telomere_params,
        subtelomere_params=subtelomere_params,
        centromere_params=centromere_params,
        pericentromere_params=pericentromere_params,
        boundary_params=boundary_params
    )


def print_welcome_banner():
    """Print YATOMAT welcome banner with Braille tomato art"""
    banner = """
----------------------------------------------------------------
🍅 YATOMAT - Yet Another Tool for Making Artificial genomes 🧬
================================================================
                ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⠻⣶⡆⠀⠿⠀⣶⠒⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                ⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣴⠾⠛⢹⣶⡤⢶⣿⡟⠶⠦⠄⠀⠀⠀⠀⠀⠀⠀⠀
                ⠀⠀⠀⠀⠀⣠⣶⣤⣤⣤⣤⣴⠂⠸⠋⢀⣄⡉⠓⠀⠲⣶⣾⣿⣷⣄⠀⠀⠀⠀
                ⠀⠀⠀⢀⣾⡿⠋⠁⣠⣤⣿⡟⢀⣠⣾⣿⣿⣿⣷⣶⣤⣼⣿⣿⣿⣿⣆⠀⠀⠀
                ⠀⠀⠀⣾⡟⠀⣰⣿⣿⣿⣿⣷⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡄⠀⠀
                ⠀⠀⢸⡿⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀
                ⠀⠀⢸⡇⢰⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀
                ⠀⠀⢸⣿⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀
                ⠀⠀⠸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠁⠀⠀
                ⠀⠀⠀⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⠀⠀
                ⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⠁⠀⠀⠀⠀
                ⠀⠀⠀⠀⠀⠀⠉⠛⠿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠟⠋⠀⠀⠀⠀⠀⠀⠀
                ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀    
================================================================
                Artificial Genome Generator                        
----------------------------------------------------------------
"""
    print(banner)

def main():
    print_welcome_banner()  # Add this line at the start of main()
    parser = argparse.ArgumentParser(description='YATOMAT - Yet Another Tool for Making Artificial genomes')
    parser.add_argument('config', type=Path, help='Path to configuration file')
    parser.add_argument('--output-dir', type=Path, default=Path('output'),
                        help='Output directory (default: output)')
    parser.add_argument('--prefix', type=str, default='chr',
                        help='Prefix for output files (default: chr)')
    parser.add_argument('--compress', action='store_true',
                        help='Compress output files')
    parser.add_argument('--seed', type=int, help='Random seed')

    args = parser.parse_args()

    if args.seed is not None:
        import numpy as np
        np.random.seed(args.seed)
        import random
        random.seed(args.seed)

    try:
        config = load_config(args.config)
        chromosome_params = create_chromosome_params(config)

        output_config = OutputConfig(
            output_dir=args.output_dir,
            prefix=args.prefix,
            compress=args.compress
        )

        assembly = ChromosomeAssembly(chromosome_params)
        sequence, features = assembly.generate()

        write_fasta(sequence, output_config.prefix, output_config.fasta_path,
                   compress=output_config.compress)
        write_gff(features, output_config.prefix, output_config.gff_path,
                 compress=output_config.compress)

        logger.info(f"Successfully generated chromosome of length {len(sequence)} "
                   f"with {len(features)} features")

    except Exception as e:
        logger.error(f"Error during chromosome generation: {e}", exc_info=True)
        sys.exit(1)

if __name__ == '__main__':
    main()
