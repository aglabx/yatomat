import random
# Remove argparse import
from pathlib import Path
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from gff_handling import GFFWriter, GenomeFeature

class Scenario(ABC):
    """Base class for genome generation scenarios"""
    
    @classmethod
    @abstractmethod
    def name(cls) -> str:
        """Return scenario name"""
        pass
    
    @classmethod
    @abstractmethod
    def description(cls) -> str:
        """Return scenario description"""
        pass
    
    @abstractmethod
    def generate(self, generator: 'GenomeGenerator'):
        """Execute the scenario"""
        pass

class SimpleChromosomeScenario(Scenario):
    """Basic scenario with decreasing chromosome sizes"""
    
    @classmethod
    def name(cls) -> str:
        return "simple"
    
    @classmethod
    def description(cls) -> str:
        return "Generate simple chromosomes with decreasing sizes"
    
    def generate(self, generator: 'GenomeGenerator'):
        sizes = generator.calculate_chromosome_sizes()
        features = []
        
        # Write FASTA
        with open("genome.fasta", 'w') as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"Processing {chr_name}...")
                
                # Add chromosome feature
                features.append(GenomeFeature(
                    seqid=chr_name,
                    source="genome_generator",
                    type="chromosome",
                    start=1,
                    end=size,
                    attributes={"ID": chr_name}
                ))
                
                f.write(f'>{chr_name}\n')
                sequence = generator.generate_sequence(size)
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')
        
        # Write GFF
        GFFWriter().write(features, "genome.gff")

class TelomericChromosomeScenario(Scenario):
    """Scenario that generates chromosomes with telomeric repeats at ends"""
    
    RIGHT_TELOMERE = "TTAGGG"  # Right end telomere
    LEFT_TELOMERE = "CCCTAA"   # Left end telomere (reverse complement)
    BASE_TELOMERE_LENGTH = 3000  # 3kb base length
    
    @classmethod
    def name(cls) -> str:
        return "telomeric"
    
    @classmethod
    def description(cls) -> str:
        return "Generate chromosomes with oriented telomeric repeats (CCCTAA)n-(sequence)-(TTAGGG)n"
    
    def generate_telomere(self, is_right_end: bool = True) -> str:
        """
        Generate telomeric sequence with random variation
        
        Args:
            is_right_end: If True, use TTAGGG, else use CCCTAA
        """
        variation = random.uniform(-0.2, 0.2)  # ±20% variation
        length = int(self.BASE_TELOMERE_LENGTH * (1 + variation))
        repeat = self.RIGHT_TELOMERE if is_right_end else self.LEFT_TELOMERE
        repeats = length // len(repeat)
        return repeat * repeats
    
    def generate(self, generator: 'GenomeGenerator'):
        sizes = generator.calculate_chromosome_sizes()
        features = []
        
        with open("genome.fasta", 'w') as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"Processing {chr_name}...")
                f.write(f'>{chr_name}\n')
                
                left_telomere = self.generate_telomere(is_right_end=False)
                right_telomere = self.generate_telomere(is_right_end=True)
                main_size = size - len(left_telomere) - len(right_telomere)
                main_sequence = generator.generate_sequence(main_size)
                
                # Add chromosome feature
                features.append(GenomeFeature(
                    seqid=chr_name,
                    source="genome_generator",
                    type="chromosome",
                    start=1,
                    end=size,
                    attributes={"ID": chr_name}
                ))
                
                # Add telomere features
                features.append(GenomeFeature(
                    seqid=chr_name,
                    source="genome_generator",
                    type="telomere",
                    start=1,
                    end=len(left_telomere),
                    attributes={
                        "ID": f"{chr_name}_tel_left",
                        "monomer": self.LEFT_TELOMERE
                    }
                ))
                
                features.append(GenomeFeature(
                    seqid=chr_name,
                    source="genome_generator",
                    type="telomere",
                    start=size - len(right_telomere) + 1,
                    end=size,
                    attributes={
                        "ID": f"{chr_name}_tel_right",
                        "monomer": self.RIGHT_TELOMERE
                    }
                ))
                
                full_sequence = left_telomere + main_sequence + right_telomere
                for i in range(0, len(full_sequence), 60):
                    f.write(full_sequence[i:i+60] + '\n')
        
        # Write GFF
        GFFWriter().write(features, "genome.gff")

class CentromericChromosomeScenario(Scenario):
    """Scenario that generates chromosomes with telomeres and centromeres"""
    
    RIGHT_TELOMERE = "TTAGGG"
    LEFT_TELOMERE = "CCCTAA"
    BASE_TELOMERE_LENGTH = 3000
    BASE_CENTROMERE_LENGTH = 1_000_000  # 1Mb
    ALPHA_MONOMER_LENGTH = 171
    
    def __init__(self):
        # Generate one alpha satellite monomer to use for all centromeres
        bases = ['A', 'T', 'G', 'C']
        self.alpha_monomer = ''.join(random.choice(bases) for _ in range(self.ALPHA_MONOMER_LENGTH))
    
    @classmethod
    def name(cls) -> str:
        return "centromeric"
    
    @classmethod
    def description(cls) -> str:
        return "Generate chromosomes with telomeres and alpha-satellite centromeres"
    
    def generate_telomere(self, is_right_end: bool = True) -> str:
        """Generate telomeric sequence with random variation"""
        variation = random.uniform(-0.2, 0.2)
        length = int(self.BASE_TELOMERE_LENGTH * (1 + variation))
        repeat = self.RIGHT_TELOMERE if is_right_end else self.LEFT_TELOMERE
        repeats = length // len(repeat)
        return repeat * repeats
    
    def generate_centromere(self) -> str:
        """Generate centromeric sequence with random size variation"""
        variation = random.uniform(-0.2, 0.2)
        length = int(self.BASE_CENTROMERE_LENGTH * (1 + variation))
        repeats = length // len(self.alpha_monomer)
        return self.alpha_monomer * repeats
    
    def generate(self, generator: 'GenomeGenerator'):
        sizes = generator.calculate_chromosome_sizes()
        features = []
        
        with open("genome.fasta", 'w') as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"\nProcessing {chr_name}...")
                f.write(f'>{chr_name}\n')
                
                # Generate components
                left_telomere = self.generate_telomere(is_right_end=False)
                right_telomere = self.generate_telomere(is_right_end=True)
                centromere = self.generate_centromere()
                
                # Calculate positions
                left_tel_end = len(left_telomere)
                remaining_size = size - left_tel_end - len(right_telomere) - len(centromere)
                left_arm_size = remaining_size // 2
                right_arm_size = remaining_size - left_arm_size
                
                cent_start = left_tel_end + left_arm_size + 1
                cent_end = cent_start + len(centromere) - 1
                right_tel_start = cent_end + right_arm_size + 1
                
                # Generate arms
                left_arm = generator.generate_sequence(left_arm_size)
                right_arm = generator.generate_sequence(right_arm_size)
                
                # Add features
                features.extend([
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="chromosome",
                        start=1,
                        end=size,
                        attributes={"ID": chr_name}
                    ),
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="telomere",
                        start=1,
                        end=left_tel_end,
                        attributes={
                            "ID": f"{chr_name}_tel_left",
                            "monomer": self.LEFT_TELOMERE
                        }
                    ),
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="centromere",
                        start=cent_start,
                        end=cent_end,
                        attributes={
                            "ID": f"{chr_name}_cen",
                            "monomer": self.alpha_monomer
                        }
                    ),
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="telomere",
                        start=right_tel_start,
                        end=size,
                        attributes={
                            "ID": f"{chr_name}_tel_right",
                            "monomer": self.RIGHT_TELOMERE
                        }
                    )
                ])
                
                # Write sequence
                full_sequence = (left_telomere + left_arm + 
                               centromere + right_arm + right_telomere)
                for i in range(0, len(full_sequence), 60):
                    f.write(full_sequence[i:i+60] + '\n')
        
        # Write GFF
        GFFWriter().write(features, "genome.gff")

class GenomeGenerator:
    """
    A class for generating artificial chromosomes with specified size distributions.
    
    Creates a genome with multiple chromosomes where each subsequent chromosome
    is approximately 10% smaller than the previous one. The chromosomes are
    composed of random DNA sequences (A, T, G, C).
    """
    
    # Register scenarios here
    SCENARIOS = {
        SimpleChromosomeScenario.name(): SimpleChromosomeScenario,
        TelomericChromosomeScenario.name(): TelomericChromosomeScenario,
        CentromericChromosomeScenario.name(): CentromericChromosomeScenario
    }

    def __init__(self, total_size_mb=100, num_chromosomes=6):
        """
        Initialize the genome generator.
        
        Args:
            total_size_mb (int): Total genome size in megabases
            num_chromosomes (int): Number of chromosomes to generate
        """
        self.total_size = total_size_mb * 1_000_000  # Convert to bases
        self.num_chromosomes = num_chromosomes
        self.bases = ['A', 'T', 'G', 'C']
        print(f"Initializing genome generator:")
        print(f"- Total genome size: {total_size_mb} Mb")
        print(f"- Number of chromosomes: {num_chromosomes}")
    
    @classmethod
    def list_scenarios(cls):
        """Print all available scenarios"""
        print("\nAvailable scenarios:")
        for name, scenario_class in cls.SCENARIOS.items():
            print(f"- {name}: {scenario_class.description()}")
    
    def execute_scenario(self, scenario_name: str):
        """Execute a specific scenario"""
        if scenario_name not in self.SCENARIOS:
            raise ValueError(f"Unknown scenario: {scenario_name}")
        
        scenario = self.SCENARIOS[scenario_name]()
        scenario.generate(self)
    
    def calculate_chromosome_sizes(self):
        """
        Calculate sizes for each chromosome using a decreasing distribution.
        
        Returns:
            list: List of chromosome sizes in base pairs
        """
        print("\nCalculating chromosome sizes...")
        sizes = []
        remaining_size = self.total_size
        
        for i in range(self.num_chromosomes):
            if i == self.num_chromosomes - 1:
                sizes.append(remaining_size)
                print(f"- Chr{i+1}: {remaining_size:,} bp (remaining bases)")
                break
                
            size = int(remaining_size * 0.26)  # First chr is ~26% of remaining
            sizes.append(size)
            remaining_size -= size
            print(f"- Chr{i+1}: {size:,} bp")
            
        return sizes
    
    def generate_sequence(self, length):
        """
        Generate a random DNA sequence of specified length.
        
        Args:
            length (int): Length of sequence to generate
            
        Returns:
            str: Random DNA sequence
        """
        print(f"Generating sequence of length {length:,} bp...")
        return ''.join(random.choice(self.bases) for _ in range(length))
    
    def write_fasta(self, output_path="genome.fasta"):
        """
        Write generated chromosomes to a FASTA file.
        
        Args:
            output_path (str): Path to output FASTA file
        """
        print(f"\nWriting sequences to {output_path}")
        sizes = self.calculate_chromosome_sizes()
        
        with open(output_path, 'w') as f:
            for chr_num, size in enumerate(sizes, 1):
                print(f"Processing chromosome {chr_num}...")
                f.write(f'>chr{chr_num}\n')
                sequence = self.generate_sequence(size)
                
                # Write sequence in lines of 60 characters
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')
        
        print(f"Finished writing genome to {output_path}")
