import random
from pathlib import Path
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from .gff_handling import GFFWriter, GenomeFeature
import math
import numpy as np
from multiprocessing import Pool
from itertools import islice
import mmap


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
    def generate(self, generator: "GenomeGenerator"):
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

    def generate(self, generator: "GenomeGenerator"):
        sizes = generator.calculate_chromosome_sizes()
        features = []

        # Write FASTA
        with open("genome.fasta", "w") as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"Processing {chr_name}...")

                # Add chromosome feature
                features.append(
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="chromosome",
                        start=1,
                        end=size,
                        attributes={"ID": chr_name},
                    )
                )

                f.write(f">{chr_name}\n")
                sequence = generator.generate_sequence(size)
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i : i + 60] + "\n")

        # Write GFF
        GFFWriter().write(features, "genome.gff")


class TelomericChromosomeScenario(Scenario):
    """Scenario that generates chromosomes with telomeric repeats at ends"""

    RIGHT_TELOMERE = "TTAGGG"  # Right end telomere
    LEFT_TELOMERE = "CCCTAA"  # Left end telomere (reverse complement)
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

    def generate(self, generator: "GenomeGenerator"):
        sizes = generator.calculate_chromosome_sizes()
        features = []

        with open("genome.fasta", "w") as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"Processing {chr_name}...")
                f.write(f">{chr_name}\n")

                left_telomere = self.generate_telomere(is_right_end=False)
                right_telomere = self.generate_telomere(is_right_end=True)
                main_size = size - len(left_telomere) - len(right_telomere)
                main_sequence = generator.generate_sequence(main_size)

                # Add chromosome feature
                features.append(
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="chromosome",
                        start=1,
                        end=size,
                        attributes={"ID": chr_name},
                    )
                )

                # Add telomere features
                features.append(
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="telomere",
                        start=1,
                        end=len(left_telomere),
                        attributes={
                            "ID": f"{chr_name}_tel_left",
                            "monomer": self.LEFT_TELOMERE,
                        },
                    )
                )

                features.append(
                    GenomeFeature(
                        seqid=chr_name,
                        source="genome_generator",
                        type="telomere",
                        start=size - len(right_telomere) + 1,
                        end=size,
                        attributes={
                            "ID": f"{chr_name}_tel_right",
                            "monomer": self.RIGHT_TELOMERE,
                        },
                    )
                )

                full_sequence = left_telomere + main_sequence + right_telomere
                for i in range(0, len(full_sequence), 60):
                    f.write(full_sequence[i : i + 60] + "\n")

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
        bases = ["A", "T", "G", "C"]
        self.alpha_monomer = "".join(
            random.choice(bases) for _ in range(self.ALPHA_MONOMER_LENGTH)
        )

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

    def generate(self, generator: "GenomeGenerator"):
        sizes = generator.calculate_chromosome_sizes()
        features = []

        with open("genome.fasta", "w") as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"\nProcessing {chr_name}...")
                f.write(f">{chr_name}\n")

                # Generate components
                left_telomere = self.generate_telomere(is_right_end=False)
                right_telomere = self.generate_telomere(is_right_end=True)
                centromere = self.generate_centromere()

                # Calculate positions
                left_tel_end = len(left_telomere)
                remaining_size = (
                    size - left_tel_end - len(right_telomere) - len(centromere)
                )
                left_arm_size = remaining_size // 2
                right_arm_size = remaining_size - left_arm_size

                cent_start = left_tel_end + left_arm_size + 1
                cent_end = cent_start + len(centromere) - 1
                right_tel_start = cent_end + right_arm_size + 1

                # Generate arms
                left_arm = generator.generate_sequence(left_arm_size)
                right_arm = generator.generate_sequence(right_arm_size)

                # Add features
                features.extend(
                    [
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="chromosome",
                            start=1,
                            end=size,
                            attributes={"ID": chr_name},
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="telomere",
                            start=1,
                            end=left_tel_end,
                            attributes={
                                "ID": f"{chr_name}_tel_left",
                                "monomer": self.LEFT_TELOMERE,
                            },
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="centromere",
                            start=cent_start,
                            end=cent_end,
                            attributes={
                                "ID": f"{chr_name}_cen",
                                "monomer": self.alpha_monomer,
                            },
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="telomere",
                            start=right_tel_start,
                            end=size,
                            attributes={
                                "ID": f"{chr_name}_tel_right",
                                "monomer": self.RIGHT_TELOMERE,
                            },
                        ),
                    ]
                )

                # Write sequence
                full_sequence = (
                    left_telomere + left_arm + centromere + right_arm + right_telomere
                )
                for i in range(0, len(full_sequence), 60):
                    f.write(full_sequence[i : i + 60] + "\n")

        # Write GFF
        GFFWriter().write(features, "genome.gff")


class CentromericFactoringScenario(CentromericChromosomeScenario):
    """Scenario that generates chromosomes with mutated centromeres"""

    @classmethod
    def name(cls) -> str:
        return "centromeric_factoring"

    def _mutate_sequence(self, sequence: str) -> str:
        sequence = list(sequence)
        length = len(sequence)

        # Substitutions (3%)
        sub_positions = random.sample(range(length), int(length * 0.03))
        bases = ["A", "T", "G", "C"]
        for pos in sub_positions:
            current = sequence[pos]
            available = [b for b in bases if b != current]
            sequence[pos] = random.choice(available)

        # Insertions (2%)
        ins_positions = random.sample(range(length), int(length * 0.02))
        for pos in sorted(ins_positions, reverse=True):
            sequence.insert(pos, random.choice(bases))

        # Deletions (2%)
        del_positions = random.sample(range(len(sequence)), int(length * 0.02))
        for pos in sorted(del_positions, reverse=True):
            sequence.pop(pos)

        # Random fragmentation (0.01%)
        frag_positions = random.sample(range(len(sequence)), int(length * 0.0001))
        for pos in sorted(frag_positions, reverse=True):
            sequence[pos : pos + 1] = ["N"] * random.randint(1, 10)

        return "".join(sequence)

    def _generate_centromere(self) -> str:
        # Generate base centromere using parent class
        centromere = super().generate_centromere()
        # Apply mutations
        return self._mutate_sequence(centromere)

    def generate(self, generator: "GenomeGenerator"):
        """Execute the scenario with mutated centromeres"""
        print("\nGenerating genome with mutated centromeres...")

        # Use parent class to generate basic structure
        sizes = generator.calculate_chromosome_sizes()
        features = []

        with open("genome.fasta", "w") as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"\nProcessing {chr_name}...")
                f.write(f">{chr_name}\n")

                # Generate components
                left_telomere = self.generate_telomere(is_right_end=False)
                right_telomere = self.generate_telomere(is_right_end=True)
                # Use modified centromere generation
                centromere = self._generate_centromere()

                # Calculate positions
                left_tel_end = len(left_telomere)
                remaining_size = (
                    size - left_tel_end - len(right_telomere) - len(centromere)
                )
                left_arm_size = remaining_size // 2
                right_arm_size = remaining_size - left_arm_size

                cent_start = left_tel_end + left_arm_size + 1
                cent_end = cent_start + len(centromere) - 1
                right_tel_start = cent_end + right_arm_size + 1

                # Generate arms
                left_arm = generator.generate_sequence(left_arm_size)
                right_arm = generator.generate_sequence(right_arm_size)

                # Add features
                features.extend(
                    [
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="chromosome",
                            start=1,
                            end=size,
                            attributes={"ID": chr_name},
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="telomere",
                            start=1,
                            end=left_tel_end,
                            attributes={
                                "ID": f"{chr_name}_tel_left",
                                "monomer": self.LEFT_TELOMERE,
                            },
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="centromere",
                            start=cent_start,
                            end=cent_end,
                            attributes={
                                "ID": f"{chr_name}_cen",
                                "monomer": self.alpha_monomer,
                            },
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="telomere",
                            start=right_tel_start,
                            end=size,
                            attributes={
                                "ID": f"{chr_name}_tel_right",
                                "monomer": self.RIGHT_TELOMERE,
                            },
                        ),
                    ]
                )

                # Write sequence
                full_sequence = (
                    left_telomere + left_arm + centromere + right_arm + right_telomere
                )
                for i in range(0, len(full_sequence), 60):
                    f.write(full_sequence[i : i + 60] + "\n")

        # Write GFF
        GFFWriter().write(features, "genome.gff")


class CentromericArrayInsertionScenario(CentromericFactoringScenario):
    ARRAY_MOTIF = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGG"  # Example array motif

    def __init__(self):
        super().__init__()
        self.array_size_range = (100, 1000)  # Size of each array insertion
        self.num_insertions_range = (10, 30)  # Number of insertions per region

    @classmethod
    def name(cls):
        return "centromeric_array_insertion"

    def generate_random_array(self):
        """Generate random array of ARRAY_MOTIF"""
        size = random.randint(*self.array_size_range)
        repeat_count = size // len(self.ARRAY_MOTIF)
        return self.ARRAY_MOTIF * repeat_count

    def insert_arrays_in_region(self, sequence, start, end):
        """Insert random arrays and return their positions"""
        modified_seq = list(sequence)
        num_insertions = random.randint(*self.num_insertions_range)
        insertions = []  # Track (position, length) of arrays

        for _ in range(num_insertions):
            insert_pos = random.randint(start, end)
            array = self.generate_random_array()
            modified_seq[insert_pos:insert_pos] = array
            insertions.append((insert_pos, len(array)))

        return "".join(modified_seq), insertions

    def generate(self, generator: "GenomeGenerator"):
        """Execute scenario with array insertions"""
        print("\nGenerating genome with centromeric mutations and array insertions...")

        sizes = generator.calculate_chromosome_sizes()
        features = []

        with open("genome.fasta", "w") as f:
            for chr_num, size in enumerate(sizes, 1):
                chr_name = f"chr{chr_num}"
                print(f"\nProcessing {chr_name}...")
                f.write(f">{chr_name}\n")

                # Generate basic components
                left_telomere = self.generate_telomere(is_right_end=False)
                right_telomere = self.generate_telomere(is_right_end=True)
                centromere = self._generate_centromere()

                # Insert arrays in centromere
                centromere, cen_arrays = self.insert_arrays_in_region(
                    centromere, start=len(centromere) // 4, end=3 * len(centromere) // 4
                )

                # Calculate initial positions
                left_tel_end = len(left_telomere)
                remaining_size = (
                    size - left_tel_end - len(right_telomere) - len(centromere)
                )
                left_arm_size = remaining_size // 2
                right_arm_size = remaining_size - left_arm_size

                # Generate and modify arms
                left_arm = generator.generate_sequence(left_arm_size)
                right_arm = generator.generate_sequence(right_arm_size)

                # Insert arrays in peritelomeric regions
                left_arm, left_arrays = self.insert_arrays_in_region(
                    left_arm, start=0, end=min(5000, len(left_arm))
                )
                right_arm, right_arrays = self.insert_arrays_in_region(
                    right_arm, start=max(0, len(right_arm) - 5000), end=len(right_arm)
                )

                # Calculate final positions
                cent_start = left_tel_end + len(left_arm) + 1
                cent_end = cent_start + len(centromere) - 1
                right_tel_start = cent_end + len(right_arm) + 1

                # Add features
                features.extend(
                    [
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="chromosome",
                            start=1,
                            end=size,
                            attributes={"ID": chr_name},
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="telomere",
                            start=1,
                            end=left_tel_end,
                            attributes={
                                "ID": f"{chr_name}_tel_left",
                                "monomer": self.LEFT_TELOMERE,
                            },
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="centromere",
                            start=cent_start,
                            end=cent_end,
                            attributes={
                                "ID": f"{chr_name}_cen",
                                "monomer": self.alpha_monomer,
                                "has_arrays": "true",
                            },
                        ),
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="telomere",
                            start=right_tel_start,
                            end=size,
                            attributes={
                                "ID": f"{chr_name}_tel_right",
                                "monomer": self.RIGHT_TELOMERE,
                            },
                        ),
                    ]
                )

                # Add array features
                for pos, length in left_arrays:
                    abs_pos = left_tel_end + pos
                    features.append(
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="array",
                            start=abs_pos,
                            end=abs_pos + length,
                            attributes={
                                "ID": f"{chr_name}_array_left_{pos}",
                                "motif": self.ARRAY_MOTIF,
                                "region": "peritelomeric_left",
                            },
                        )
                    )

                for pos, length in cen_arrays:
                    abs_pos = cent_start + pos
                    features.append(
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="array",
                            start=abs_pos,
                            end=abs_pos + length,
                            attributes={
                                "ID": f"{chr_name}_array_cen_{pos}",
                                "motif": self.ARRAY_MOTIF,
                                "region": "centromeric",
                            },
                        )
                    )

                for pos, length in right_arrays:
                    abs_pos = cent_end + pos + 1
                    features.append(
                        GenomeFeature(
                            seqid=chr_name,
                            source="genome_generator",
                            type="array",
                            start=abs_pos,
                            end=abs_pos + length,
                            attributes={
                                "ID": f"{chr_name}_array_right_{pos}",
                                "motif": self.ARRAY_MOTIF,
                                "region": "peritelomeric_right",
                            },
                        )
                    )

                # Write sequence
                full_sequence = (
                    left_telomere + left_arm + centromere + right_arm + right_telomere
                )
                for i in range(0, len(full_sequence), 60):
                    f.write(full_sequence[i : i + 60] + "\n")

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
        CentromericChromosomeScenario.name(): CentromericChromosomeScenario,
        CentromericFactoringScenario.name(): CentromericFactoringScenario,
        CentromericArrayInsertionScenario.name(): CentromericArrayInsertionScenario,
    }

    def __init__(self, total_size_mb=100, num_chromosomes=6, fastq_params=None):
        """
        Initialize the genome generator.

        Args:
            total_size_mb (int): Total genome size in megabases
            num_chromosomes (int): Number of chromosomes to generate
        """
        self.total_size = total_size_mb * 1_000_000  # Convert to bases
        self.num_chromosomes = num_chromosomes
        self.bases = ["A", "T", "G", "C"]
        self.fastq_params = fastq_params or {}
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
        return "".join(random.choice(self.bases) for _ in range(length))

    def write_fasta(self, output_path="genome.fasta"):
        """
        Write generated chromosomes to a FASTA file.

        Args:
            output_path (str): Path to output FASTA file
        """
        print(f"\nWriting sequences to {output_path}")
        sizes = self.calculate_chromosome_sizes()

        with open(output_path, "w") as f:
            for chr_num, size in enumerate(sizes, 1):
                print(f"Processing chromosome {chr_num}...")
                f.write(f">chr{chr_num}\n")
                sequence = self.generate_sequence(size)

                # Write sequence in lines of 60 characters
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i : i + 60] + "\n")

        print(f"Finished writing genome to {output_path}")

    def generate_fastq(self):
        """Generate FASTQ files from genome"""
        if not self.fastq_params:
            return

        generator = FastqGenerator(**self.fastq_params)
        generator.generate_fastq()


class FastqGenerator:
    def __init__(self, read_length=151, coverage=30, is_paired=True, batch_size=10000):
        self.read_length = read_length
        self.coverage = coverage
        self.is_paired = is_paired
        self.error_rate = 0.01
        self.batch_size = batch_size
        self.complement = str.maketrans("ATGC", "TAGC")

    def _process_batch(self, args):
        sequence, start_idx, num_reads = args
        reads = []
        # Pre-generate random positions for efficiency
        positions = np.random.randint(0, len(sequence) - self.read_length, num_reads)
        error_masks = np.random.random((num_reads, self.read_length)) < self.error_rate
        quality = chr(40) * self.read_length  # Fixed quality for speed

        for i, pos in enumerate(positions):
            read = sequence[pos : pos + self.read_length]
            if error_masks[i].any():
                read = list(read)
                for j in np.where(error_masks[i])[0]:
                    read[j] = np.random.choice(list("ATGC".replace(read[j], "")))
                read = "".join(read)
            reads.append((read, quality))

            if self.is_paired:
                mate_pos = min(pos + 500, len(sequence) - self.read_length)
                mate = sequence[mate_pos : mate_pos + self.read_length]
                mate = mate.translate(self.complement)[::-1]
                reads.append((mate, quality))

        return reads

    def generate_fastq(self, fasta_path="genome.fasta", output_prefix="genome"):
        # Memory-map the FASTA file
        with open(fasta_path, "r") as f:
            sequence = "".join(line.strip() for line in f if not line.startswith(">"))

        num_reads = (len(sequence) * self.coverage) // self.read_length
        batches = [
            (sequence, i, min(self.batch_size, num_reads - i))
            for i in range(0, num_reads, self.batch_size)
        ]

        # Process batches in parallel
        with Pool() as pool:
            if self.is_paired:
                with open(
                    f"{output_prefix}_R1.fastq", "w", buffering=8192 * 1024
                ) as f1, open(
                    f"{output_prefix}_R2.fastq", "w", buffering=8192 * 1024
                ) as f2:
                    for batch_num, reads in enumerate(
                        pool.imap(self._process_batch, batches)
                    ):
                        for i in range(0, len(reads), 2):
                            read1, qual1 = reads[i]
                            read2, qual2 = reads[i + 1]
                            read_num = batch_num * self.batch_size + i // 2
                            f1.write(f"@read_{read_num}/1\n{read1}\n+\n{qual1}\n")
                            f2.write(f"@read_{read_num}/2\n{read2}\n+\n{qual2}\n")
            else:
                with open(f"{output_prefix}.fastq", "w", buffering=8192 * 1024) as f:
                    for batch_num, reads in enumerate(
                        pool.imap(self._process_batch, batches)
                    ):
                        for i, (read, qual) in enumerate(reads):
                            read_num = batch_num * self.batch_size + i
                            f.write(f"@read_{read_num}\n{read}\n+\n{qual}\n")
