import argparse
import sys
from pathlib import Path
from .genome_generator import GenomeGenerator
import os


def create_parser():
    parser = argparse.ArgumentParser(
        description="YAToMat - Yet Another Tool for Making arTificial genomes"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Generate genome command
    gen_parser = subparsers.add_parser("generate", help="Generate artificial genome")
    gen_parser.add_argument(
        "--scenario", type=str, default="simple", help="Generation scenario to use"
    )
    gen_parser.add_argument(
        "--size", type=int, default=100, help="Total genome size in Mb"
    )
    gen_parser.add_argument(
        "--chromosomes", type=int, default=6, help="Number of chromosomes"
    )
    gen_parser.add_argument(
        "--output-dir", type=str, default=".", help="Output directory"
    )

    # List scenarios command
    list_parser = subparsers.add_parser(
        "list-scenarios", help="List available genome generation scenarios"
    )

    gen_parser.add_argument("--fastq", action="store_true", help="Generate FASTQ files")
    gen_parser.add_argument(
        "--read-length", type=int, default=151, help="Read length for FASTQ generation"
    )
    gen_parser.add_argument(
        "--coverage", type=int, default=30, help="Coverage depth for FASTQ generation"
    )
    gen_parser.add_argument(
        "--single-end",
        action="store_true",
        help="Generate single-end reads instead of paired-end",
    )
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    if args.command == "list-scenarios":
        GenomeGenerator.list_scenarios()
        return 0

    if args.command == "generate":
        print(f"=== Starting YAToMat Genome Generation ===")
        print(f"Using scenario: {args.scenario}")

        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Change to output directory
        original_dir = Path.cwd()
        os.chdir(output_dir)

        try:
            # Add FASTQ parameters
            fastq_params = None
            if args.fastq:
                fastq_params = {
                    "read_length": args.read_length,
                    "coverage": args.coverage,
                    "is_paired": not args.single_end,
                }

            generator = GenomeGenerator(
                total_size_mb=args.size,
                num_chromosomes=args.chromosomes,
                fastq_params=fastq_params,
            )

            # Execute scenario
            generator.execute_scenario(args.scenario)

            # Generate FASTQ if requested
            if args.fastq:
                print("Generating FASTQ files...")
                generator.generate_fastq()

        except ValueError as e:
            print(f"Error: {e}")
            GenomeGenerator.list_scenarios()
            return 1
        finally:
            os.chdir(original_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
