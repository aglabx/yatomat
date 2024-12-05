import argparse
import sys
from pathlib import Path
from genome_generator import GenomeGenerator
import os

def create_parser():
    parser = argparse.ArgumentParser(
        description="YAToMat - Yet Another Tool for Making arTificial genomes"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Generate genome command
    gen_parser = subparsers.add_parser('generate', help='Generate artificial genome')
    gen_parser.add_argument('--scenario', type=str, default='simple',
                           help='Generation scenario to use')
    gen_parser.add_argument('--size', type=int, default=100,
                           help='Total genome size in Mb')
    gen_parser.add_argument('--chromosomes', type=int, default=6,
                           help='Number of chromosomes')
    gen_parser.add_argument('--output-dir', type=str, default='.',
                           help='Output directory')
    
    # List scenarios command
    list_parser = subparsers.add_parser('list-scenarios', 
                                       help='List available genome generation scenarios')
    
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    if args.command == 'list-scenarios':
        GenomeGenerator.list_scenarios()
        return 0
    
    if args.command == 'generate':
        print(f"=== Starting YAToMat Genome Generation ===")
        print(f"Using scenario: {args.scenario}")
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Change to output directory
        original_dir = Path.cwd()
        os.chdir(output_dir)
        
        try:
            generator = GenomeGenerator(
                total_size_mb=args.size,
                num_chromosomes=args.chromosomes
            )
            generator.execute_scenario(args.scenario)
            
        except ValueError as e:
            print(f"Error: {e}")
            GenomeGenerator.list_scenarios()
            return 1
            
        finally:
            # Return to original directory
            os.chdir(original_dir)
        
        print("\n=== Generation Complete ===")
        print(f"Output files are in: {output_dir}")
        return 0
    
    return 1

if __name__ == "__main__":
    sys.exit(main())