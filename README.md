# YATOMAT - Yet Another Tool for Making Artificial Genomes

![YATOMAT Logo](https://private-user-images.githubusercontent.com/142793/381877894-fcb3d1b4-53cb-4fac-899a-e0c411b89c9b.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3MzA0MTc0OTYsIm5iZiI6MTczMDQxNzE5NiwicGF0aCI6Ii8xNDI3OTMvMzgxODc3ODk0LWZjYjNkMWI0LTUzY2ItNGZhYy04OTlhLWUwYzQxMWI4OWM5Yi5wbmc_WC1BbXotQWxnb3JpdGhtPUFXUzQtSE1BQy1TSEEyNTYmWC1BbXotQ3JlZGVudGlhbD1BS0lBVkNPRFlMU0E1M1BRSzRaQSUyRjIwMjQxMDMxJTJGdXMtZWFzdC0xJTJGczMlMkZhd3M0X3JlcXVlc3QmWC1BbXotRGF0ZT0yMDI0MTAzMVQyMzI2MzZaJlgtQW16LUV4cGlyZXM9MzAwJlgtQW16LVNpZ25hdHVyZT1mODZhZWEzZmNiYWQzNzZiMmZhYzk3ODUwOWFlYTcyMzkxMTY3ZDUxODEyMGY3NTNkOWE0NzVjYmIxNGUzMDU4JlgtQW16LVNpZ25lZEhlYWRlcnM9aG9zdCJ9.nJ0n1Jc5j0FuRg7QXO22GZHI8Hf4SsP3Mcq95wHST78)

## Overview

YATOMAT is a comprehensive tool designed to generate realistic artificial genomes. It allows users to create chromosome sequences with detailed annotations, including telomeres, centromeres, pericentromeres, and other chromosomal regions. The tool supports various configurations and outputs in FASTA and GFF3 formats.

## Features

- **Flexible Configuration**: Customize chromosome parameters through a JSON configuration file.
- **Realistic Sequence Generation**: Generate sequences with realistic GC content, repeat densities, and chromatin states.
- **Detailed Annotations**: Output detailed annotations in GFF3 format, including telomeres, centromeres, and other features.
- **Compression Support**: Optionally compress output files in gzip format.
- **Reproducibility**: Set a random seed for reproducible results.
- **Logging**: Detailed logging for tracking the generation process.

## Installation

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/yourusername/yatomat.git
cd yatomat
pip install -r requirements.txt
```

## Usage

### Command Line Interface

Run the main script with the required configuration file and optional parameters:

```bash
python yatomat.py config.json --output-dir output --prefix chr1 --compress --seed 42
```

### Configuration File

The configuration file is a JSON file that specifies the parameters for chromosome generation. Below is an example configuration file (`example_config.json`):

```json
{
    "total_length": 80000000,
    "p_arm": {
        "length": 30000000,
        "gc_content": 0.45,
        "gc_std": 0.02,
        "repeat_density": 0.3,
        "tad_organization": 0.8,
        "chromatin_state": "EUCHROMATIN"
    },
    "q_arm": {
        "length": 50000000,
        "gc_content": 0.48,
        "gc_std": 0.02,
        "repeat_density": 0.35,
        "tad_organization": 0.75,
        "chromatin_state": "EUCHROMATIN"
    },
    "telomere": {
        "max_length": 10000
    },
    "subtelomere": {
        "max_length": 20000
    },
    "centromere": {
        "max_core_length": 4000000
    },
    "pericentromere": {
        "min_total_length": 1000000
    },
    "boundary": {
        "transition_length": 50000
    }
}
```

### Output

The tool generates two main output files:
- **FASTA File**: Contains the generated chromosome sequence.
- **GFF3 File**: Contains annotations for the generated features.

## Modules

### Core Modules

- **core.py**: Basic functions for sequence generation and mutation handling.
- **repeats.py**: Functions for handling various types of repeats and HOR structures.
- **regions/common.py**: Base classes and interfaces for chromosome regions.
- **regions/telomeres.py**: Functions for generating telomeric and subtelomeric regions.
- **regions/centromeres.py**: Functions for generating centromeric regions.
- **regions/pericentromeres.py**: Functions for generating pericentromeric regions.
- **regions/boundaries.py**: Functions for generating boundary regions between different chromatin states.

### Main Script

- **yatomat.py**: The main script that integrates all modules and generates the final output.

## Testing

To ensure the correctness of the generated sequences and annotations, comprehensive test scripts are provided for each module. Run the tests using:

```bash
python -m unittest discover tests
```

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests.

## License

This project is licensed under the MIT License.

## Contact

For any questions or issues, please open an issue on GitHub or contact the maintainer at ad3002@gmail.com.
