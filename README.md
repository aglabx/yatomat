# YATOMAT - Yet Another Tool for Making Artificial Genomes

<img src="https://raw.githubusercontent.com/aglabx/yatomat/refs/heads/main/yatomat.png" alt="YATOMAT Logo" width="350">

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

### From PyPI

You can install YATOMAT directly from PyPI:

```bash
pip install yatomat
```

### From Source

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/yourusername/yatomat.git
cd yatomat
pip install .
```

## Usage

After installation, YATOMAT can be run from anywhere using the command line:

```bash
yatomat config.json --output-dir output --prefix chr1 --compress --seed 42
```

Available options:
- `config.json`: Path to configuration file (required)
- `--output-dir`: Output directory (default: output)
- `--prefix`: Prefix for output files (default: chr)
- `--compress`: Compress output files (optional)
- `--seed`: Random seed for reproducibility (optional)

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
pytest .
```

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests.

## License

This project is licensed under the MIT License.

## Contact

For any questions or issues, please open an issue on GitHub or contact the maintainer at ad3002@gmail.com.
