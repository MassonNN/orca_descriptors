# ORCA Descriptors

A Python library for automatic calculation of quantum chemical descriptors for QSAR analysis using ORCA quantum chemistry software.

## Installation

### Using pip

```bash
pip install orca-descriptors
```

### Using Poetry (development)

```bash
poetry install
```

## Usage

### As a Python Library

```python
from orca_descriptors import Orca, Molecule

# Initialize ORCA calculator
orca = Orca(
    script_path="orca",
    functional="PBE0",
    basis_set="def2-SVP",
    method_type="Opt",
    dispersion_correction="D3BJ",
    solvation_model="COSMO(Water)",
    n_processors=8,
)

# Create molecule from SMILES
mol = Molecule.from_smiles("C1=CC=CC=C1")

# Calculate descriptors
homo = orca.homo_energy(mol)
lumo = orca.lumo_energy(mol)
gap = orca.gap_energy(mol)
```

### As a Command-Line Utility

After installation, you can use `orca_descriptors` as a command-line tool:

#### Run Benchmark

Calibrate time estimation by running a benchmark calculation. The benchmark uses benzene (C1=CC=CC=C1) as a standard test molecule for machine calibration:

```bash
orca_descriptors run_benchmark
```

#### Estimate Calculation Time

Estimate calculation time for a molecule without running the actual calculation:

```bash
orca_descriptors approximate_time --molecule C1=CC=CC=C1
```

**Automatic Parameter Scaling**: The time estimation automatically scales benchmark data for different parameters (number of processors, functional, basis set). You don't need to re-run the benchmark if you change these parameters - the system will automatically recalculate the estimated time based on the existing benchmark data.

For example:
- If benchmark was run with 1 processor, estimation for 4 processors will automatically account for parallel efficiency
- If benchmark used `def2-SVP`, estimation for `def2-TZVP` will scale based on basis set size (O(N^3.5) scaling)
- Different functionals are scaled based on their relative computational costs

#### Available Parameters

All parameters from the `Orca` class are available as command-line arguments:

- `--script_path`: Path to ORCA executable (default: 'orca')
- `--working_dir`: Working directory for calculations (default: current directory)
- `--output_dir`: Directory for output files (default: current directory)
- `--functional`: DFT functional (default: PBE0)
- `--basis_set`: Basis set (default: def2-SVP)
- `--method_type`: Calculation type: Opt, SP, or Freq (default: Opt)
- `--dispersion_correction`: Dispersion correction, e.g., D3BJ (default: D3BJ). Use 'None' to disable.
- `--solvation_model`: Solvation model, e.g., 'COSMO(Water)' (default: None). Use 'None' to disable.
- `--n_processors`: Number of processors (default: 1)
- `--max_scf_cycles`: Maximum SCF cycles (default: 100)
- `--scf_convergence`: SCF convergence threshold (default: 1e-6)
- `--charge`: Molecular charge (default: 0)
- `--multiplicity`: Spin multiplicity (default: 1)
- `--cache_dir`: Directory for caching results (default: output_dir/.orca_cache)
- `--log_level`: Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)
- `--max_wait`: Maximum time to wait for output file creation in seconds (default: 300)

#### Example Commands

```bash
# Run benchmark with custom parameters (uses benzene as standard test molecule)
orca_descriptors run_benchmark \
    --functional PBE0 \
    --basis_set def2-SVP \
    --n_processors 4 \
    --working_dir ./calculations

# Estimate time for optimization calculation
orca_descriptors approximate_time \
    --molecule CCO \
    --method_type Opt \
    --n_opt_steps 20 \
    --functional PBE0 \
    --basis_set def2-TZVP \
    --n_processors 8
```

## Requirements

- Python >= 3.10
- ORCA 6.0.1 installed and available in PATH
- RDKit >= 2023.0.0

## License

See LICENSE.md

