#!/usr/bin/env python3
"""Script to run ORCA benchmark for time estimation."""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from orca_descriptors import Orca
from rdkit.Chem import MolFromSmiles, AddHs
import logging

logging.basicConfig(level=logging.INFO)

def main():
    """Run benchmark calculation."""
    # Initialize ORCA with default settings
    orca = Orca(
        script_path="orca",
        working_dir=".",
        functional="PBE0",
        basis_set="def2-SVP",
        n_processors=1,
    )
    
    # Use benzene as benchmark molecule (or provide custom molecule)
    if len(sys.argv) > 1:
        smiles = sys.argv[1]
        print(f"Running benchmark with molecule: {smiles}")
        mol = MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        mol = AddHs(mol)
    else:
        print("Running benchmark with default molecule: benzene (C1=CC=CC=C1)")
        mol = None  # Will use default benzene
    
    # Run benchmark
    benchmark_data = orca.run_benchmark(mol)
    
    print("\n" + "=" * 70)
    print("Benchmark Summary:")
    print("=" * 70)
    print(f"Functional: {benchmark_data['functional']}")
    print(f"Basis set: {benchmark_data['basis_set']}")
    print(f"Number of processors: {benchmark_data['n_processors']}")
    print(f"Number of atoms: {benchmark_data['n_atoms']}")
    print(f"Number of basis functions: {benchmark_data['n_basis']}")
    print(f"SCF cycle time: {benchmark_data['scf_time']:.2f} seconds")
    print(f"Number of SCF cycles: {benchmark_data['n_scf_cycles']}")
    print(f"Total calculation time: {benchmark_data['total_time']:.2f} seconds")
    print(f"\nBenchmark data saved to: {orca.time_estimator.benchmark_file}")
    print("=" * 70)

if __name__ == "__main__":
    main()

