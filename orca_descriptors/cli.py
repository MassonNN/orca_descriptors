"""Command-line interface for orca_descriptors."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from rdkit.Chem import MolFromSmiles, AddHs

from orca_descriptors import Orca


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser with all ORCA parameters."""
    parser = argparse.ArgumentParser(
        prog="orca_descriptors",
        description="ORCA descriptors library for QSAR analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Common arguments for all commands
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--script_path",
        type=str,
        default="orca",
        help="Path to ORCA executable (default: 'orca')",
    )
    common_parser.add_argument(
        "--working_dir",
        type=str,
        default=".",
        help="Working directory for ORCA calculations (default: current directory)",
    )
    common_parser.add_argument(
        "--output_dir",
        type=str,
        default=".",
        help="Directory for output files (default: current directory)",
    )
    common_parser.add_argument(
        "--functional",
        type=str,
        default="PBE0",
        help="DFT functional (default: PBE0)",
    )
    common_parser.add_argument(
        "--basis_set",
        type=str,
        default="def2-SVP",
        help="Basis set (default: def2-SVP)",
    )
    common_parser.add_argument(
        "--method_type",
        type=str,
        default="Opt",
        choices=["Opt", "SP", "Freq"],
        help="Calculation type: Opt (optimization), SP (single point), Freq (frequency) (default: Opt)",
    )
    common_parser.add_argument(
        "--dispersion_correction",
        type=str,
        default="D3BJ",
        help="Dispersion correction, e.g., D3BJ (default: D3BJ). Use 'None' to disable.",
    )
    common_parser.add_argument(
        "--solvation_model",
        type=str,
        default=None,
        help="Solvation model, e.g., 'COSMO(Water)' (default: None). Use 'None' to disable.",
    )
    common_parser.add_argument(
        "--n_processors",
        type=int,
        default=1,
        help="Number of processors (default: 1)",
    )
    common_parser.add_argument(
        "--max_scf_cycles",
        type=int,
        default=100,
        help="Maximum SCF cycles (default: 100)",
    )
    common_parser.add_argument(
        "--scf_convergence",
        type=float,
        default=1e-6,
        help="SCF convergence threshold (default: 1e-6)",
    )
    common_parser.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Molecular charge (default: 0)",
    )
    common_parser.add_argument(
        "--multiplicity",
        type=int,
        default=1,
        help="Spin multiplicity (default: 1)",
    )
    common_parser.add_argument(
        "--cache_dir",
        type=str,
        default=None,
        help="Directory for caching results (default: output_dir/.orca_cache)",
    )
    common_parser.add_argument(
        "--log_level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO)",
    )
    common_parser.add_argument(
        "--max_wait",
        type=int,
        default=300,
        help="Maximum time to wait for output file creation in seconds (default: 300)",
    )
    
    # run_benchmark command
    benchmark_parser = subparsers.add_parser(
        "run_benchmark",
        parents=[common_parser],
        help="Run ORCA benchmark calculation to calibrate time estimation. Uses benzene (C1=CC=CC=C1) as a standard test molecule.",
    )
    
    # approximate_time command
    time_parser = subparsers.add_parser(
        "approximate_time",
        parents=[common_parser],
        help="Estimate calculation time for a molecule without running the calculation",
    )
    time_parser.add_argument(
        "--molecule",
        type=str,
        required=True,
        help="SMILES string of the molecule to estimate time for",
    )
    time_parser.add_argument(
        "--n_opt_steps",
        type=int,
        default=None,
        help="Expected number of optimization steps (for Opt method, default: auto-estimate)",
    )
    
    return parser


def parse_log_level(level_str: str) -> int:
    """Convert log level string to logging constant."""
    level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    return level_map.get(level_str.upper(), logging.INFO)


def parse_dispersion_correction(value: Optional[str]) -> Optional[str]:
    """Parse dispersion correction value."""
    if value is None or value.lower() == "none":
        return None
    return value


def parse_solvation_model(value: Optional[str]) -> Optional[str]:
    """Parse solvation model value."""
    if value is None or value.lower() == "none":
        return None
    return value


def cmd_run_benchmark(args: argparse.Namespace) -> int:
    """Run benchmark calculation.
    
    Uses benzene (C1=CC=CC=C1) as a standard test molecule for machine calibration.
    """
    # Use benzene as standard benchmark molecule for machine calibration
    benchmark_smiles = "C1=CC=CC=C1"
    try:
        mol = MolFromSmiles(benchmark_smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {benchmark_smiles}")
        mol = AddHs(mol)
    except Exception as e:
        print(f"ERROR: Failed to parse benchmark molecule SMILES '{benchmark_smiles}': {e}", file=sys.stderr)
        return 1
    
    log_level = parse_log_level(args.log_level)
    dispersion = parse_dispersion_correction(args.dispersion_correction)
    solvation = parse_solvation_model(args.solvation_model)
    
    orca = Orca(
        script_path=args.script_path,
        working_dir=args.working_dir,
        output_dir=args.output_dir,
        functional=args.functional,
        basis_set=args.basis_set,
        method_type="SP",  # Benchmark uses single point
        dispersion_correction=dispersion,
        solvation_model=solvation,
        n_processors=args.n_processors,
        max_scf_cycles=args.max_scf_cycles,
        scf_convergence=args.scf_convergence,
        charge=args.charge,
        multiplicity=args.multiplicity,
        cache_dir=args.cache_dir,
        log_level=log_level,
        max_wait=args.max_wait,
    )
    
    try:
        benchmark_data = orca.run_benchmark(mol)
        
        print("\n" + "=" * 70)
        print("Benchmark Summary:")
        print("=" * 70)
        for key, value in benchmark_data.items():
            if isinstance(value, float):
                unit = " seconds" if "time" in key else ""
                print(f"{key.replace('_', ' ').capitalize()}: {value:.2f}{unit}")
            else:
                print(f"{key.replace('_', ' ').capitalize()}: {value}")
        print(f"\nBenchmark data saved to: {orca.time_estimator.benchmark_file}")
        print("=" * 70 + "\n")
        
        return 0
    except Exception as e:
        print(f"ERROR: Benchmark calculation failed: {e}", file=sys.stderr)
        return 1


def cmd_approximate_time(args: argparse.Namespace) -> int:
    """Estimate calculation time without running the calculation."""
    try:
        mol = MolFromSmiles(args.molecule)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {args.molecule}")
        mol = AddHs(mol)
    except Exception as e:
        print(f"ERROR: Failed to parse molecule SMILES '{args.molecule}': {e}", file=sys.stderr)
        return 1
    
    log_level = parse_log_level(args.log_level)
    dispersion = parse_dispersion_correction(args.dispersion_correction)
    solvation = parse_solvation_model(args.solvation_model)
    
    orca = Orca(
        script_path=args.script_path,
        working_dir=args.working_dir,
        output_dir=args.output_dir,
        functional=args.functional,
        basis_set=args.basis_set,
        method_type=args.method_type,
        dispersion_correction=dispersion,
        solvation_model=solvation,
        n_processors=args.n_processors,
        max_scf_cycles=args.max_scf_cycles,
        scf_convergence=args.scf_convergence,
        charge=args.charge,
        multiplicity=args.multiplicity,
        cache_dir=args.cache_dir,
        log_level=log_level,
        max_wait=args.max_wait,
    )
    
    try:
        estimated_time = orca.estimate_calculation_time(mol, n_opt_steps=args.n_opt_steps)
        
        if estimated_time > 0:
            hours = int(estimated_time // 3600)
            minutes = int((estimated_time % 3600) // 60)
            seconds = estimated_time % 60
            
            print("\n" + "=" * 70)
            print("Estimated Calculation Time:")
            print("=" * 70)
            print(f"Total: {estimated_time:.2f} seconds")
            if hours > 0:
                print(f"  ({hours} hours, {minutes} minutes, {seconds:.1f} seconds)")
            elif minutes > 0:
                print(f"  ({minutes} minutes, {seconds:.1f} seconds)")
            else:
                print(f"  ({seconds:.1f} seconds)")
            print("=" * 70 + "\n")
            return 0
        else:
            print("WARNING: Could not estimate calculation time.", file=sys.stderr)
            print("Please run 'orca_descriptors run_benchmark' first.", file=sys.stderr)
            return 1
    except Exception as e:
        print(f"ERROR: Time estimation failed: {e}", file=sys.stderr)
        return 1


def main() -> int:
    """Main entry point for CLI."""
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    if args.command == "run_benchmark":
        return cmd_run_benchmark(args)
    elif args.command == "approximate_time":
        return cmd_approximate_time(args)
    else:
        print(f"ERROR: Unknown command: {args.command}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

