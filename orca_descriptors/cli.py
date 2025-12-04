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
        default="AM1",
        help="DFT functional or semi-empirical method (default: AM1)",
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
        "--cache_server_url",
        type=str,
        default=None,
        help="URL of the remote cache server (e.g., http://localhost:3000). Requires --cache_api_token.",
    )
    common_parser.add_argument(
        "--cache_api_token",
        type=str,
        default=None,
        help="API token for remote cache authentication. Requires --cache_server_url.",
    )
    common_parser.add_argument(
        "--cache_timeout",
        type=int,
        default=30,
        help="Timeout for remote cache requests in seconds (default: 30)",
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
    common_parser.add_argument(
        "--use_mpirun",
        action="store_true",
        help="Use mpirun for parallel execution (default: False)",
    )
    common_parser.add_argument(
        "--mpirun_path",
        type=str,
        default=None,
        help="Path to mpirun executable (default: None, will search in PATH)",
    )
    common_parser.add_argument(
        "--extra_env",
        type=str,
        nargs="*",
        default=None,
        metavar="KEY=VALUE",
        help="Additional environment variables to pass to ORCA process (format: KEY=VALUE). Can be specified multiple times.",
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
    
    # clear command
    clear_parser = subparsers.add_parser(
        "clear",
        parents=[common_parser],
        help="Remove all ORCA files in the working directory (if not removed due to error)",
    )
    
    # purge_cache command
    purge_parser = subparsers.add_parser(
        "purge_cache",
        parents=[common_parser],
        help="Remove ORCA cache",
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


def parse_extra_env(env_list: Optional[list[str]]) -> Optional[dict]:
    """Parse extra environment variables from list of KEY=VALUE strings.
    
    Args:
        env_list: List of strings in format "KEY=VALUE"
        
    Returns:
        Dictionary with environment variables, or None if env_list is None/empty
    """
    if not env_list:
        return None
    
    env_dict = {}
    for env_var in env_list:
        if "=" not in env_var:
            raise ValueError(f"Invalid environment variable format: {env_var}. Expected KEY=VALUE")
        key, value = env_var.split("=", 1)
        env_dict[key] = value
    
    return env_dict


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
    extra_env = parse_extra_env(args.extra_env)
    
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
        use_mpirun=args.use_mpirun,
        mpirun_path=args.mpirun_path,
        extra_env=extra_env,
        cache_server_url=args.cache_server_url,
        cache_api_token=args.cache_api_token,
        cache_timeout=args.cache_timeout,
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
    extra_env = parse_extra_env(args.extra_env)
    
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
        use_mpirun=args.use_mpirun,
        mpirun_path=args.mpirun_path,
        extra_env=extra_env,
        cache_server_url=args.cache_server_url,
        cache_api_token=args.cache_api_token,
        cache_timeout=args.cache_timeout,
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


def cmd_clear(args: argparse.Namespace) -> int:
    """Remove all ORCA files in the working directory."""
    working_dir = Path(args.working_dir)
    if not working_dir.exists():
        print(f"ERROR: Working directory does not exist: {working_dir}", file=sys.stderr)
        return 1
    
    orca_extensions = [
        ".inp",
        ".out",
        ".log",
        ".smd.out",
        ".gbw",
        ".densities",
        ".densitiesinfo",
        ".ges",
        ".property.txt",
        ".bibtex",
        ".cpcm",
        ".cpcm_corr",
        ".engrad",
        ".opt",
        ".xyz",
        "_trj.xyz",
        ".molden",
        ".mkl",
        ".tmp",
        ".int.tmp",
        ".cube",
        ".prop",
        ".prop.scf",
    ]
    
    removed_count = 0
    removed_size = 0
    errors = []
    processed_files = set()
    
    # First, remove all files starting with "orca_" (covers orca_<hash>.tmp0, orca_<hash>.tmp1, etc.)
    for file_path in working_dir.glob("orca_*"):
        try:
            if file_path.is_file() and file_path not in processed_files:
                size = file_path.stat().st_size
                file_path.unlink()
                removed_count += 1
                removed_size += size
                processed_files.add(file_path)
        except Exception as e:
            errors.append(f"Failed to remove {file_path}: {e}")
    
    # Then, remove files by extension patterns
    for ext in orca_extensions:
        patterns = [
            f"orca_*{ext}",
            f"*{ext}",
        ]
        
        for pattern in patterns:
            for file_path in working_dir.glob(pattern):
                try:
                    if file_path.is_file() and file_path not in processed_files:
                        size = file_path.stat().st_size
                        file_path.unlink()
                        removed_count += 1
                        removed_size += size
                        processed_files.add(file_path)
                except Exception as e:
                    errors.append(f"Failed to remove {file_path}: {e}")
    
    print("\n" + "=" * 70)
    print("ORCA Files Cleanup Summary:")
    print("=" * 70)
    print(f"Files removed: {removed_count}")
    print(f"Total size freed: {removed_size / (1024 * 1024):.2f} MB")
    
    if errors:
        print(f"\nErrors encountered: {len(errors)}")
        for error in errors[:10]:
            print(f"  {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more errors")
    
    print("=" * 70 + "\n")
    
    return 0 if not errors else 1


def cmd_purge_cache(args: argparse.Namespace) -> int:
    """Remove ORCA cache."""
    from orca_descriptors.cache import CacheManager
    
    if args.cache_dir:
        cache_dir = Path(args.cache_dir)
    else:
        output_dir = Path(args.output_dir)
        cache_dir = output_dir / ".orca_cache"
    
    if not cache_dir.exists():
        print(f"Cache directory does not exist: {cache_dir}")
        return 0
    
    try:
        cache = CacheManager(cache_dir)
        cache.clear()
        print(f"\nCache cleared successfully: {cache_dir}\n")
        return 0
    except Exception as e:
        print(f"ERROR: Failed to clear cache: {e}", file=sys.stderr)
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
    elif args.command == "clear":
        return cmd_clear(args)
    elif args.command == "purge_cache":
        return cmd_purge_cache(args)
    else:
        print(f"ERROR: Unknown command: {args.command}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

