"""Time estimation for ORCA calculations."""

import json
import logging
from pathlib import Path
from typing import Optional

from rdkit.Chem import Mol

from orca_descriptors.output_parser import ORCAOutputParser

logger = logging.getLogger(__name__)


class ORCATimeEstimator:
    """Estimate ORCA calculation time based on benchmark data."""
    
    def __init__(self, working_dir: Path):
        """Initialize time estimator.
        
        Args:
            working_dir: Working directory where benchmark data is stored
        """
        self.working_dir = Path(working_dir)
        self.benchmark_file = self.working_dir / ".orca_benchmark.json"
        self.parser = ORCAOutputParser()
    
    def run_benchmark(
        self,
        mol: Mol,
        functional: str = "PBE0",
        basis_set: str = "def2-SVP",
        script_path: str = "orca",
        n_processors: int = 1,
        use_mpirun: bool = False,
        mpirun_path: Optional[str] = None,
        extra_env: Optional[dict] = None,
    ) -> dict:
        """Run benchmark calculation to calibrate time estimation.
        
        Args:
            mol: Test molecule (should be small, e.g., benzene)
            functional: DFT functional to use
            basis_set: Basis set to use
            script_path: Path to ORCA executable
            n_processors: Number of processors
            use_mpirun: Whether to use mpirun for parallel execution (default: False)
            mpirun_path: Path to mpirun executable (default: None, will search in PATH)
            extra_env: Additional environment variables to pass to ORCA process (default: None)
            
        Returns:
            Dictionary with benchmark data
        """
        import subprocess
        import os
        import time
        import hashlib
        from rdkit.Chem import MolToXYZBlock, AllChem
        
        logger.info("=" * 70)
        logger.info("Running ORCA benchmark calculation...")
        logger.info("=" * 70)
        
        # Generate 3D coordinates if needed
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol)
        
        # Generate input file for SP calculation
        mol_hash = hashlib.sha256(MolToXYZBlock(mol).encode()).hexdigest()
        input_file = self.working_dir / f"benchmark_{mol_hash}.inp"
        
        # Create SP input (single point, no optimization)
        input_lines = [
            f"! SP {functional} {basis_set}",
            "",
            "%pal",
            f"  nprocs {n_processors}",
            "end",
            "",
            "%scf",
            "  MaxIter 100",
            "  ConvForced true",
            "end",
            "",
            "%output",
            "  PrintLevel 2",
            "end",
            "",
            "* xyz 0 1",
        ]
        
        # Add coordinates
        xyz_block = MolToXYZBlock(mol)
        xyz_lines = xyz_block.strip().split("\n")
        if len(xyz_lines) > 2:
            for line in xyz_lines[2:]:
                input_lines.append(line)
        
        input_lines.append("*")
        input_file.write_text("\n".join(input_lines))
        
        # Run ORCA
        if not os.path.isabs(script_path):
            import shutil
            orca_path = shutil.which(script_path)
            if not orca_path:
                raise RuntimeError(f"ORCA executable not found: {script_path}")
        else:
            orca_path = script_path
        
        # Build command with optional mpirun
        if use_mpirun:
            import shutil
            if mpirun_path:
                mpirun = mpirun_path
            else:
                mpirun = shutil.which("mpirun")
                if not mpirun:
                    raise RuntimeError(
                        "mpirun not found in PATH. Please specify mpirun_path or ensure mpirun is in PATH."
                    )
            cmd = [mpirun, "-np", str(n_processors), orca_path, input_file.name]
        else:
            cmd = [orca_path, input_file.name]
        
        # Build environment
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(n_processors)
        if extra_env:
            env.update(extra_env)
        
        # Define base_name before using it
        base_name = f"benchmark_{mol_hash}"
        
        start_time = time.time()
        
        # Save stdout to log file (ORCA may not create .out file for SP calculations)
        log_file_path = self.working_dir / f"{base_name}.log"
        
        with open(log_file_path, "w") as log_file:
            process = subprocess.run(
                cmd,
                cwd=str(self.working_dir),
                stdout=log_file,
                stderr=subprocess.STDOUT,
                env=env,
                universal_newlines=True,
            )
        
        total_time = time.time() - start_time
        
        # Use the log file we created, or find ORCA's output file
        output_file = None
        
        # First, check if our log file exists and has content
        log_file_path = self.working_dir / f"{base_name}.log"
        if log_file_path.exists() and log_file_path.stat().st_size > 0:
            output_file = log_file_path
        else:
            # Try to find ORCA's output files
            possible_outputs = [
                self.working_dir / f"{base_name}.out",
                self.working_dir / f"{base_name}.log",
                self.working_dir / f"{base_name}.smd.out",
            ]
            
            for possible_output in possible_outputs:
                if possible_output.exists():
                    output_file = possible_output
                    break
        
        if not output_file:
            raise RuntimeError(
                f"Benchmark output file not found for {base_name}.\n"
                f"Return code: {process.returncode}\n"
                f"Expected log file: {log_file_path}\n"
                f"Log file exists: {log_file_path.exists() if log_file_path else False}"
            )
        
        # Parse output to get benchmark data
        content = output_file.read_text()
        
        # Extract number of basis functions
        n_basis = self._extract_n_basis(content)
        
        # Extract SCF cycle time
        scf_time = self._extract_scf_time(content, total_time)
        
        # Extract number of atoms
        n_atoms = mol.GetNumAtoms()
        
        benchmark_data = {
            "functional": functional,
            "basis_set": basis_set,
            "n_processors": n_processors,
            "n_atoms": n_atoms,
            "n_basis": n_basis,
            "scf_time": scf_time,
            "total_time": total_time,
            "n_scf_cycles": self._extract_n_scf_cycles(content),
        }
        
        # Save benchmark data
        self._save_benchmark(benchmark_data)
        
        # Clean up temporary files
        self._cleanup_temp_files(base_name, output_file)
        
        logger.info("=" * 70)
        logger.info("Benchmark completed successfully!")
        logger.info(f"Number of basis functions: {n_basis}")
        logger.info(f"SCF cycle time: {scf_time:.2f} seconds")
        logger.info(f"Total calculation time: {total_time:.2f} seconds")
        logger.info("=" * 70)
        
        return benchmark_data
    
    def _cleanup_temp_files(self, base_name: str, output_file: Path):
        """Remove temporary ORCA files after successful calculation.
        
        Args:
            base_name: Base name for ORCA files (without extension)
            output_file: Path to the main output file to keep
        """
        # Files to keep (main results)
        files_to_keep = {
            output_file.name,
            f"{base_name}.inp",  # Keep input file for reference
        }
        
        # Also keep alternative output files if they exist
        for ext in ['.out', '.log', '.smd.out']:
            alt_file = self.working_dir / f"{base_name}{ext}"
            if alt_file.exists():
                files_to_keep.add(alt_file.name)
        
        # Temporary file extensions to remove
        temp_extensions = [
            '.gbw',      # Wavefunction file
            '.densities', # Density files
            '.densitiesinfo',
            '.ges',      # Geometry file
            '.property.txt',  # Property file
            '.bibtex',   # Bibliography
            '.cpcm',     # CPCM files
            '.cpcm_corr',
            '.engrad',   # Energy gradient
            '.opt',      # Optimization file
            '.xyz',      # XYZ trajectory (if not needed)
            '_trj.xyz',  # Trajectory file
            '.molden',  # Molden file
            '.mkl',     # MKL file
            '.tmp',     # Temporary files
            '.int.tmp', # Integral temporary files
        ]
        
        removed_count = 0
        removed_size = 0
        
        # Remove temporary files
        for ext in temp_extensions:
            # Handle both patterns: base_name.ext and base_name_pattern.ext
            patterns = [
                f"{base_name}{ext}",
                f"{base_name}*{ext}",
            ]
            
            for pattern in patterns:
                for temp_file in self.working_dir.glob(pattern):
                    if temp_file.name not in files_to_keep and temp_file.is_file():
                        try:
                            file_size = temp_file.stat().st_size
                            temp_file.unlink()
                            removed_count += 1
                            removed_size += file_size
                            logger.debug(f"Removed temporary file: {temp_file.name}")
                        except Exception as e:
                            logger.warning(f"Failed to remove temporary file {temp_file.name}: {e}")
        
        if removed_count > 0:
            size_mb = removed_size / (1024 * 1024)
            logger.info(f"Cleaned up {removed_count} temporary files ({size_mb:.2f} MB)")
    
    def estimate_time(
        self,
        mol: Mol,
        method_type: str = "Opt",
        functional: str = "PBE0",
        basis_set: str = "def2-SVP",
        n_processors: int = 1,
        n_opt_steps: Optional[int] = None,
    ) -> float:
        """Estimate calculation time for a molecule.
        
        Automatically scales benchmark data if parameters differ (e.g., different
        number of processors). This allows estimation without re-running benchmark.
        
        Args:
            mol: Target molecule
            method_type: Calculation type ("Opt", "Freq", "SP")
            functional: DFT functional
            basis_set: Basis set
            n_processors: Number of processors
            n_opt_steps: Expected number of optimization steps (for Opt)
            
        Returns:
            Estimated time in seconds
        """
        benchmark = self._load_benchmark()
        
        if not benchmark:
            logger.warning("No benchmark data found. Run benchmark first.")
            return 0.0
        
        # Track which parameters differ
        params_differ = []
        if benchmark["functional"] != functional:
            params_differ.append(f"functional ({benchmark['functional']} -> {functional})")
        if benchmark["basis_set"] != basis_set:
            params_differ.append(f"basis_set ({benchmark['basis_set']} -> {basis_set})")
        if benchmark["n_processors"] != n_processors:
            params_differ.append(f"n_processors ({benchmark['n_processors']} -> {n_processors})")
        
        # If parameters differ, apply scaling factors
        scaling_factor = 1.0
        
        if params_differ:
            logger.info(
                f"Benchmark parameters differ. Scaling estimation: {', '.join(params_differ)}"
            )
            
            # Scale for different number of processors
            # Parallel efficiency: typically 0.7-0.9 for good scaling
            # Time scales approximately as 1/n_proc with efficiency factor
            if benchmark["n_processors"] != n_processors:
                n_proc_old = benchmark["n_processors"]
                n_proc_new = n_processors
                
                # Parallel efficiency factor (accounts for communication overhead)
                # Better scaling for more processors, but diminishing returns
                if n_proc_new > n_proc_old:
                    # Scaling up: efficiency decreases
                    efficiency = 0.85 - 0.05 * min(4, (n_proc_new / n_proc_old - 1))
                else:
                    # Scaling down: efficiency improves
                    efficiency = 0.90 + 0.05 * min(2, (n_proc_old / n_proc_new - 1))
                
                proc_scaling = (n_proc_old / n_proc_new) * efficiency
                scaling_factor *= proc_scaling
                logger.info(
                    f"Processor scaling: {n_proc_old} -> {n_proc_new} processors "
                    f"(efficiency: {efficiency:.2f}, factor: {proc_scaling:.3f})"
                )
            
            # Scale for different functional (rough estimates)
            # Different functionals have different computational costs
            functional_costs = {
                "PBE0": 1.0,      # Reference
                "B3LYP": 1.0,     # Similar to PBE0
                "M06": 1.2,       # Slightly more expensive
                "M06-2X": 1.3,
                "wB97X-D": 1.4,   # More expensive
                "CAM-B3LYP": 1.1,
                "PBE": 0.9,       # Slightly cheaper (no exact exchange)
                "BLYP": 0.9,
            }
            
            if benchmark["functional"] != functional:
                cost_old = functional_costs.get(benchmark["functional"], 1.0)
                cost_new = functional_costs.get(functional, 1.0)
                func_scaling = cost_new / cost_old
                scaling_factor *= func_scaling
                logger.info(
                    f"Functional scaling: {benchmark['functional']} -> {functional} "
                    f"(factor: {func_scaling:.3f})"
                )
            
            # Scale for different basis set
            # Larger basis sets are more expensive (roughly O(N^3.5) scaling)
            basis_sizes = {
                "STO-3G": 0.3,
                "3-21G": 0.5,
                "6-31G": 0.7,
                "6-31G*": 0.8,
                "6-31G**": 0.9,
                "def2-SVP": 1.0,      # Reference
                "def2-TZVP": 2.5,
                "def2-TZVPP": 3.0,
                "def2-QZVP": 5.0,
                "cc-pVDZ": 1.0,
                "cc-pVTZ": 2.8,
                "cc-pVQZ": 5.5,
            }
            
            if benchmark["basis_set"] != basis_set:
                size_old = basis_sizes.get(benchmark["basis_set"], 1.0)
                size_new = basis_sizes.get(basis_set, 1.0)
                # Use power law scaling for basis set size
                basis_scaling = (size_new / size_old) ** 3.5
                scaling_factor *= basis_scaling
                logger.info(
                    f"Basis set scaling: {benchmark['basis_set']} -> {basis_set} "
                    f"(factor: {basis_scaling:.3f})"
                )
        
        # Get target molecule size
        n_atoms_target = mol.GetNumAtoms()
        n_atoms_benchmark = benchmark["n_atoms"]
        n_basis_benchmark = benchmark["n_basis"]
        
        # Estimate basis functions for target molecule
        atoms_per_basis = n_basis_benchmark / n_atoms_benchmark if n_atoms_benchmark > 0 else 1
        n_basis_target = int(n_atoms_target * atoms_per_basis)
        
        # Use more realistic scaling: O(N^2.5) for basis functions
        # This is more accurate for typical DFT calculations
        # The exponent 2.5 accounts for:
        # - O(N^2) for matrix operations
        # - O(N^0.5) overhead for larger systems
        scaling_exponent = 2.5
        basis_ratio = (n_basis_target / n_basis_benchmark) ** scaling_exponent if n_basis_benchmark > 0 else 1.0
        
        # Apply scaling factor for different parameters (processors, functional, basis set)
        size_scaling = basis_ratio * scaling_factor
        
        # Estimate total time based on method type
        # Use benchmark total_time as base, scaled by molecule size
        if method_type == "SP":
            # Single point: scale benchmark total_time by size
            estimated_time = benchmark["total_time"] * size_scaling
        elif method_type == "Opt":
            # Optimization: scale benchmark total_time, but account for optimization overhead
            # Opt typically takes 2-3x longer than SP for the same molecule
            # But benchmark is SP, so we need to account for Opt overhead
            opt_overhead = 2.5  # Opt is typically 2-3x slower than SP
            estimated_time = benchmark["total_time"] * size_scaling * opt_overhead
            
            # Optionally adjust for expected number of optimization steps
            # But in practice, the size scaling already accounts for this
            if n_opt_steps is not None:
                # If user provides n_opt_steps, we can fine-tune
                # But typically the size scaling is sufficient
                pass
        elif method_type == "Freq":
            # Frequency calculation: much more expensive
            # Typically 10-20x slower than SP due to Hessian calculation
            freq_overhead = 15.0
            estimated_time = benchmark["total_time"] * size_scaling * freq_overhead
        else:
            # Default: assume similar to Opt
            opt_overhead = 2.5
            estimated_time = benchmark["total_time"] * size_scaling * opt_overhead
        
        if params_differ:
            logger.info(
                f"Estimated time with parameter scaling (total scaling factor: {scaling_factor:.3f}, "
                f"size scaling: {basis_ratio:.3f})"
            )
        
        return estimated_time
    
    def _extract_n_basis(self, content: str) -> int:
        """Extract number of basis functions from ORCA output."""
        import re
        
        # Look for "Number of basis functions" or similar
        patterns = [
            r"Number of basis functions\s*\.\.\.\s*(\d+)",
            r"Basis functions\s*:\s*(\d+)",
            r"NBasis\s*=\s*(\d+)",
            r"Number of contracted basis functions\s*:\s*(\d+)",
            r"Total number of basis functions\s*:\s*(\d+)",
            r"Number of contracted shells\s*:\s*(\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return int(match.group(1))
        
        # Try to find in summary section
        # Look for lines like "Basis: def2-SVP (XX functions)"
        match = re.search(r"Basis:.*?\((\d+)\s+functions?\)", content, re.IGNORECASE)
        if match:
            return int(match.group(1))
        
        # Fallback: estimate from atoms (rough approximation)
        # For def2-SVP: ~15-20 functions per heavy atom, ~1 per H
        logger.warning("Could not extract number of basis functions from output")
        return 0
    
    def _extract_scf_time(self, content: str, total_time: float) -> float:
        """Extract average SCF cycle time from ORCA output."""
        import re
        
        # Look for SCF timing information in timings section
        # Format: "SCF iterations                   ...       66.924 sec (=   1.115 min)  52.4 %"
        patterns = [
            r"SCF iterations\s+\.\.\.\s+([\d.]+)\s+sec",
            r"SCF\s+\.\.\.\s+([\d.]+)\s+sec",
            r"SCF iterations\s*:\s*([\d.]+)\s*sec",
        ]
        
        scf_time_total = None
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                scf_time_total = float(match.group(1))
                break
        
        # Extract number of SCF cycles
        n_cycles = self._extract_n_scf_cycles(content)
        
        if scf_time_total and n_cycles > 0:
            return scf_time_total / n_cycles
        
        # Fallback: estimate from total time
        if n_cycles > 0:
            # Assume SCF takes ~80% of total time for SP calculation
            return (total_time * 0.8) / n_cycles
        
        # Last resort: estimate
        return total_time / 20  # Assume 20 cycles
    
    def _extract_n_scf_cycles(self, content: str) -> int:
        """Extract number of SCF cycles from ORCA output."""
        import re
        
        patterns = [
            r"SCF\s+ITERATIONS\s*:\s*(\d+)",
            r"Total number of iterations\s*:\s*(\d+)",
            r"Number of SCF iterations\s*:\s*(\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return int(match.group(1))
        
        # Count SCF iteration lines
        iterations = re.findall(r"ITERATION\s+(\d+)", content, re.IGNORECASE)
        if iterations:
            return len(set(iterations))
        
        return 20  # Default estimate
    
    def _save_benchmark(self, data: dict):
        """Save benchmark data to file."""
        self.working_dir.mkdir(parents=True, exist_ok=True)
        with open(self.benchmark_file, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Benchmark data saved to {self.benchmark_file}")
    
    def _load_benchmark(self) -> Optional[dict]:
        """Load benchmark data from file."""
        if not self.benchmark_file.exists():
            return None
        
        try:
            with open(self.benchmark_file, "r") as f:
                return json.load(f)
        except Exception as e:
            logger.warning(f"Failed to load benchmark data: {e}")
            return None

