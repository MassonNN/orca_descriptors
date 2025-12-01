"""Calculation execution methods for ORCA."""

import logging
import os
import re
import subprocess
from pathlib import Path
from typing import Any

from rdkit.Chem import Mol

logger = logging.getLogger(__name__)

if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
logger.propagate = False


class CalculationMixin:
    """Mixin class providing calculation execution methods for Orca."""
    
    def _build_command(self, orca_path: str, input_filename: str) -> list[str]:
        """Build command to run ORCA, optionally with mpirun.
        
        Args:
            orca_path: Path to ORCA executable
            input_filename: Name of input file
            
        Returns:
            Command list for subprocess
        """
        if self.use_mpirun:
            import shutil
            if self.mpirun_path:
                mpirun = self.mpirun_path
            else:
                mpirun = shutil.which("mpirun")
                if not mpirun:
                    raise RuntimeError(
                        "mpirun not found in PATH. Please specify mpirun_path or ensure mpirun is in PATH."
                    )
            cmd = [mpirun, "-np", str(self.n_processors), orca_path, input_filename]
        else:
            cmd = [orca_path, input_filename]
        return cmd
    
    def _build_environment(self) -> dict:
        """Build environment variables for ORCA process.
        
        Returns:
            Dictionary with environment variables
        """
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(self.n_processors)
        
        if self.extra_env:
            env.update(self.extra_env)
        
        return env
    
    def _check_orca_errors(self, content: str) -> list[str]:
        """Check for errors in ORCA output.
        
        Args:
            content: ORCA output content
            
        Returns:
            List of error messages found
        """
        errors = []
        lines = content.split('\n')
        
        error_patterns = [
            (r'INPUT ERROR', 'Input error detected'),
            (r'FATAL ERROR', 'Fatal error detected'),
            (r'ABORTING', 'Calculation aborted'),
            (r'TERMINATED ABNORMALLY', 'Calculation terminated abnormally'),
            (r'SCF NOT CONVERGED', 'SCF did not converge'),
            (r'GEOMETRY OPTIMIZATION FAILED', 'Geometry optimization failed'),
            (r'UNRECOGNIZED.*KEYWORD', 'Unrecognized keyword'),
            (r'DUPLICATED.*KEYWORD', 'Duplicated keyword'),
            (r'Unknown identifier', 'Unknown identifier in input block'),
            (r'Invalid assignment', 'Invalid assignment in input block'),
            (r'Cannot open', 'Cannot open file'),
            (r'File not found', 'File not found'),
            (r'ERROR.*finished by error', 'Error termination'),
            (r'ERROR.*aborting', 'Error aborting'),
            (r'ERROR.*termination', 'Error termination'),
            (r'ERROR(?!.*Last (MAX-Density|RMS-Density|DIIS Error|Orbital))', 'Error detected'),
        ]
        
        for i, line in enumerate(lines):
            line_upper = line.upper()
            
            if any(keyword in line_upper for keyword in [
                'LAST MAX-DENSITY', 'LAST RMS-DENSITY', 
                'LAST DIIS ERROR', 'LAST ORBITAL GRADIENT', 'LAST ORBITAL ROTATION',
                'TOLERANCE :'
            ]):
                continue
            
            if 'ERROR DETECTED' in line_upper or 'ERROR:' in line_upper:
                is_scf_info = False
                for j in range(i + 1, min(i + 5, len(lines))):
                    next_line_upper = lines[j].upper()
                    if any(keyword in next_line_upper for keyword in [
                        'LAST MAX-DENSITY', 'LAST RMS-DENSITY', 
                        'LAST DIIS ERROR', 'LAST ORBITAL', 'TOLERANCE'
                    ]):
                        is_scf_info = True
                        break
                if is_scf_info:
                    continue
            
            for pattern, error_type in error_patterns:
                if re.search(pattern, line, re.IGNORECASE):
                    context_start = max(0, i - 2)
                    context_end = min(len(lines), i + 3)
                    context = '\n'.join(lines[context_start:context_end])
                    errors.append(f"{error_type}:\n{context}")
                    break
        
        return errors
    
    def _run_calculation(self, mol: Mol) -> Path:
        """Run ORCA calculation and return path to output file.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Path to output file (cached if available)
        """
        if self.pre_optimize:
            mol = self._pre_optimize_geometry(mol)
        
        mol_hash = self._get_molecule_hash(mol)
        
        cached_output = self.cache.get(mol_hash)
        if cached_output and cached_output.exists():
            return cached_output
        
        input_content = self.input_generator.generate(
            mol=mol,
            functional=self.functional,
            basis_set=self.basis_set,
            method_type=self.method_type,
            dispersion_correction=self.dispersion_correction,
            solvation_model=self.solvation_model,
            n_processors=self.n_processors,
            max_scf_cycles=self.max_scf_cycles,
            scf_convergence=self.scf_convergence,
            charge=self.charge,
            multiplicity=self.multiplicity,
        )
        
        input_file = self.working_dir / f"orca_{mol_hash}.inp"
        base_name = f"orca_{mol_hash}"
        
        input_file.write_text(input_content)
        
        if os.path.isabs(self.script_path):
            orca_path = self.script_path
        else:
            import shutil
            orca_path = shutil.which(self.script_path)
            if not orca_path:
                raise RuntimeError(f"ORCA executable not found: {self.script_path}")
        
        input_filename = input_file.name
        cmd = self._build_command(orca_path, input_filename)
        env = self._build_environment()
        
        log_file_path = self.working_dir / f"orca_{mol_hash}.log"
        
        estimated_time = self.time_estimator.estimate_time(
            mol=mol,
            method_type=self.method_type,
            functional=self.functional,
            basis_set=self.basis_set,
            n_processors=self.n_processors,
        )
        
        if estimated_time > 0:
            hours = int(estimated_time // 3600)
            minutes = int((estimated_time % 3600) // 60)
            seconds = int(estimated_time % 60)
            if hours > 0:
                time_str = f"{hours}h {minutes}m {seconds}s"
            elif minutes > 0:
                time_str = f"{minutes}m {seconds}s"
            else:
                time_str = f"{seconds}s"
            logger.debug(f"Estimated calculation time: ~{time_str}")
        
        logger.debug("=" * 70)
        logger.debug(f"Running ORCA calculation: {input_filename}")
        logger.debug(f"Working directory: {self.working_dir}")
        logger.debug("=" * 70)
        
        with open(log_file_path, "w") as log_file:
            process = subprocess.Popen(
                cmd,
                cwd=str(self.working_dir),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                env=env,
                universal_newlines=True,
                bufsize=1,
                errors='ignore',
            )
            
            output_lines = []
            errors_detected = []
            
            try:
                for line in process.stdout:
                    if line:
                        output_lines.append(line)
                        line_upper = line.upper()
                        error_keywords = [
                            'INPUT ERROR', 'FATAL ERROR', 
                            'ABORTING', 'TERMINATED ABNORMALLY',
                            'SCF NOT CONVERGED', 'GEOMETRY OPTIMIZATION FAILED',
                            'UNKNOWN IDENTIFIER', 'UNRECOGNIZED', 'DUPLICATED',
                            'INVALID ASSIGNMENT'
                        ]
                        if 'ERROR' in line_upper:
                            if ('TERMINATED NORMALLY' not in line_upper and 
                                'NO ERROR' not in line_upper and
                                'LAST MAX-DENSITY' not in line_upper and
                                'LAST RMS-DENSITY' not in line_upper and
                                'LAST DIIS ERROR' not in line_upper):
                                errors_detected.append(line.strip())
                        
                        log_file.write(line)
                        log_file.flush()
                        
                        if logger.isEnabledFor(logging.DEBUG):
                            print(line, end='')
                
                process.wait()
                
                if process.returncode != 0:
                    error_summary = '\n'.join(errors_detected[:5])
                    raise RuntimeError(
                        f"ORCA calculation failed with return code {process.returncode}.\n"
                        f"Errors detected:\n{error_summary}"
                    )
                
                output_content = ''.join(output_lines)
                errors = self._check_orca_errors(output_content)
                
                if errors:
                    error_summary = '\n'.join(errors[:3])
                    raise RuntimeError(
                        f"ORCA calculation completed with errors:\n{error_summary}"
                    )
                
            except KeyboardInterrupt:
                process.terminate()
                process.wait()
                raise
            except Exception as e:
                if process.poll() is None:
                    process.terminate()
                    process.wait()
                raise
        
        possible_outputs = [
            self.working_dir / f"{base_name}.out",
            self.working_dir / f"{base_name}.log",
            self.working_dir / f"{base_name}.smd.out",
        ]
        
        output_file = None
        for possible_output in possible_outputs:
            if possible_output.exists():
                output_file = possible_output
                break
        
        if output_file is None:
            raise FileNotFoundError(
                f"ORCA output file not found. Expected one of: {[str(p) for p in possible_outputs]}"
            )
        
        cached_file = self.cache.store(mol_hash, output_file)
        
        self._cleanup_temp_files(base_name, output_file)
        
        if cached_file and cached_file.exists():
            return cached_file
        else:
            return output_file
    
    def _cleanup_temp_files(self, base_name: str, output_file: Path):
        """Remove all ORCA files after successful calculation.
        
        All files are removed since results are cached. This includes input files,
        output files, and all temporary files created by ORCA.
        
        Args:
            base_name: Base name for ORCA files (without extension)
            output_file: Path to the main output file (will be removed)
        """
        orca_extensions = [
            '.inp',
            '.out',
            '.log',
            '.smd.out',
            '.gbw',
            '.densities',
            '.densitiesinfo',
            '.ges',
            '.property.txt',
            '.bibtex',
            '.cpcm',
            '.cpcm_corr',
            '.engrad',
            '.opt',
            '.xyz',
            '_trj.xyz',
            '.molden',
            '.mkl',
            '.tmp',
            '.int.tmp',
        ]
        
        removed_count = 0
        removed_size = 0
        
        for ext in orca_extensions:
            patterns = [
                f"{base_name}{ext}",
                f"{base_name}*{ext}",
            ]
            
            for pattern in patterns:
                for orca_file in self.working_dir.glob(pattern):
                    if orca_file.is_file():
                        try:
                            file_size = orca_file.stat().st_size
                            orca_file.unlink()
                            removed_count += 1
                            removed_size += file_size
                            logger.debug(f"Removed ORCA file: {orca_file.name}")
                        except Exception as e:
                            logger.debug(f"Failed to remove ORCA file {orca_file.name}: {e}")
        
        if removed_count > 0:
            size_mb = removed_size / (1024 * 1024)
            logger.debug(f"Cleaned up {removed_count} ORCA files ({size_mb:.2f} MB)")
    
    def _get_output(self, mol: Mol) -> dict[str, Any]:
        """Get parsed output from ORCA calculation.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary with parsed ORCA output data
        """
        try:
            from orca_descriptors.batch_processing import XMolecule, DescriptorCall
            if isinstance(mol, XMolecule):
                raise ValueError("XMolecule cannot be used in _get_output")
        except ImportError:
            pass
        
        output_file = self._run_calculation(mol)
        return self.output_parser.parse(output_file, mol)

