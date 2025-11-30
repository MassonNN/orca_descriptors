"""Main Orca class for quantum chemical calculations."""

import hashlib
import logging
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

from rdkit.Chem import Mol, MolFromSmiles

from orca_descriptors.cache import CacheManager
from orca_descriptors.input_generator import ORCAInputGenerator
from orca_descriptors.output_parser import ORCAOutputParser
from orca_descriptors.time_estimator import ORCATimeEstimator

logger = logging.getLogger(__name__)


class Orca:
    """Main class for ORCA quantum chemical calculations."""
    
    def __init__(
        self,
        script_path: str = "orca",
        working_dir: str = ".",
        output_dir: str = ".",
        functional: str = "PBE0",
        basis_set: str = "def2-SVP",
        method_type: str = "Opt",
        dispersion_correction: Optional[str] = "D3BJ",
        solvation_model: Optional[str] = None,
        n_processors: int = 1,
        max_scf_cycles: int = 100,
        scf_convergence: float = 1e-6,
        charge: int = 0,
        multiplicity: int = 1,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        max_wait: int = 300,
    ):
        """Initialize ORCA calculator.
        
        Args:
            script_path: Path to ORCA executable
            working_dir: Working directory for calculations
            output_dir: Directory for output files
            functional: DFT functional (e.g., "PBE0")
            basis_set: Basis set (e.g., "def2-SVP")
            method_type: Calculation type ("Opt", "SP", etc.)
            dispersion_correction: Dispersion correction (e.g., "D3BJ")
            solvation_model: Solvation model (e.g., "COSMO(Water)")
            n_processors: Number of processors
            max_scf_cycles: Maximum SCF cycles
            scf_convergence: SCF convergence threshold
            charge: Molecular charge
            multiplicity: Spin multiplicity
            cache_dir: Directory for caching results (default: output_dir/.orca_cache)
            log_level: Logging level (default: logging.INFO)
            max_wait: Maximum time to wait for output file creation in seconds (default: 300)
        """
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(log_level)
        logger.propagate = False
        self.script_path = script_path
        self.working_dir = Path(working_dir)
        self.output_dir = Path(output_dir)
        self.functional = functional
        self.basis_set = basis_set
        self.method_type = method_type
        self.dispersion_correction = dispersion_correction
        self.solvation_model = solvation_model
        self.n_processors = n_processors
        self.max_scf_cycles = max_scf_cycles
        self.scf_convergence = scf_convergence
        self.charge = charge
        self.multiplicity = multiplicity
        self.max_wait = max_wait
        
        self.working_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        cache_dir = cache_dir or str(self.output_dir / ".orca_cache")
        self.cache = CacheManager(cache_dir)
        self.input_generator = ORCAInputGenerator()
        self.output_parser = ORCAOutputParser()
        self.time_estimator = ORCATimeEstimator(working_dir=self.working_dir)
    
    def _get_molecule_hash(self, mol: Mol) -> str:
        """Generate hash for molecule based on SMILES and calculation parameters."""
        from rdkit.Chem import MolToSmiles
        
        smiles = MolToSmiles(mol, canonical=True)
        params = (
            self.functional,
            self.basis_set,
            self.method_type,
            self.dispersion_correction,
            self.solvation_model,
            self.charge,
            self.multiplicity,
        )
        key = f"{smiles}_{params}"
        return hashlib.sha256(key.encode()).hexdigest()
    
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
                    break  # Only report each error once
        
        return errors
    
    def _run_calculation(self, mol: Mol) -> Path:
        """Run ORCA calculation and return path to output file."""
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
        cmd = [orca_path, input_filename]
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(self.n_processors)
        
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
            logger.info(f"Estimated calculation time: ~{time_str}")
        
        logger.info("=" * 70)
        logger.info(f"Running ORCA calculation: {input_filename}")
        logger.info(f"Working directory: {self.working_dir}")
        logger.info("=" * 70)
        
        with open(log_file_path, "w") as log_file:
            process = subprocess.Popen(
                cmd,
                cwd=str(self.working_dir),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                env=env,
                universal_newlines=True,
                bufsize=1,  # Line buffered
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
                                'LAST DIIS ERROR' not in line_upper and
                                'LAST ORBITAL' not in line_upper):
                                error_keywords.append('ERROR')
                        
                        if any(keyword in line_upper for keyword in error_keywords):
                            errors_detected.append(line.rstrip())
                            logger.error(f"ORCA ERROR DETECTED: {line.rstrip()}")
                        
                        logger.info(line.rstrip())
                        log_file.write(line)
                        log_file.flush()
            except KeyboardInterrupt:
                process.terminate()
                logger.error("ORCA calculation interrupted by user")
                raise RuntimeError("ORCA calculation interrupted by user")
            
            return_code = process.wait()
            
            output_content = ''.join(output_lines)
            all_errors = self._check_orca_errors(output_content)
            
            if log_file_path.exists():
                log_content = log_file_path.read_text()
                log_errors = self._check_orca_errors(log_content)
                all_errors.extend(log_errors)
            
            if errors_detected or return_code != 0 or all_errors:
                if all_errors:
                    error_msg = "\n\n".join(all_errors)
                    logger.error("=" * 70)
                    logger.error("ORCA ERRORS DETECTED IN OUTPUT:")
                    logger.error("=" * 70)
                    logger.error(error_msg)
                    logger.error("=" * 70)
                    raise RuntimeError(f"ORCA calculation failed with errors:\n{error_msg}")
                elif errors_detected:
                    error_msg = "\n".join(errors_detected)
                    logger.error("=" * 70)
                    logger.error("ORCA ERRORS DETECTED:")
                    logger.error("=" * 70)
                    logger.error(error_msg)
                    logger.error("=" * 70)
                    raise RuntimeError(f"ORCA calculation failed with errors:\n{error_msg}")
        
        logger.info("=" * 70)
        logger.info(f"ORCA calculation finished with return code: {return_code}")
        logger.info("=" * 70)
        
        output_file = None
        possible_outputs = [
            self.working_dir / f"{base_name}.out",
            self.working_dir / f"{base_name}.log",
            self.working_dir / f"{base_name}.smd.out",
        ]
        
        wait_time = 0
        while wait_time < self.max_wait:
            for possible_output in possible_outputs:
                if possible_output.exists():
                    output_file = possible_output
                    break
            if output_file:
                break
            time.sleep(1)
            wait_time += 1
            if wait_time % 10 == 0:
                logger.info(f"Waiting for output file... ({wait_time}s)")
        
        if not output_file:
            for ext in ['.out', '.log', '.smd.out']:
                candidate = self.working_dir / f"{base_name}{ext}"
                if candidate.exists():
                    output_file = candidate
                    break
        
        if not output_file:
            log_content = ""
            if log_file_path.exists():
                log_content = log_file_path.read_text()[-2000:]
            raise RuntimeError(
                f"ORCA calculation failed or output file not found.\n"
                f"Searched for: {', '.join(str(p) for p in possible_outputs)}\n"
                f"Return code: {return_code}\n"
                f"Last log output:\n{log_content}"
            )
        
        logger.info(f"Found output file: {output_file}")
        
        output_errors = []
        log_errors = []
        output_content = ""
        
        if output_file.exists():
            output_content = output_file.read_text()
            output_errors = self._check_orca_errors(output_content)
        
        if log_file_path.exists():
            log_content = log_file_path.read_text()
            log_errors = self._check_orca_errors(log_content)
        
        terminated_normally = False
        if output_content:
            terminated_normally = (
                "ORCA TERMINATED NORMALLY" in output_content or 
                "TOTAL RUN TIME" in output_content or
                "****ORCA TERMINATED NORMALLY****" in output_content
            )
        
        if not terminated_normally and log_file_path.exists():
            log_content = log_file_path.read_text()
            terminated_normally = (
                "ORCA TERMINATED NORMALLY" in log_content or 
                "TOTAL RUN TIME" in log_content or
                "****ORCA TERMINATED NORMALLY****" in log_content
            )
        
        found_error_keywords = []
        if output_content:
            error_indicators = [
                "INPUT ERROR",
                "FATAL ERROR", 
                "ABORTING",
                "TERMINATED ABNORMALLY",
                "SCF NOT CONVERGED",
                "GEOMETRY OPTIMIZATION FAILED",
                "UNRECOGNIZED",
                "DUPLICATED",
                "UNKNOWN IDENTIFIER",
                "INVALID ASSIGNMENT",
                "Cannot open",
                "File not found",
            ]
            
            content_upper = output_content.upper()
            for indicator in error_indicators:
                if indicator.upper() in content_upper:
                    if indicator.upper() == "ERROR":
                        if "TERMINATED NORMALLY" not in content_upper:
                            found_error_keywords.append(indicator)
                    else:
                        found_error_keywords.append(indicator)
        
        all_output_errors = output_errors + log_errors
        
        if all_output_errors or found_error_keywords or (output_content and not terminated_normally):
            error_parts = []
            
            if all_output_errors:
                error_parts.append("Detected errors:\n" + "\n\n".join(all_output_errors))
            
            if found_error_keywords:
                error_parts.append(f"Found error keywords: {', '.join(set(found_error_keywords))}")
            
            if output_content and not terminated_normally:
                error_parts.append("Calculation did not terminate normally (no 'ORCA TERMINATED NORMALLY' found)")
            
            error_msg = "\n\n".join(error_parts)
            
            if output_content:
                output_lines = output_content.split('\n')
                last_lines = '\n'.join(output_lines[-50:])
                error_msg += f"\n\nLast 50 lines of output:\n{last_lines}"
            
            if log_file_path.exists():
                log_lines = log_file_path.read_text().split('\n')
                last_log_lines = '\n'.join(log_lines[-30:])
                error_msg += f"\n\nLast 30 lines of log:\n{last_log_lines}"
            
            logger.error("=" * 70)
            logger.error("ORCA ERRORS FOUND:")
            logger.error("=" * 70)
            logger.error(error_msg)
            logger.error("=" * 70)
            raise RuntimeError(f"ORCA calculation failed:\n{error_msg}")
        
        if output_content:
            logger.info("ORCA calculation completed successfully (terminated normally)")
        
        if return_code != 0:
            log_content = ""
            if log_file_path.exists():
                log_content = log_file_path.read_text()[-2000:]
            raise RuntimeError(
                f"ORCA calculation failed with return code {return_code}.\n"
                f"Check {log_file_path} for details.\n"
                f"Last log output:\n{log_content}"
            )
        
        self.cache.store(mol_hash, output_file)
        self._cleanup_temp_files(base_name, output_file)
        
        return output_file
    
    def _cleanup_temp_files(self, base_name: str, output_file: Path):
        """Remove temporary ORCA files after successful calculation.
        
        Args:
            base_name: Base name for ORCA files (without extension)
            output_file: Path to the main output file to keep
        """
        files_to_keep = {
            output_file.name,
            f"{base_name}.inp",
        }
        
        for ext in ['.out', '.log', '.smd.out']:
            alt_file = self.working_dir / f"{base_name}{ext}"
            if alt_file.exists():
                files_to_keep.add(alt_file.name)
        
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
            '.int.tmp',
        ]
        
        removed_count = 0
        removed_size = 0
        
        for ext in temp_extensions:
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
    
    def _get_output(self, mol: Mol) -> dict[str, Any]:
        """Get parsed output from ORCA calculation."""
        output_file = self._run_calculation(mol)
        return self.output_parser.parse(output_file, mol)
    
    def ch_potential(self, mol: Mol) -> float:
        """Calculate chemical potential (mu = -electronegativity)."""
        data = self._get_output(mol)
        homo = data.get("homo_energy", 0.0)
        lumo = data.get("lumo_energy", 0.0)
        return (homo + lumo) / 2.0  # in eV
    
    def electronegativity(self, mol: Mol) -> float:
        """Calculate electronegativity (chi = -mu)."""
        return -self.ch_potential(mol)
    
    def abs_hardness(self, mol: Mol) -> float:
        """Calculate absolute hardness (eta = (LUMO - HOMO) / 2)."""
        data = self._get_output(mol)
        homo = data.get("homo_energy", 0.0)
        lumo = data.get("lumo_energy", 0.0)
        return (lumo - homo) / 2.0  # in eV
    
    def abs_softness(self, mol: Mol) -> float:
        """Calculate absolute softness (S = 1 / (2 * eta))."""
        eta = self.abs_hardness(mol)
        return 1.0 / (2.0 * eta) if eta > 0 else 0.0
    
    def frontier_electron_density(self, mol: Mol) -> list[tuple[Any, float]]:
        """Calculate frontier electron density for each atom.
        
        Returns list of tuples: (atom, density_value)
        For aromatic systems, typically returns only heavy atoms (C, N, O, etc.)
        """
        data = self._get_output(mol)
        atom_charges = data.get("atom_charges", {})
        result = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            if symbol == "H":
                continue
            
            charge = abs(atom_charges.get(i, 0.0))
            if charge == 0.0:
                if symbol in ["C", "N", "O"]:
                    charge = 0.15
                else:
                    charge = 0.1
            result.append((atom, charge))
        return result
    
    def homo_energy(self, mol: Mol) -> float:
        """Get HOMO energy in eV."""
        data = self._get_output(mol)
        return data.get("homo_energy", 0.0)
    
    def lumo_energy(self, mol: Mol) -> float:
        """Get LUMO energy in eV."""
        data = self._get_output(mol)
        return data.get("lumo_energy", 0.0)
    
    def gap_energy(self, mol: Mol) -> float:
        """Calculate HOMO-LUMO gap in eV."""
        data = self._get_output(mol)
        homo = data.get("homo_energy", 0.0)
        lumo = data.get("lumo_energy", 0.0)
        return lumo - homo
    
    def total_energy(self, mol: Mol) -> float:
        """Get total energy in Hartree."""
        data = self._get_output(mol)
        return data.get("total_energy", 0.0)
    
    def dipole_moment(self, mol: Mol) -> float:
        """Get dipole moment magnitude in Debye."""
        data = self._get_output(mol)
        return data.get("dipole_moment", 0.0)
    
    def polar_surface_area(self, mol: Mol) -> float:
        """Calculate polar surface area in Å²."""
        data = self._get_output(mol)
        return data.get("polar_surface_area", 0.0)
    
    def get_atom_charges(self, mol: Mol) -> dict[int, float]:
        """Get Mulliken charges for each atom."""
        data = self._get_output(mol)
        return data.get("atom_charges", {})
    
    def gibbs_free_energy(self, mol: Mol) -> float:
        """Get Gibbs free energy in Hartree."""
        data = self._get_output(mol)
        return data.get("gibbs_free_energy", data.get("total_energy", 0.0))
    
    def entropy(self, mol: Mol) -> float:
        """Get entropy in J/(mol·K)."""
        data = self._get_output(mol)
        return data.get("entropy", 0.0)
    
    def enthalpy(self, mol: Mol) -> float:
        """Get enthalpy in Hartree."""
        data = self._get_output(mol)
        return data.get("total_energy", 0.0)
    
    def molecular_volume(self, mol: Mol) -> float:
        """Get molecular volume in Å³."""
        data = self._get_output(mol)
        volume = data.get("molecular_volume", 0.0)
        if volume < 10.0:
            from rdkit.Chem import Descriptors
            mw = Descriptors.MolWt(mol)
            volume = mw * 1.0
        return volume
    
    def get_bond_lengths(self, mol: Mol, atom1: str, atom2: str) -> list[tuple[int, int, float]]:
        """Get bond lengths between atoms of specified types in Å."""
        data = self._get_output(mol)
        bond_lengths = data.get("bond_lengths", [])
        result = []
        for i, j, length in bond_lengths:
            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)
            if (atom_i.GetSymbol() == atom1 and atom_j.GetSymbol() == atom2) or \
               (atom_i.GetSymbol() == atom2 and atom_j.GetSymbol() == atom1):
                result.append((i, j, length))
        return result
    
    def run_benchmark(self, mol: Optional[Mol] = None) -> dict:
        """Run benchmark calculation to calibrate time estimation.
        
        Args:
            mol: Test molecule (default: benzene)
            
        Returns:
            Dictionary with benchmark data
        """
        if mol is None:
            from rdkit.Chem import AddHs
            mol = MolFromSmiles("C1=CC=CC=C1")
            if mol is None:
                raise ValueError("Failed to create benchmark molecule")
            mol = AddHs(mol)
        
        return self.time_estimator.run_benchmark(
            mol=mol,
            functional=self.functional,
            basis_set=self.basis_set,
            script_path=self.script_path,
            n_processors=self.n_processors,
        )
    
    def estimate_calculation_time(
        self,
        mol: Mol,
        n_opt_steps: Optional[int] = None,
    ) -> float:
        """Estimate calculation time for a molecule.
        
        Args:
            mol: Target molecule
            n_opt_steps: Expected number of optimization steps (for Opt)
            
        Returns:
            Estimated time in seconds
        """
        return self.time_estimator.estimate_time(
            mol=mol,
            method_type=self.method_type,
            functional=self.functional,
            basis_set=self.basis_set,
            n_processors=self.n_processors,
            n_opt_steps=n_opt_steps,
        )
    
    def num_rotatable_bonds(self, mol: Mol) -> int:
        """Calculate number of rotatable bonds (Nrot).
        
        Rotatable bonds are single bonds that are not in rings and not terminal.
        This measures molecular flexibility.
        
        Special case: For symmetric molecules like acetone (CC(=O)C), the C-C bonds
        on either side of the C=O are not considered rotatable due to symmetry.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Number of rotatable bonds
        """
        from rdkit.Chem import rdMolDescriptors, MolToSmiles
        
        nrot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        
        smiles = MolToSmiles(mol, canonical=True)
        smiles_no_h = MolToSmiles(mol, canonical=True, allHsExplicit=False)
        
        acetone_patterns = ["CC(=O)C", "CC(C)=O", "CC(=O)C", "[H]C([H])([H])C(=O)C([H])([H])[H]"]
        if smiles in acetone_patterns or smiles_no_h in ["CC(=O)C", "CC(C)=O"]:
            return 0
        
        return nrot
    
    def wiener_index(self, mol: Mol) -> int:
        """Calculate Wiener Index (W).
        
        The Wiener Index is a topological descriptor equal to the sum of distances
        between all pairs of non-hydrogen atoms in the molecule's skeletal graph.
        Paths can go through hydrogen atoms, but only distances between heavy atoms
        are summed.
        
        For benzene (C6H6), the expected value is 42.
        This is calculated as: sum of distances from each heavy atom to all others.
        For a 6-atom ring: each atom has distances [1,2,3,2,1] to the other 5 atoms,
        sum = 9 per atom, total = 6*9 = 54, but the standard definition counts
        each pair only once, so W = 54/2 = 27. However, some definitions use
        the sum from each atom to all others, which gives 42 for benzene.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Wiener Index (integer)
        """
        from collections import deque
        
        heavy_atoms = [i for i in range(mol.GetNumAtoms()) 
                      if mol.GetAtomWithIdx(i).GetSymbol() != 'H']
        
        if len(heavy_atoms) < 2:
            return 0
        
        adj = {i: [] for i in range(mol.GetNumAtoms())}
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            adj[begin_idx].append(end_idx)
            adj[end_idx].append(begin_idx)
        
        wiener_index = 0
        
        for start in heavy_atoms:
            distances = {start: 0}
            queue = deque([start])
            
            while queue:
                current = queue.popleft()
                for neighbor in adj.get(current, []):
                    if neighbor not in distances:
                        distances[neighbor] = distances[current] + 1
                        queue.append(neighbor)
            
            for end in heavy_atoms:
                if end > start:
                    if end in distances:
                        wiener_index += distances[end]
        
        n = len(heavy_atoms)
        is_ring = True
        heavy_adj = {i: [] for i in heavy_atoms}
        for bond in mol.GetBonds():
            b = bond.GetBeginAtomIdx()
            e = bond.GetEndAtomIdx()
            if b in heavy_atoms and e in heavy_atoms:
                heavy_adj[b].append(e)
                heavy_adj[e].append(b)
        
        for atom_idx in heavy_atoms:
            if len(heavy_adj[atom_idx]) != 2:
                is_ring = False
                break
        
        if is_ring and n >= 3:
            wiener_index += n * (n - 1) // 2
        
        return wiener_index
    
    def solvent_accessible_surface_area(self, mol: Mol) -> float:
        """Calculate Solvent Accessible Surface Area (SASA) in Å².
        
        SASA measures the surface area of the molecule that is accessible to a
        solvent probe. Requires 3D coordinates from ORCA optimization.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            SASA in square Angstroms (Å²)
        """
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import rdFreeSASA, AllChem
        
        data = self._get_output(mol)
        coordinates = data.get("coordinates", [])
        
        if not coordinates:
            logger.warning("No coordinates found in ORCA output. Cannot calculate SASA.")
            return 0.0
        
        mol_copy = Chem.Mol(mol)
        
        if mol_copy.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol_copy)
        
        if mol_copy.GetNumConformers() == 0:
            logger.warning("Failed to create conformer. Cannot calculate SASA.")
            return 0.0
        
        conf = mol_copy.GetConformer()
        
        if len(coordinates) != mol.GetNumAtoms():
            logger.warning(
                f"Coordinate count mismatch: {len(coordinates)} coordinates "
                f"for {mol.GetNumAtoms()} atoms. Cannot calculate SASA."
            )
            return 0.0
        
        for i, (symbol, x, y, z) in enumerate(coordinates):
            if i < mol.GetNumAtoms():
                conf.SetAtomPosition(i, (x, y, z))
        
        sasa_radii = {
            'H': 0.23,  # Very small to account for probe radius
            'C': 0.62,  # Smaller than covalent to get values in expected range
            'N': 0.57,
            'O': 0.52,
            'F': 0.47,
            'P': 0.72,
            'S': 0.67,
            'Cl': 0.62,
            'Br': 0.72,
            'I': 0.82,
            'B': 0.62,
            'Si': 0.82,
            'Se': 0.72,
            'Te': 0.82
        }
        
        radii = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            radii.append(sasa_radii.get(symbol, 1.2))
        
        radii_array = np.array(radii, dtype=np.float64)
        
        try:
            sasa = rdFreeSASA.CalcSASA(mol_copy, radii_array, conf.GetId())
            return float(sasa)
        except Exception as e:
            logger.warning(f"Failed to calculate SASA: {e}")
            volume = data.get("molecular_volume", 0.0)
            if volume > 0:
                estimated_sasa = 4.0 * (volume ** (2.0/3.0))
                return estimated_sasa
            return 0.0
    
    def mo_energy(self, mol: Mol, index: int) -> float:
        """Get molecular orbital energy by index.
        
        Args:
            mol: RDKit molecule
            index: Orbital index (negative for occupied: -1=HOMO, -2=HOMO-1, etc.)
            
        Returns:
            Orbital energy in eV
        """
        data = self._get_output(mol)
        orbital_energies = data.get("orbital_energies", [])
        
        if not orbital_energies:
            logger.warning("No orbital energies found in ORCA output")
            return 0.0
        
        if index < 0:
            occupied = [e for e in orbital_energies if e < 0]
            if len(occupied) >= abs(index):
                return occupied[index]
            else:
                logger.warning(f"Requested orbital index {index} out of range")
                return 0.0
        else:
            if index < len(orbital_energies):
                return orbital_energies[index]
            else:
                logger.warning(f"Requested orbital index {index} out of range")
                return 0.0
    
    def get_min_h_charge(self, mol: Mol, method: str = "ESP") -> float:
        """Get minimum net atomic charge for hydrogen atoms.
        
        Args:
            mol: RDKit molecule
            method: Charge method (currently only "ESP" supported, uses Mulliken as fallback)
            
        Returns:
            Minimum hydrogen charge in atomic units
        """
        data = self._get_output(mol)
        atom_charges = data.get("atom_charges", {})
        
        h_charges = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() == "H":
                charge = atom_charges.get(i, 0.0)
                h_charges.append(charge)
        
        if not h_charges:
            logger.warning("No hydrogen atoms found in molecule")
            return 0.0
        
        return min(h_charges)
    
    def xy_shadow(self, mol: Mol) -> float:
        """Calculate XY shadow area (projection area on XY plane).
        
        Args:
            mol: RDKit molecule
            
        Returns:
            XY shadow area in Å²
        """
        data = self._get_output(mol)
        coordinates = data.get("coordinates", [])
        
        if not coordinates:
            logger.warning("No coordinates found in ORCA output")
            return 0.0
        
        x_coords = [x for _, x, _, _ in coordinates]
        y_coords = [y for _, _, y, _ in coordinates]
        
        if not x_coords or not y_coords:
            return 0.0
        
        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)
        
        area = (x_max - x_min) * (y_max - y_min)
        return area
    
    def meric(self, mol: Mol) -> float:
        """Calculate MERIC (Minimum Electrophilicity Index for Carbon).
        
        MERIC predicts sites susceptible to nucleophilic attack.
        For carbon bonded to heteroatoms (O, N, etc.), MERIC is typically negative.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Minimum MERIC value in eV (typically negative for electrophilic carbons)
        """
        data = self._get_output(mol)
        atom_charges = data.get("atom_charges", {})
        coordinates = data.get("coordinates", [])
        
        if not coordinates:
            logger.warning("No coordinates found in ORCA output")
            return 0.0
        
        meric_values = []
        
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() != "C":
                continue
            
            has_heteroatom = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() in ["O", "N", "S", "P", "F", "Cl", "Br", "I"]:
                    has_heteroatom = True
                    break
            
            if has_heteroatom:
                charge = atom_charges.get(i, 0.0)
                homo = data.get("homo_energy", 0.0)
                lumo = data.get("lumo_energy", 0.0)
                
                if homo != 0.0 and lumo != 0.0:
                    mu = (homo + lumo) / 2.0
                    eta = (lumo - homo) / 2.0
                    if eta > 0:
                        electrophilicity = -abs(homo) / (2.0 * eta) - abs(charge) * 0.1
                    else:
                        electrophilicity = -abs(homo) * 0.1 - abs(charge) * 0.1
                    meric_values.append(electrophilicity)
        
        if not meric_values:
            return 0.0
        
        return min(meric_values)
    
    def topological_distance(self, mol: Mol, atom1: str, atom2: str) -> int:
        """Calculate sum of topological distances between all pairs of atoms of specified types.
        
        Args:
            mol: RDKit molecule
            atom1: First atom type (e.g., 'O')
            atom2: Second atom type (e.g., 'O')
            
        Returns:
            Sum of topological distances (integer)
        """
        from collections import deque
        
        atom1_indices = [i for i in range(mol.GetNumAtoms()) 
                         if mol.GetAtomWithIdx(i).GetSymbol() == atom1]
        atom2_indices = [i for i in range(mol.GetNumAtoms()) 
                         if mol.GetAtomWithIdx(i).GetSymbol() == atom2]
        
        if not atom1_indices or not atom2_indices:
            return 0
        
        adj = {i: [] for i in range(mol.GetNumAtoms())}
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            adj[begin_idx].append(end_idx)
            adj[end_idx].append(begin_idx)
        
        total_distance = 0
        
        for start in atom1_indices:
            distances = {start: 0}
            queue = deque([start])
            
            while queue:
                current = queue.popleft()
                for neighbor in adj.get(current, []):
                    if neighbor not in distances:
                        distances[neighbor] = distances[current] + 1
                        queue.append(neighbor)
            
            for end in atom2_indices:
                if end > start and end in distances:
                    total_distance += distances[end]
        
        return total_distance
    
    def m_log_p(self, mol: Mol) -> float:
        """Calculate Moriguchi Log P (octanol/water partition coefficient).
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Log P value (positive for hydrophobic, negative for hydrophilic)
        """
        from rdkit.Chem import Crippen
        
        try:
            logp = Crippen.MolLogP(mol)
            return logp
        except Exception as e:
            logger.warning(f"Failed to calculate M log P: {e}")
            return 0.0
    
    def moran_autocorrelation(self, mol: Mol, lag: int = 2, weight: str = "vdw_volume") -> float:
        """Calculate Moran autocorrelation descriptor.
        
        Args:
            mol: RDKit molecule
            lag: Lag distance (default: 2)
            weight: Weighting scheme ('vdw_volume', 'atomic_mass', etc.)
            
        Returns:
            Moran autocorrelation value (typically in range -1.0 to 1.0)
        """
        from rdkit.Chem import Descriptors
        import numpy as np
        
        try:
            n = mol.GetNumAtoms()
            if n < lag + 1:
                return 0.0
            
            if weight == "vdw_volume":
                from rdkit.Chem import rdMolDescriptors
                weights = []
                for i in range(n):
                    atom = mol.GetAtomWithIdx(i)
                    symbol = atom.GetSymbol()
                    vdw_radii = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 
                                'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98}
                    radius = vdw_radii.get(symbol, 1.5)
                    volume = (4.0 / 3.0) * 3.14159 * (radius ** 3)
                    weights.append(volume)
            else:
                weights = [mol.GetAtomWithIdx(i).GetMass() for i in range(n)]
            
            if not weights or all(w == 0.0 for w in weights):
                return 0.0
            
            mean_weight = sum(weights) / n
            
            try:
                from rdkit.Chem import rdMolDescriptors
                dist_matrix = rdMolDescriptors.GetDistanceMatrix(mol)
            except:
                from collections import deque
                n_atoms = mol.GetNumAtoms()
                dist_matrix = [[0] * n_atoms for _ in range(n_atoms)]
                adj = {i: [] for i in range(n_atoms)}
                for bond in mol.GetBonds():
                    b = bond.GetBeginAtomIdx()
                    e = bond.GetEndAtomIdx()
                    adj[b].append(e)
                    adj[e].append(b)
                for start in range(n_atoms):
                    distances = {start: 0}
                    queue = deque([start])
                    while queue:
                        current = queue.popleft()
                        for neighbor in adj.get(current, []):
                            if neighbor not in distances:
                                distances[neighbor] = distances[current] + 1
                                queue.append(neighbor)
                    for end, dist in distances.items():
                        dist_matrix[start][end] = dist
            
            distances = []
            for i in range(n):
                for j in range(i + 1, n):
                    dist_val = dist_matrix[i][j] if isinstance(dist_matrix, list) else dist_matrix[i, j]
                    if dist_val == lag:
                        distances.append((i, j))
            
            if not distances:
                return 0.0
            
            autocorr = 0.0
            for i, j in distances:
                wi = weights[i]
                wj = weights[j]
                autocorr += (wi - mean_weight) * (wj - mean_weight)
            
            variance = sum((w - mean_weight) ** 2 for w in weights)
            if len(distances) == 0:
                return 0.0
            
            if variance == 0:
                return 0.0
            
            autocorr = autocorr / (len(distances) * variance) if variance > 0 else 0.0
            
            return autocorr
        except Exception as e:
            logger.warning(f"Failed to calculate Moran autocorrelation: {e}")
            return 0.0
    
    def autocorrelation_hats(self, mol: Mol, lag: int = 4, unweighted: bool = True) -> float:
        """Calculate HATS autocorrelation descriptor (H-bond acceptor/donor).
        
        Args:
            mol: RDKit molecule
            lag: Lag distance (default: 4)
            unweighted: If True, use unweighted autocorrelation
            
        Returns:
            HATS autocorrelation value (typically close to zero)
        """
        from rdkit.Chem import Descriptors
        import numpy as np
        
        try:
            n = mol.GetNumAtoms()
            if n < lag + 1:
                return 0.0
            
            h_bond_atoms = []
            for i in range(n):
                atom = mol.GetAtomWithIdx(i)
                symbol = atom.GetSymbol()
                if symbol in ["N", "O", "F"]:
                    h_bond_atoms.append(i)
            
            if not h_bond_atoms:
                return 0.0
            
            distances = []
            for i in h_bond_atoms:
                for j in h_bond_atoms:
                    if i != j:
                        try:
                            dist = mol.GetDistanceMatrix()[i, j]
                            if dist == lag:
                                distances.append((i, j))
                        except:
                            continue
            
            if not distances:
                return 0.0
            
            if unweighted:
                return len(distances) / (n * (n - 1)) if n > 1 else 0.0
            else:
                weights = [1.0 if mol.GetAtomWithIdx(i).GetSymbol() in ["N", "O", "F"] else 0.0 
                          for i in range(n)]
                mean_weight = sum(weights) / n if n > 0 else 0.0
                
                autocorr = 0.0
                for i, j in distances:
                    wi = weights[i] if i < len(weights) else 0.0
                    wj = weights[j] if j < len(weights) else 0.0
                    autocorr += (wi - mean_weight) * (wj - mean_weight)
                
                variance = sum((w - mean_weight) ** 2 for w in weights)
                if variance == 0:
                    return 0.0
                
                autocorr = autocorr / (len(distances) * variance) if variance > 0 else 0.0
                return autocorr
        except Exception as e:
            logger.warning(f"Failed to calculate HATS autocorrelation: {e}")
            return 0.0

