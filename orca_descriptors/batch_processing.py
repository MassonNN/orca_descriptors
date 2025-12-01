"""Batch processing for ORCA calculations with pandas compatibility."""

import hashlib
import logging
import multiprocessing
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Optional, Union

from rdkit.Chem import Mol, MolFromSmiles, AddHs

from orca_descriptors.orca import Orca

logger = logging.getLogger(__name__)


class DescriptorCall:
    """Represents a descriptor method call with its parameters.
    
    This class is used to capture descriptor method calls with their parameters
    when using the x_molecule() API.
    """
    def __init__(self, method_name: str, args: tuple = (), kwargs: dict = None):
        """Initialize descriptor call.
        
        Args:
            method_name: Name of the descriptor method
            args: Positional arguments (excluding the molecule)
            kwargs: Keyword arguments
        """
        self.method_name = method_name
        self.args = args
        self.kwargs = kwargs or {}
    
    def __repr__(self):
        args_str = ", ".join(repr(a) for a in self.args)
        kwargs_str = ", ".join(f"{k}={v!r}" for k, v in self.kwargs.items())
        params = ", ".join(filter(None, [args_str, kwargs_str]))
        return f"DescriptorCall({self.method_name!r}, {params})"


class XMolecule:
    """Mock molecule for defining descriptors with parameters.
    
    This class is used as a placeholder molecule when calling descriptor methods
    on an Orca instance. When a descriptor method is called with an XMolecule
    as the first argument, it returns a DescriptorCall object instead of
    performing the actual calculation.
    
    Example::
    
        x = batch_processing.x_molecule()
        desc = orca.ch_potential(x)  # Returns DescriptorCall('ch_potential', (), {})
        desc = orca.mo_energy(x, -3)  # Returns DescriptorCall('mo_energy', (-3,), {})
    """
    def __init__(self, orca: Optional[Orca] = None):
        """Initialize X molecule.
        
        Args:
            orca: Optional Orca instance (not currently used, but kept for API compatibility)
        """
        self.orca = orca


def _calculate_worker_multiprocessing(args: tuple) -> tuple[int, dict[str, Any]]:
    """Worker function for multiprocessing (top-level function for pickling).
    
    Args:
        args: Tuple of (idx, smiles, descriptors, total, orca_params)
              where orca_params is a dict with all Orca initialization parameters
              and descriptors can be list of strings or list of DescriptorCall objects
        
    Returns:
        Tuple of (idx, result_dict)
    """
    from orca_descriptors.batch_processing import DescriptorCall
    
    idx, smiles, descriptors, total, orca_params = args
    orca = Orca(**orca_params)
    try:
        mol = MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"[{idx + 1}/{total}] Failed to parse SMILES: {smiles}")
            if descriptors and isinstance(descriptors[0], DescriptorCall):
                keys = []
                for d in descriptors:
                    param_str = ""
                    if d.args:
                        param_str = "_" + "_".join(str(a) for a in d.args)
                    if d.kwargs:
                        param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(d.kwargs.items()))
                    keys.append(f"{d.method_name}{param_str}")
            else:
                keys = descriptors if descriptors else []
            return idx, {desc: None for desc in keys}
        mol = AddHs(mol)
    except Exception as e:
        logger.debug(f"[{idx + 1}/{total}] Failed to parse SMILES '{smiles}': {e}")
        if descriptors and isinstance(descriptors[0], DescriptorCall):
            keys = []
            for d in descriptors:
                param_str = ""
                if d.args:
                    param_str = "_" + "_".join(str(a) for a in d.args)
                if d.kwargs:
                    param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(d.kwargs.items()))
                keys.append(f"{d.method_name}{param_str}")
        else:
            keys = descriptors if descriptors else []
        return idx, {desc: None for desc in keys}
    
    result_descriptors = {}
    special_descriptors = {
        'moran_autocorrelation': {'lag': 2, 'weight': 'vdw_volume'},
        'autocorrelation_hats': {'lag': 4, 'unweighted': True},
    }
    
    for desc in descriptors:
        try:
            if isinstance(desc, DescriptorCall):
                method_name = desc.method_name
                method = getattr(orca, method_name, None)
                
                if method is None or not callable(method):
                    logger.debug(f"[{idx + 1}/{total}] Method '{method_name}' not found or not callable")
                    param_str = ""
                    if desc.args:
                        param_str = "_" + "_".join(str(a) for a in desc.args)
                    if desc.kwargs:
                        param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(desc.kwargs.items()))
                    result_key = f"{method_name}{param_str}"
                    result_descriptors[result_key] = None
                    continue
                
                result = method(mol, *desc.args, **desc.kwargs)
                
                param_str = ""
                if desc.args:
                    param_str = "_" + "_".join(str(a) for a in desc.args)
                if desc.kwargs:
                    param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(desc.kwargs.items()))
                result_key = f"{method_name}{param_str}"
                
                if method_name == 'frontier_electron_density' and isinstance(result, list):
                    result = max(charge for _, charge in result) if result else 0.0
                
                result_descriptors[result_key] = result
            else:
                desc_name = desc
                method = getattr(orca, desc_name, None)
                
                if method is None or not callable(method):
                    logger.debug(f"[{idx + 1}/{total}] Method '{desc_name}' not found or not callable")
                    result_descriptors[desc_name] = None
                    continue
                
                if desc_name in special_descriptors:
                    result = method(mol, **special_descriptors[desc_name])
                elif desc_name == 'frontier_electron_density':
                    frontier_density = method(mol)
                    result = max(charge for _, charge in frontier_density) if frontier_density else 0.0
                else:
                    result = method(mol)
                
                result_descriptors[desc_name] = result
        except Exception as e:
            error_msg = str(e)
            brief_error = error_msg.split('\n')[0] if error_msg else str(e)
            
            if isinstance(desc, DescriptorCall):
                desc_name = desc.method_name
                param_str = ""
                if desc.args:
                    param_str = "_" + "_".join(str(a) for a in desc.args)
                if desc.kwargs:
                    param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(desc.kwargs.items()))
                result_key = f"{desc_name}{param_str}"
            else:
                desc_name = desc
                result_key = desc_name
            
            logger.info(f"[{idx + 1}/{total}] Error calculating '{desc_name}' for SMILES {smiles}: {brief_error}")
            logger.debug(f"[{idx + 1}/{total}] Full error: {error_msg}")
            result_descriptors[result_key] = None
    
    try:
        mol_hash = orca._get_molecule_hash(mol)
        base_name = f"orca_{mol_hash}"
        possible_outputs = [
            orca.working_dir / f"{base_name}.out",
            orca.working_dir / f"{base_name}.log",
            orca.working_dir / f"{base_name}.smd.out",
        ]
        for output_file in possible_outputs:
            if output_file.exists():
                orca._cleanup_temp_files(base_name, output_file)
                break
    except Exception as cleanup_error:
        logger.debug(f"Failed to cleanup files for {smiles}: {cleanup_error}")
    
    return idx, result_descriptors


class ORCABatchProcessing:
    """Batch processing class for ORCA calculations with pandas compatibility.
    
    This class provides efficient batch processing of molecular descriptors
    with support for parallelization via mpirun (ORCA's built-in parallelization)
    or Python multiprocessing.
    """
    
    def __init__(
        self,
        orca: Optional[Orca] = None,
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
        use_mpirun: bool = False,
        mpirun_path: Optional[str] = None,
        extra_env: Optional[dict] = None,
        parallel_mode: str = "sequential",
        n_workers: Optional[int] = None,
        pre_optimize: bool = True,
    ):
        """Initialize ORCA batch processor.
        
        Args:
            orca: Existing Orca instance to use (if provided, other parameters are ignored)
            script_path: Path to ORCA executable
            working_dir: Working directory for calculations
            output_dir: Directory for output files
            functional: DFT functional (e.g., "PBE0")
            basis_set: Basis set (e.g., "def2-SVP")
            method_type: Calculation type ("Opt", "SP", etc.)
            dispersion_correction: Dispersion correction (e.g., "D3BJ")
            solvation_model: Solvation model (e.g., "COSMO(Water)")
            n_processors: Number of processors per ORCA calculation
            max_scf_cycles: Maximum SCF cycles
            scf_convergence: SCF convergence threshold
            charge: Molecular charge
            multiplicity: Spin multiplicity
            cache_dir: Directory for caching results (default: output_dir/.orca_cache)
            log_level: Logging level (default: logging.INFO)
            max_wait: Maximum time to wait for output file creation in seconds (default: 300)
            use_mpirun: Whether to use mpirun for parallel execution within each ORCA calculation
            mpirun_path: Path to mpirun executable (default: None, will search in PATH)
            extra_env: Additional environment variables to pass to ORCA process (default: None)
            parallel_mode: Parallelization mode: "sequential", "multiprocessing", or "mpirun"
                          "sequential": Process molecules one by one
                          "multiprocessing": Use Python multiprocessing to run multiple ORCA calculations in parallel
                          "mpirun": Use mpirun for each ORCA calculation (ORCA's built-in parallelization)
            n_workers: Number of parallel workers for multiprocessing mode (default: number of CPUs)
            pre_optimize: Whether to pre-optimize geometry with MMFF94 before ORCA calculation (default: True)
        """
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            handler.setLevel(log_level)
            logger.addHandler(handler)
        else:
            for handler in logger.handlers:
                handler.setLevel(log_level)
        logger.setLevel(log_level)
        logger.propagate = False
        
        if orca is not None:
            self.orca = orca
        else:
            self.orca = Orca(
                script_path=script_path,
                working_dir=working_dir,
                output_dir=output_dir,
                functional=functional,
                basis_set=basis_set,
                method_type=method_type,
                dispersion_correction=dispersion_correction,
                solvation_model=solvation_model,
                n_processors=n_processors,
                max_scf_cycles=max_scf_cycles,
                scf_convergence=scf_convergence,
                charge=charge,
                multiplicity=multiplicity,
                cache_dir=cache_dir,
                log_level=log_level,
                max_wait=max_wait,
                use_mpirun=use_mpirun,
                mpirun_path=mpirun_path,
                extra_env=extra_env,
                pre_optimize=pre_optimize,
            )
        
        self.parallel_mode = parallel_mode
        if n_workers is None:
            self.n_workers = multiprocessing.cpu_count()
        else:
            self.n_workers = n_workers
        
        self.orca.working_dir.mkdir(parents=True, exist_ok=True)
        self.orca.output_dir.mkdir(parents=True, exist_ok=True)
    
    def x_molecule(self) -> XMolecule:
        """Create a mock molecule for defining descriptors with parameters.
        
        This method returns a special XMolecule object that can be used to
        define descriptors with their parameters. When you call a descriptor
        method on the Orca instance with this X molecule, it returns a
        DescriptorCall object that captures the method name and parameters.
        
        Example::
        
            x = batch_processing.x_molecule()
            descriptors = [
                orca.ch_potential(x),
                orca.topological_distance(x, 'O', 'O'),
                orca.mo_energy(x, -3)
            ]
            result = batch_processing.calculate_descriptors(smiles_list, descriptors=descriptors)
        
        Returns:
            XMolecule instance for defining descriptors
        """
        return XMolecule(self.orca)
    
    def _prepare_molecule(self, smiles: str) -> Optional[Mol]:
        """Prepare molecule from SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            RDKit molecule with hydrogens added, or None if parsing failed
        """
        try:
            mol = MolFromSmiles(smiles)
            if mol is None:
                return None
            return AddHs(mol)
        except Exception as e:
            logger.debug(f"Failed to parse SMILES '{smiles}': {e}")
            return None
    
    def _estimate_times(self, smiles_list: list[str]) -> list[float]:
        """Estimate calculation times for all molecules.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            List of estimated times in seconds
        """
        estimated_times = []
        for smiles in smiles_list:
            try:
                mol = self._prepare_molecule(smiles)
                if mol is not None:
                    estimated_time = self.orca.time_estimator.estimate_time(
                        mol=mol,
                        method_type=self.orca.method_type,
                        functional=self.orca.functional,
                        basis_set=self.orca.basis_set,
                        n_processors=self.orca.n_processors,
                    )
                    estimated_times.append(estimated_time)
                else:
                    estimated_times.append(0.0)
            except Exception as e:
                logger.debug(f"Failed to estimate time for SMILES '{smiles}': {e}")
                estimated_times.append(0.0)
        return estimated_times
    
    def _get_parallel_efficiency(self, n_workers: int) -> float:
        """Get parallel efficiency factor based on number of workers.
        
        Args:
            n_workers: Number of parallel workers
            
        Returns:
            Efficiency factor (0.0-1.0), typically 0.7-0.9
        """
        if n_workers <= 4:
            return 0.85
        elif n_workers <= 8:
            return 0.80
        elif n_workers <= 16:
            return 0.75
        else:
            return 0.70
    
    def _format_time(self, seconds: float) -> str:
        """Format time in seconds to human-readable string.
        
        Args:
            seconds: Time in seconds
            
        Returns:
            Formatted time string (e.g., "2h 30m 15s", "45m 30s", "30s")
        """
        if seconds <= 0:
            return "0s"
        
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = int(seconds % 60)
        
        if hours > 0:
            return f"{hours}h {minutes}m {secs}s"
        elif minutes > 0:
            return f"{minutes}m {secs}s"
        else:
            return f"{secs}s"
    
    def _parse_orca_error(self, content: str) -> tuple[str, str]:
        """Parse ORCA error from output content.
        
        Extracts a brief error summary for INFO logging and full error details
        for DEBUG logging.
        
        Args:
            content: ORCA output content
            
        Returns:
            Tuple of (brief_error, detailed_error)
        """
        errors = []
        lines = content.split('\n')
        
        error_patterns = [
            (r'INPUT ERROR', 'Input error'),
            (r'FATAL ERROR', 'Fatal error'),
            (r'ABORTING', 'Calculation aborted'),
            (r'TERMINATED ABNORMALLY', 'Terminated abnormally'),
            (r'SCF NOT CONVERGED', 'SCF not converged'),
            (r'GEOMETRY OPTIMIZATION FAILED', 'Geometry optimization failed'),
            (r'UNRECOGNIZED.*KEYWORD', 'Unrecognized keyword'),
            (r'DUPLICATED.*KEYWORD', 'Duplicated keyword'),
            (r'Unknown identifier', 'Unknown identifier'),
            (r'Invalid assignment', 'Invalid assignment'),
            (r'Cannot open', 'Cannot open file'),
            (r'File not found', 'File not found'),
            (r'ERROR.*finished by error', 'Error termination'),
            (r'ERROR.*aborting', 'Error aborting'),
            (r'ERROR.*termination', 'Error termination'),
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
        
        if errors:
            brief_error = errors[0].split('\n')[0]
            detailed_error = "\n\n".join(errors)
            return brief_error, detailed_error
        
        terminated_normally = (
            "ORCA TERMINATED NORMALLY" in content or 
            "TOTAL RUN TIME" in content or
            "****ORCA TERMINATED NORMALLY****" in content
        )
        
        if not terminated_normally:
            return "Calculation did not terminate normally", content[-500:]
        
        return "", ""
    
    def _cleanup_temp_files(self, base_name: str, output_file: Path):
        """Remove all ORCA files after calculation.
        
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
                for orca_file in self.orca.working_dir.glob(pattern):
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
    
    def _calculate_single_molecule(
        self,
        smiles: str,
        descriptors: Union[list[str], list[DescriptorCall]],
        idx: int,
        total: int,
    ) -> dict[str, Any]:
        """Calculate descriptors for a single molecule.
        
        Args:
            smiles: SMILES string
            descriptors: List of descriptor names (str) or DescriptorCall objects
            idx: Current molecule index (0-based)
            total: Total number of molecules
            
        Returns:
            Dictionary with descriptor values
        """
        start_time = time.time()
        
        try:
            mol = self._prepare_molecule(smiles)
            if mol is None:
                logger.warning(f"[{idx + 1}/{total}] Failed to parse SMILES: {smiles}")
                if descriptors and isinstance(descriptors[0], DescriptorCall):
                    keys = [f"{d.method_name}_{'_'.join(str(a) for a in d.args)}" for d in descriptors]
                else:
                    keys = descriptors if descriptors else []
                return {desc: None for desc in keys}
            
            result_descriptors = {}
            special_descriptors = {
                'moran_autocorrelation': {'lag': 2, 'weight': 'vdw_volume'},
                'autocorrelation_hats': {'lag': 4, 'unweighted': True},
            }
            
            for desc in descriptors:
                try:
                    if isinstance(desc, DescriptorCall):
                        method_name = desc.method_name
                        method = getattr(self.orca, method_name, None)
                        
                        if method is None or not callable(method):
                            logger.debug(f"[{idx + 1}/{total}] Method '{method_name}' not found or not callable")
                            result_key = f"{method_name}_{'_'.join(str(a) for a in desc.args)}"
                            result_descriptors[result_key] = None
                            continue
                        
                        result = method(mol, *desc.args, **desc.kwargs)
                        
                        param_str = ""
                        if desc.args:
                            param_str = "_" + "_".join(str(a) for a in desc.args)
                        if desc.kwargs:
                            param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(desc.kwargs.items()))
                        result_key = f"{method_name}{param_str}"
                        
                        if method_name == 'frontier_electron_density' and isinstance(result, list):
                            result = max(charge for _, charge in result) if result else 0.0
                        
                        result_descriptors[result_key] = result
                    else:
                        desc_name = desc
                        method = getattr(self.orca, desc_name, None)
                        
                        if method is None or not callable(method):
                            logger.debug(f"[{idx + 1}/{total}] Method '{desc_name}' not found or not callable")
                            result_descriptors[desc_name] = None
                            continue
                        
                        if desc_name in special_descriptors:
                            result = method(mol, **special_descriptors[desc_name])
                        elif desc_name == 'frontier_electron_density':
                            frontier_density = method(mol)
                            result = max(charge for _, charge in frontier_density) if frontier_density else 0.0
                        else:
                            result = method(mol)
                        
                        result_descriptors[desc_name] = result
                except Exception as e:
                    error_msg = str(e)
                    brief_error, detailed_error = self._parse_orca_error(error_msg)
                    
                    if isinstance(desc, DescriptorCall):
                        desc_name = desc.method_name
                        param_str = ""
                        if desc.args:
                            param_str = "_" + "_".join(str(a) for a in desc.args)
                        if desc.kwargs:
                            param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(desc.kwargs.items()))
                        result_key = f"{desc_name}{param_str}"
                    else:
                        desc_name = desc
                        result_key = desc_name
                    
                    if brief_error:
                        logger.info(f"[{idx + 1}/{total}] Error calculating '{desc_name}' for SMILES {smiles}: {brief_error}")
                        logger.debug(f"[{idx + 1}/{total}] Detailed error for '{desc_name}': {detailed_error}")
                    else:
                        logger.info(f"[{idx + 1}/{total}] Error calculating '{desc_name}' for SMILES {smiles}: {e}")
                        logger.debug(f"[{idx + 1}/{total}] Full error traceback: {error_msg}")
                    
                    result_descriptors[result_key] = None
            
            try:
                mol_hash = self.orca._get_molecule_hash(mol)
                base_name = f"orca_{mol_hash}"
                possible_outputs = [
                    self.orca.working_dir / f"{base_name}.out",
                    self.orca.working_dir / f"{base_name}.log",
                    self.orca.working_dir / f"{base_name}.smd.out",
                ]
                for output_file in possible_outputs:
                    if output_file.exists():
                        self._cleanup_temp_files(base_name, output_file)
                        break
            except Exception as cleanup_error:
                logger.debug(f"Failed to cleanup files for {smiles}: {cleanup_error}")
            
            actual_time = time.time() - start_time
            logger.debug(f"[{idx + 1}/{total}] Completed in {actual_time:.2f}s")
            
            return result_descriptors
            
        except Exception as e:
            error_msg = str(e)
            brief_error, detailed_error = self._parse_orca_error(error_msg)
            
            if brief_error:
                logger.info(f"[{idx + 1}/{total}] Error processing SMILES {smiles}: {brief_error}")
                logger.debug(f"[{idx + 1}/{total}] Detailed error: {detailed_error}")
            else:
                logger.info(f"[{idx + 1}/{total}] Error processing SMILES {smiles}: {e}")
                logger.debug(f"[{idx + 1}/{total}] Full error traceback: {error_msg}")
            
            # Create keys for result dict
            if descriptors and isinstance(descriptors[0], DescriptorCall):
                keys = []
                for d in descriptors:
                    param_str = ""
                    if d.args:
                        param_str = "_" + "_".join(str(a) for a in d.args)
                    if d.kwargs:
                        param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(d.kwargs.items()))
                    keys.append(f"{d.method_name}{param_str}")
            else:
                keys = descriptors if descriptors else []
            return {desc: None for desc in keys}
    
    
    def _get_available_descriptors(self) -> list[str]:
        """Get list of all available descriptor method names.
        
        Returns:
            List of descriptor method names
        """
        excluded = {
            'run_benchmark', 'estimate_calculation_time', 'calculate_descriptors',
            '_get_available_descriptors', 'get_atom_charges', 'get_bond_lengths',
        }
        
        descriptor_methods = []
        for name in dir(self.orca):
            if not name.startswith('_') and name not in excluded:
                method = getattr(self.orca, name, None)
                if callable(method):
                    descriptor_methods.append(name)
        
        return sorted(descriptor_methods)
    
    def calculate_descriptors(
        self,
        smiles_column: Union[Any, list[str]],
        descriptors: Optional[list[str]] = None,
        progress: bool = True,
    ) -> Any:
        """Calculate descriptors for molecules from SMILES and add to DataFrame.
        
        This method provides pandas compatibility. It accepts a pandas Series or
        DataFrame column with SMILES strings and returns a DataFrame with added
        descriptor columns. Only descriptor columns are added - no new 'smiles'
        column is created.
        
        - If input is a Series: returns DataFrame with only descriptor columns
        - If input is a DataFrame: returns DataFrame with original columns + descriptor columns
        
        By default, calculates all available descriptors. Use the `descriptors`
        parameter to specify a subset of descriptors to calculate.
        
        Args:
            smiles_column: pandas Series/DataFrame column with SMILES strings,
                          or a list of SMILES strings
            descriptors: Optional list of descriptor names to calculate.
                       If None, calculates all available descriptors.
                       Descriptor names correspond to method names of the Orca class.
            progress: Whether to show progress (default: True)
            
        Returns:
            DataFrame with descriptor columns added (if pandas available).
            For Series input: only descriptor columns.
            For DataFrame input: original columns + descriptor columns.
            For list input: list of dictionaries (if pandas not available)
            
        Raises:
            ImportError: If pandas is not installed and a pandas object is passed
            ValueError: If an invalid descriptor name is provided
        """
        try:
            import pandas as pd
            pandas_available = True
        except ImportError:
            pandas_available = False
        
        if pandas_available:
            if isinstance(smiles_column, pd.Series):
                # For Series, we don't create a DataFrame with 'smiles' column
                # We'll only return descriptors
                df = None
                smiles_list = smiles_column.tolist()
            elif isinstance(smiles_column, pd.DataFrame):
                if 'smiles' not in smiles_column.columns:
                    raise ValueError("DataFrame must contain a 'smiles' column")
                df = smiles_column.copy()
                smiles_list = smiles_column['smiles'].tolist()
            elif isinstance(smiles_column, list):
                df = None
                smiles_list = smiles_column
            else:
                raise TypeError(
                    "smiles_column must be a pandas Series, DataFrame, or a list of SMILES strings. "
                    f"Got {type(smiles_column)}"
                )
        else:
            if isinstance(smiles_column, list):
                df = None
                smiles_list = smiles_column
            else:
                raise ImportError(
                    "pandas is required for DataFrame/Series input. "
                    "Install pandas with: pip install 'orca-descriptors[pandas]' "
                    "or pass a list of SMILES strings."
                )
        
        if descriptors is not None:
            descriptor_calls = []
            descriptor_names = []
            
            for desc in descriptors:
                if isinstance(desc, DescriptorCall):
                    descriptor_calls.append(desc)
                    param_str = ""
                    if desc.args:
                        param_str = "_" + "_".join(str(a) for a in desc.args)
                    if desc.kwargs:
                        param_str += "_" + "_".join(f"{k}_{v}" for k, v in sorted(desc.kwargs.items()))
                    descriptor_names.append(f"{desc.method_name}{param_str}")
                elif isinstance(desc, str):
                    descriptor_names.append(desc)
                else:
                    raise ValueError(
                        f"Invalid descriptor type: {type(desc)}. "
                        f"Expected str or DescriptorCall, got {desc}"
                    )
            
            if descriptor_calls:
                descriptors_to_calculate = descriptor_calls
            else:
                available = self._get_available_descriptors()
                invalid_descriptors = [d for d in descriptor_names if d not in available]
                if invalid_descriptors:
                    raise ValueError(
                        f"Invalid descriptor names: {invalid_descriptors}. "
                        f"Available descriptors: {', '.join(available)}"
                    )
                descriptors_to_calculate = descriptor_names
        else:
            descriptors_to_calculate = self._get_available_descriptors()
        
        total = len(smiles_list)
        if total == 0:
            if pandas_available and df is not None:
                return df
            else:
                return []
        
        estimated_times = []
        if progress and total > 1:
            logger.info("Estimating calculation times...")
            estimated_times = self._estimate_times(smiles_list)
            total_estimated = sum(estimated_times)
            
            if self.parallel_mode == "multiprocessing" and total > 1 and self.n_workers > 1:
                efficiency = self._get_parallel_efficiency(self.n_workers)
                parallel_estimated = total_estimated / (self.n_workers * efficiency)
                logger.info(
                    f"Total estimated time (sequential): {self._format_time(total_estimated)}"
                )
                logger.info(
                    f"Total estimated time (parallel, {self.n_workers} workers, "
                    f"efficiency {efficiency:.0%}): {self._format_time(parallel_estimated)}"
                )
                total_estimated = parallel_estimated
            elif total_estimated > 0:
                logger.info(f"Total estimated time: {self._format_time(total_estimated)}")
        
        if self.parallel_mode == "multiprocessing" and total > 1:
            efficiency = self._get_parallel_efficiency(self.n_workers)
            logger.info(
                f"Processing {total} molecules using multiprocessing with {self.n_workers} workers "
                f"(parallel efficiency: {efficiency:.0%})"
            )
            
            # Prepare Orca parameters for worker processes
            orca_params = {
                'script_path': self.orca.script_path,
                'working_dir': str(self.orca.working_dir),
                'output_dir': str(self.orca.output_dir),
                'functional': self.orca.functional,
                'basis_set': self.orca.basis_set,
                'method_type': self.orca.method_type,
                'dispersion_correction': self.orca.dispersion_correction,
                'solvation_model': self.orca.solvation_model,
                'n_processors': self.orca.n_processors,
                'max_scf_cycles': self.orca.max_scf_cycles,
                'scf_convergence': self.orca.scf_convergence,
                'charge': self.orca.charge,
                'multiplicity': self.orca.multiplicity,
                'cache_dir': str(self.orca.cache.cache_dir) if hasattr(self.orca.cache, 'cache_dir') else None,
                'log_level': logging.DEBUG,
                'max_wait': self.orca.max_wait,
                'use_mpirun': self.orca.use_mpirun,
                'mpirun_path': self.orca.mpirun_path,
                'extra_env': self.orca.extra_env,
                'pre_optimize': self.orca.pre_optimize,
            }
            
            with multiprocessing.Pool(processes=self.n_workers) as pool:
                args_list = [
                    (idx, smiles, descriptors_to_calculate, total, orca_params)
                    for idx, smiles in enumerate(smiles_list)
                ]
                results = pool.map(_calculate_worker_multiprocessing, args_list)
            
            results.sort(key=lambda x: x[0])
            descriptors_list = [result[1] for result in results]
        else:
            descriptors_list = []
            actual_times = []
            
            for idx, smiles in enumerate(smiles_list):
                is_cached = False
                try:
                    mol = self._prepare_molecule(smiles)
                    if mol is not None:
                        mol_hash = self.orca._get_molecule_hash(mol)
                        cached_output = self.orca.cache.get(mol_hash)
                        is_cached = cached_output is not None and cached_output.exists()
                except Exception:
                    pass
                
                if progress and total > 1:
                    remaining = total - idx
                    
                    if is_cached:
                        if actual_times:
                            avg_actual = sum(actual_times) / len(actual_times)
                            logger.info(
                                f"Processing molecule {idx + 1}/{total} (remaining: {remaining}, "
                                f"CACHED, avg: {avg_actual:.1f}s/molecule)"
                            )
                        else:
                            logger.info(f"Processing molecule {idx + 1}/{total} (remaining: {remaining}, CACHED)")
                    else:
                        if actual_times and len(actual_times) > 0:
                            alpha = 0.3
                            if len(actual_times) == 1:
                                avg_time = actual_times[0]
                            else:
                                avg_time = actual_times[0]
                                for t in actual_times[1:]:
                                    avg_time = alpha * t + (1 - alpha) * avg_time
                            remaining_estimated = avg_time * remaining
                        else:
                            remaining_estimated = sum(estimated_times[idx:]) if estimated_times else 0.0
                        
                        if remaining_estimated > 0:
                            time_str = f"~{self._format_time(remaining_estimated)}"
                            
                            if actual_times:
                                avg_actual = sum(actual_times) / len(actual_times)
                                logger.info(
                                    f"Processing molecule {idx + 1}/{total} (remaining: {remaining}, "
                                    f"estimated time: {time_str}, avg: {avg_actual:.1f}s/molecule)"
                                )
                            else:
                                logger.info(f"Processing molecule {idx + 1}/{total} (remaining: {remaining}, estimated time: {time_str})")
                        else:
                            logger.info(f"Processing molecule {idx + 1}/{total} (remaining: {remaining})")
                
                start_time = time.time()
                result = self._calculate_single_molecule(smiles, descriptors_to_calculate, idx, total)
                descriptors_list.append(result)
                actual_times.append(time.time() - start_time)
        
        if pandas_available:
            descriptor_df = pd.DataFrame(descriptors_list)
            
            if df is not None:
                result_df = pd.concat([df.reset_index(drop=True), descriptor_df.reset_index(drop=True)], axis=1)
                return result_df
            else:
                return descriptor_df
        else:
            return descriptors_list

