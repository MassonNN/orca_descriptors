"""Base class and utility methods for ORCA calculations."""

import hashlib
import logging
from pathlib import Path
from typing import Any, Optional

from rdkit.Chem import Mol, MolFromSmiles, AddHs

from orca_descriptors.cache import CacheManager
from orca_descriptors.input_generator import ORCAInputGenerator
from orca_descriptors.output_parser import ORCAOutputParser
from orca_descriptors.time_estimator import ORCATimeEstimator

logger = logging.getLogger(__name__)


class OrcaBase:
    """Base class for ORCA calculations with common utility methods."""
    
    def _is_semi_empirical(self) -> bool:
        """Check if the current method is semi-empirical.
        
        Returns:
            True if semi-empirical, False otherwise
        """
        return self.input_generator._is_semi_empirical(self.functional)
    
    def _pre_optimize_geometry(self, mol: Mol) -> Mol:
        """Pre-optimize molecule geometry using MMFF94 force field.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Molecule with optimized geometry (new conformer)
        """
        from rdkit.Chem import AllChem
        
        mol_copy = Mol(mol)
        
        if mol_copy.GetNumConformers() == 0:
            status = AllChem.EmbedMolecule(
                mol_copy,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )
            if status != 0:
                status = AllChem.EmbedMolecule(mol_copy)
                if status != 0:
                    logger.warning("Failed to embed molecule for MMFF94 optimization")
                    return mol_copy
        
        try:
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol_copy)
            if mmff_props is None:
                logger.warning("MMFF94 properties could not be initialized, skipping pre-optimization")
                return mol_copy
            
            conf = mol_copy.GetConformer()
            optimizer = AllChem.MMFFGetMoleculeForceField(mol_copy, mmff_props, confId=conf.GetId())
            
            if optimizer is None:
                logger.warning("MMFF94 force field could not be initialized, skipping pre-optimization")
                return mol_copy
            
            result = optimizer.Minimize(maxIts=1000)
            if result != 0:
                logger.debug(f"MMFF94 optimization converged with status: {result}")
            else:
                logger.debug("MMFF94 optimization completed successfully")
        except Exception as e:
            logger.warning(f"MMFF94 optimization failed: {e}, using original geometry")
        
        return mol_copy
    
    def _get_molecule_hash(self, mol: Mol) -> str:
        """Generate hash for molecule based on SMILES and calculation parameters.
        
        For semi-empirical methods, basis_set and dispersion_correction are
        not included in the hash since they are not used.
        """
        from rdkit.Chem import MolToSmiles
        
        smiles = MolToSmiles(mol, canonical=True)
        
        if self._is_semi_empirical():
            params = (
                self.functional,
                self.method_type,
                self.solvation_model,
                self.charge,
                self.multiplicity,
                self.pre_optimize,
            )
        else:
            params = (
                self.functional,
                self.basis_set,
                self.method_type,
                self.dispersion_correction,
                self.solvation_model,
                self.charge,
                self.multiplicity,
                self.pre_optimize,
            )
        
        key = f"{smiles}_{params}"
        return hashlib.sha256(key.encode()).hexdigest()
    
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
        for name in dir(self):
            if not name.startswith('_') and name not in excluded:
                method = getattr(self, name, None)
                if callable(method):
                    descriptor_methods.append(name)
        
        return sorted(descriptor_methods)
    
    @staticmethod
    def _format_time(seconds: float) -> str:
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
    
    @staticmethod
    def _update_time_estimates(actual_times: list[float], estimated_times: list[float], 
                               current_idx: int, total: int) -> None:
        """Update estimated times for remaining molecules based on actual performance.
        
        Args:
            actual_times: List of actual execution times for processed molecules
            estimated_times: List of estimated times (will be modified in place)
            current_idx: Current molecule index
            total: Total number of molecules
        """
        if not estimated_times or not actual_times:
            return
        
        avg_actual_time = sum(actual_times) / len(actual_times)
        
        ratios = [
            actual_times[k] / estimated_times[k] 
            for k in range(len(actual_times)) 
            if k < len(estimated_times) and estimated_times[k] > 0
        ]
        
        avg_ratio = sum(ratios) / len(ratios) if ratios else None
        
        for j in range(current_idx + 1, total):
            if j < len(estimated_times):
                original_estimate = estimated_times[j]
                if original_estimate > 0 and avg_ratio is not None:
                    estimated_times[j] = original_estimate * avg_ratio
                else:
                    estimated_times[j] = avg_actual_time
            elif j >= len(estimated_times):
                estimated_times.append(avg_actual_time)
    
    def run_benchmark(self, mol: Optional[Mol] = None) -> dict:
        """Run benchmark calculation to calibrate time estimation.
        
        Args:
            mol: Test molecule (default: benzene)
            
        Returns:
            Dictionary with benchmark data
        """
        if mol is None:
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
            use_mpirun=self.use_mpirun,
            mpirun_path=self.mpirun_path,
            extra_env=self.extra_env,
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

