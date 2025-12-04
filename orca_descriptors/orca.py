"""Main Orca class for quantum chemical calculations."""

import logging
from pathlib import Path
from typing import Any, Optional, Union

from rdkit.Chem import Mol, MolFromSmiles, AddHs

from orca_descriptors.base import OrcaBase
from orca_descriptors.cache import CacheManager
from orca_descriptors.calculation import CalculationMixin
from orca_descriptors.descriptors import (
    ElectronicDescriptorsMixin,
    EnergyDescriptorsMixin,
    StructuralDescriptorsMixin,
    TopologicalDescriptorsMixin,
    MiscDescriptorsMixin,
)
from orca_descriptors.input_generator import ORCAInputGenerator
from orca_descriptors.output_parser import ORCAOutputParser
from orca_descriptors.time_estimator import ORCATimeEstimator

logger = logging.getLogger(__name__)


class Orca(
    OrcaBase,
    CalculationMixin,
    ElectronicDescriptorsMixin,
    EnergyDescriptorsMixin,
    StructuralDescriptorsMixin,
    TopologicalDescriptorsMixin,
    MiscDescriptorsMixin,
):
    """Main class for ORCA quantum chemical calculations.
    
    Supports both DFT methods (with basis sets) and semi-empirical methods
    (AM1, PM3, PM6, PM7, etc.). For semi-empirical methods, basis_set and
    dispersion_correction parameters are ignored.
    """
    
    def __init__(
        self,
        script_path: str = "orca",
        working_dir: str = ".",
        output_dir: str = ".",
        functional: str = "AM1",
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
        pre_optimize: bool = True,
        cache_server_url: Optional[str] = None,
        cache_api_token: Optional[str] = None,
        cache_timeout: int = 30,
    ):
        """Initialize ORCA calculator.
        
        Args:
            script_path: Path to ORCA executable
            working_dir: Working directory for calculations
            output_dir: Directory for output files
            functional: DFT functional (e.g., "PBE0") or semi-empirical method
                       (e.g., "AM1", "PM3", "PM6", "PM7"). For semi-empirical methods,
                       basis_set and dispersion_correction are ignored.
            basis_set: Basis set (e.g., "def2-SVP"). Ignored for semi-empirical methods.
            method_type: Calculation type ("Opt", "SP", etc.)
            dispersion_correction: Dispersion correction (e.g., "D3BJ").
                                 Ignored for semi-empirical methods.
            solvation_model: Solvation model (e.g., "COSMO(Water)")
            n_processors: Number of processors
            max_scf_cycles: Maximum SCF cycles
            scf_convergence: SCF convergence threshold
            charge: Molecular charge
            multiplicity: Spin multiplicity
            cache_dir: Directory for caching results (default: output_dir/.orca_cache)
            log_level: Logging level (default: logging.INFO)
            max_wait: Maximum time to wait for output file creation in seconds (default: 300)
            use_mpirun: Whether to use mpirun for parallel execution (default: False)
            mpirun_path: Path to mpirun executable (default: None, will search in PATH)
            extra_env: Additional environment variables to pass to ORCA process (default: None)
            pre_optimize: Whether to pre-optimize geometry with MMFF94 before ORCA calculation (default: True)
            cache_server_url: URL of the remote cache server (e.g., "http://localhost:3000")
            cache_api_token: API token for remote cache authentication
            cache_timeout: Timeout for remote cache requests in seconds (default: 30)
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
        self.use_mpirun = use_mpirun
        self.mpirun_path = mpirun_path
        self.extra_env = extra_env or {}
        self.pre_optimize = pre_optimize
        
        self.working_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize remote cache client if credentials provided
        remote_cache_client = None
        if cache_server_url and cache_api_token:
            try:
                from orca_descriptors.remote_cache import RemoteCacheClient
                remote_cache_client = RemoteCacheClient(
                    server_url=cache_server_url,
                    api_token=cache_api_token,
                    timeout=cache_timeout,
                )
                logger.info(f"Remote cache enabled: {cache_server_url}")
            except Exception as e:
                logger.warning(
                    f"Failed to initialize remote cache client: {e}. "
                    f"Continuing with local cache only."
                )
        elif cache_server_url or cache_api_token:
            logger.warning(
                "Both cache_server_url and cache_api_token must be provided "
                "to enable remote caching. Continuing with local cache only."
            )
        
        cache_dir = cache_dir or str(self.output_dir / ".orca_cache")
        self.cache = CacheManager(cache_dir, remote_cache_client=remote_cache_client)
        self.input_generator = ORCAInputGenerator()
        self.output_parser = ORCAOutputParser()
        self.time_estimator = ORCATimeEstimator(working_dir=self.working_dir)
    
    def calculate_descriptors(
        self,
        smiles_column: Union[Any, list[str]],
        descriptors: Optional[list[str]] = None,
        progress: bool = True,
    ) -> Any:
        """Calculate descriptors for molecules from SMILES and add to DataFrame.
        
        This method provides optional pandas compatibility. If pandas is available,
        it accepts a pandas Series or DataFrame column and returns a DataFrame
        with added descriptor columns. If pandas is not available, it accepts
        a list of SMILES strings and returns a list of dictionaries.
        
        By default, calculates all available descriptors. Use the `descriptors`
        parameter to specify a subset of descriptors to calculate.
        
        Note: This method is a wrapper around ORCABatchProcessing for backward compatibility.
        For advanced features like multiprocessing, use ORCABatchProcessing directly.
        
        Args:
            smiles_column: pandas Series/DataFrame column with SMILES strings,
                          or a list of SMILES strings
            descriptors: Optional list of descriptor names to calculate.
                       If None, calculates all available descriptors.
                       Descriptor names correspond to method names of the Orca class.
            progress: Whether to show progress (default: True)
            
        Returns:
            DataFrame with descriptor columns added (if pandas available),
            or list of dictionaries (if pandas not available)
            
        Raises:
            ImportError: If pandas is not installed and a pandas object is passed
            ValueError: If an invalid descriptor name is provided
        """
        from orca_descriptors.batch_processing import ORCABatchProcessing
        
        batch_processor = ORCABatchProcessing(orca=self, parallel_mode="sequential")
        
        result = batch_processor.calculate_descriptors(
            smiles_column=smiles_column,
            descriptors=descriptors,
            progress=progress,
        )
        
        try:
            import pandas as pd
            if isinstance(result, pd.DataFrame):
                if isinstance(smiles_column, pd.Series):
                    result.insert(0, 'smiles', smiles_column.values)
                elif isinstance(smiles_column, pd.DataFrame):
                    pass
        except ImportError:
            pass
        
        return result

