import logging
import sys

import pytest

from orca_descriptors import Orca

# Configure logging for tests
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout,
    force=True
)

# Set up logger for orca_descriptors
orca_logger = logging.getLogger('orca_descriptors')
orca_logger.setLevel(logging.INFO)


@pytest.fixture(scope="session", autouse=True)
def setup_logging():
    """Setup logging for all tests."""
    yield
    # Cleanup if needed


@pytest.fixture
def orca():
    return Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="Opt",  # Optimization
        dispersion_correction="D3BJ",  # Dispersion correction
        solvation_model="COSMO(Water)",  # Solvation model
        n_processors=1,  # Use 1 processor to avoid MPI issues in tests
        max_scf_cycles=100,
        scf_convergence=1e-6,
        charge=0,
        multiplicity=1,
        log_level=logging.INFO,
    )


@pytest.fixture
def orca_gas():
    """ORCA instance without solvation model for gas-phase calculations."""
    return Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="Opt",  # Optimization
        dispersion_correction="D3BJ",  # Dispersion correction
        solvation_model=None,  # No solvation model for gas-phase calculations
        n_processors=1,  # Use 1 processor to avoid MPI issues in tests
        max_scf_cycles=100,
        scf_convergence=1e-6,
        charge=0,
        multiplicity=1,
        log_level=logging.INFO,
    )


@pytest.fixture
def orca_am1():
    """ORCA instance with AM1 semi-empirical method."""
    return Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        method_type="SP",  # Single point for faster calculations
        n_processors=1,
        pre_optimize=True,
        log_level=logging.INFO,
    )


@pytest.fixture(params=["dft", "am1"])
def orca_parametrized(request):
    """Parametrized ORCA fixture for testing with different methods.
    
    Parameters:
        - "dft": PBE0/def2-SVP with D3BJ dispersion
        - "am1": AM1 semi-empirical method
    """
    if request.param == "dft":
        return Orca(
            script_path="orca",
            working_dir="tests/test_data",
            output_dir="tests/test_data",
            functional="PBE0",
            basis_set="def2-SVP",
            method_type="SP",  # Single point for faster tests
            dispersion_correction="D3BJ",
            n_processors=1,
            pre_optimize=True,
            log_level=logging.INFO,
        )
    elif request.param == "am1":
        return Orca(
            script_path="orca",
            working_dir="tests/test_data",
            output_dir="tests/test_data",
            functional="AM1",
            method_type="SP",  # Single point for faster calculations
            n_processors=1,
            pre_optimize=True,
            log_level=logging.INFO,
        )

