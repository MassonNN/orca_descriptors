"""Test semi-empirical methods support (AM1, PM3, PM6, PM7, etc.)."""

import pytest
from rdkit.Chem import MolFromSmiles, AddHs

from orca_descriptors import Orca
from orca_descriptors.input_generator import ORCAInputGenerator


def test_semi_empirical_detection():
    """Test that semi-empirical methods are correctly detected."""
    generator = ORCAInputGenerator()
    
    # Test various semi-empirical methods
    assert generator._is_semi_empirical("AM1") is True
    assert generator._is_semi_empirical("PM3") is True
    assert generator._is_semi_empirical("PM6") is True
    assert generator._is_semi_empirical("PM7") is True
    assert generator._is_semi_empirical("RM1") is True
    assert generator._is_semi_empirical("MNDO") is True
    assert generator._is_semi_empirical("OM1") is True
    
    # Test case insensitivity
    assert generator._is_semi_empirical("am1") is True
    assert generator._is_semi_empirical("pm3") is True
    
    # Test DFT methods (should return False)
    assert generator._is_semi_empirical("PBE0") is False
    assert generator._is_semi_empirical("B3LYP") is False
    assert generator._is_semi_empirical("M06") is False


def test_semi_empirical_input_generation():
    """Test that input files for semi-empirical methods are generated correctly."""
    generator = ORCAInputGenerator()
    mol = AddHs(MolFromSmiles("CCO"))  # Ethanol
    
    # Generate input for AM1
    input_content = generator.generate(
        mol=mol,
        functional="AM1",
        basis_set="def2-SVP",  # Should be ignored
        method_type="Opt",
        dispersion_correction="D3BJ",  # Should be ignored
        n_processors=1,
    )
    
    # Check that input contains AM1 but not basis set or dispersion
    assert "AM1" in input_content
    assert "def2-SVP" not in input_content
    assert "D3BJ" not in input_content
    assert "! Opt" in input_content
    
    # Generate input for PM3
    input_content_pm3 = generator.generate(
        mol=mol,
        functional="PM3",
        basis_set="def2-TZVP",  # Should be ignored
        method_type="SP",
        n_processors=1,
    )
    
    assert "PM3" in input_content_pm3
    assert "def2-TZVP" not in input_content_pm3
    assert "! SP" in input_content_pm3


def test_semi_empirical_orca_instance():
    """Test that Orca instance correctly handles semi-empirical methods."""
    orca_am1 = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        basis_set="def2-SVP",  # Should be ignored but not cause errors
        method_type="SP",
        dispersion_correction="D3BJ",  # Should be ignored
        n_processors=1,
    )
    
    # Check that semi-empirical is detected
    assert orca_am1._is_semi_empirical() is True
    
    # Check that parameters are stored (even if ignored)
    assert orca_am1.functional == "AM1"
    assert orca_am1.basis_set == "def2-SVP"  # Stored but ignored in calculations
    assert orca_am1.dispersion_correction == "D3BJ"  # Stored but ignored


def test_semi_empirical_hash():
    """Test that hash calculation correctly handles semi-empirical methods."""
    mol = AddHs(MolFromSmiles("CCO"))
    
    # AM1 with different basis sets should have same hash (basis_set ignored)
    orca1 = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
    )
    
    orca2 = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        basis_set="def2-TZVP",  # Different basis set
        method_type="SP",
        n_processors=1,
    )
    
    hash1 = orca1._get_molecule_hash(mol)
    hash2 = orca2._get_molecule_hash(mol)
    
    # Hashes should be the same (basis_set is ignored for semi-empirical)
    assert hash1 == hash2, "Hash should be same for semi-empirical methods regardless of basis_set"
    
    # But different functionals should have different hashes
    orca3 = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PM3",  # Different method
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
    )
    
    hash3 = orca3._get_molecule_hash(mol)
    assert hash1 != hash3, "Hash should differ for different semi-empirical methods"


def test_semi_empirical_vs_dft_hash():
    """Test that semi-empirical and DFT methods have different hashes."""
    mol = AddHs(MolFromSmiles("CCO"))
    
    orca_semi = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
    )
    
    orca_dft = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
    )
    
    hash_semi = orca_semi._get_molecule_hash(mol)
    hash_dft = orca_dft._get_molecule_hash(mol)
    
    assert hash_semi != hash_dft, "Semi-empirical and DFT methods should have different hashes"


def test_semi_empirical_dispersion_ignored():
    """Test that dispersion correction is ignored for semi-empirical methods."""
    generator = ORCAInputGenerator()
    mol = AddHs(MolFromSmiles("CCO"))
    
    # Generate input with dispersion correction
    input_content = generator.generate(
        mol=mol,
        functional="PM6",
        basis_set="def2-SVP",
        method_type="Opt",
        dispersion_correction="D3BJ",  # Should be ignored
        n_processors=1,
    )
    
    # Dispersion should not appear in input
    assert "D3BJ" not in input_content
    assert "PM6" in input_content


def test_all_semi_empirical_methods():
    """Test that all supported semi-empirical methods are recognized."""
    generator = ORCAInputGenerator()
    
    supported_methods = ["AM1", "PM3", "PM6", "PM7", "RM1", "MNDO", "MNDOD", "OM1", "OM2", "OM3"]
    
    for method in supported_methods:
        assert generator._is_semi_empirical(method) is True, f"{method} should be recognized as semi-empirical"
        assert generator._is_semi_empirical(method.lower()) is True, f"{method.lower()} should be recognized (case insensitive)"


def test_semi_empirical_with_pre_optimize():
    """Test that pre_optimize works with semi-empirical methods."""
    orca = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        method_type="SP",
        n_processors=1,
        pre_optimize=True,
    )
    
    mol = AddHs(MolFromSmiles("CCO"))
    
    # Should be able to pre-optimize
    optimized = orca._pre_optimize_geometry(mol)
    assert optimized is not None
    assert optimized.GetNumConformers() > 0
    
    # Hash should include pre_optimize
    hash1 = orca._get_molecule_hash(mol)
    
    orca_no_opt = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        method_type="SP",
        n_processors=1,
        pre_optimize=False,
    )
    
    hash2 = orca_no_opt._get_molecule_hash(mol)
    assert hash1 != hash2, "Hash should differ when pre_optimize differs"


def test_semi_empirical_calculation():
    """Integration test: Run actual ORCA calculation with semi-empirical method (AM1)."""
    orca = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="AM1",
        method_type="SP",  # Single point for faster calculation
        n_processors=1,
        pre_optimize=True,
        log_level=10,  # DEBUG level
    )
    
    # Use a simple molecule (methane) for fast calculation
    mol = AddHs(MolFromSmiles("C"))
    
    # Check that it's detected as semi-empirical
    assert orca._is_semi_empirical() is True
    
    # Run a simple calculation - get total energy
    try:
        total_energy = orca.total_energy(mol)
        
        # Check that we got a result (should be a float)
        assert isinstance(total_energy, (int, float))
        assert total_energy != 0.0, "Total energy should not be zero"
        
        # For methane with AM1, energy should be negative (stable molecule)
        assert total_energy < 0, "Total energy should be negative for stable molecule"
        
        # Also test that we can get HOMO energy
        homo_energy = orca.homo_energy(mol)
        assert isinstance(homo_energy, (int, float))
        
    except Exception as e:
        pytest.fail(f"ORCA calculation with AM1 failed: {e}")

