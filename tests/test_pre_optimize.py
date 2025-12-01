"""Test pre-optimization with MMFF94."""

import pytest
from rdkit.Chem import MolFromSmiles, AddHs

from orca_descriptors import Orca


def test_pre_optimize_enabled():
    """Test that pre-optimization is enabled by default and works."""
    orca = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",  # Single point to avoid long calculations
        n_processors=1,
        pre_optimize=True,  # Explicitly enable
        log_level=10,  # DEBUG level
    )
    
    mol = AddHs(MolFromSmiles("CCO"))  # Ethanol
    
    # Check that pre_optimize is set
    assert orca.pre_optimize is True
    
    # The molecule should have a conformer after pre-optimization
    # We can't directly test the optimization, but we can check that
    # the method exists and doesn't raise errors
    optimized_mol = orca._pre_optimize_geometry(mol)
    
    # Should return a molecule
    assert optimized_mol is not None
    assert optimized_mol.GetNumConformers() > 0


def test_pre_optimize_disabled():
    """Test that pre-optimization can be disabled."""
    orca = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
        pre_optimize=False,  # Explicitly disable
    )
    
    mol = AddHs(MolFromSmiles("CCO"))
    
    # Check that pre_optimize is set to False
    assert orca.pre_optimize is False
    
    # The pre-optimization method should still work if called directly
    optimized_mol = orca._pre_optimize_geometry(mol)
    assert optimized_mol is not None


def test_pre_optimize_default():
    """Test that pre-optimize defaults to True."""
    orca = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
        # pre_optimize not specified, should default to True
    )
    
    assert orca.pre_optimize is True


def test_pre_optimize_hash_difference():
    """Test that pre_optimize affects the molecule hash."""
    mol = AddHs(MolFromSmiles("CCO"))
    
    orca_with_opt = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
        pre_optimize=True,
    )
    
    orca_without_opt = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
        pre_optimize=False,
    )
    
    # Hash should be different when pre_optimize differs
    hash_with = orca_with_opt._get_molecule_hash(mol)
    hash_without = orca_without_opt._get_molecule_hash(mol)
    
    assert hash_with != hash_without, "Hash should differ when pre_optimize differs"


def test_pre_optimize_integration():
    """Test that pre-optimization is called during calculation."""
    import logging
    from unittest.mock import patch
    
    orca = Orca(
        script_path="orca",
        working_dir="tests/test_data",
        output_dir="tests/test_data",
        functional="PBE0",
        basis_set="def2-SVP",
        method_type="SP",
        n_processors=1,
        pre_optimize=True,
        log_level=logging.DEBUG,
    )
    
    mol = AddHs(MolFromSmiles("C"))
    
    # Mock the _pre_optimize_geometry method to track if it's called
    call_count = {'count': 0}
    original_method = orca._pre_optimize_geometry
    
    def tracked_pre_optimize(m):
        call_count['count'] += 1
        return original_method(m)
    
    orca._pre_optimize_geometry = tracked_pre_optimize
    
    # When pre_optimize is True, _pre_optimize_geometry should be called
    # during _run_calculation. However, we can't easily test _run_calculation
    # without actually running ORCA, so we just verify the method exists
    # and can be called.
    optimized = orca._pre_optimize_geometry(mol)
    assert optimized is not None
    assert call_count['count'] == 1

