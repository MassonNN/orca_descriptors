from pandas import DataFrame
from rdkit.Chem import MolFromSmiles, AddHs
from orca_descriptors import Orca


def test_pandas_compatibility(orca: Orca):
    """Test pandas compatibility with specific descriptors."""
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    df = DataFrame(["C1=CC=CC=C1", "CCO", "CC(=O)C"], columns=["smiles"])
    
    # Test with specific descriptors (subset)
    selected_descriptors = [
        'homo_energy', 'lumo_energy', 'gap_energy', 'ch_potential',
        'electronegativity', 'abs_hardness', 'abs_softness',
        'frontier_electron_density', 'entropy', 'enthalpy'
    ]
    df = orca.calculate_descriptors(
        smiles_column=df['smiles'],
        descriptors=selected_descriptors
    )

    assert df.loc[0, 'homo_energy'] == orca.homo_energy(mol)
    assert df.loc[0, 'lumo_energy'] == orca.lumo_energy(mol)
    assert df.loc[0, 'gap_energy'] == orca.gap_energy(mol)
    assert df.loc[0, 'ch_potential'] == orca.ch_potential(mol)
    assert df.loc[0, 'electronegativity'] == orca.electronegativity(mol)
    assert df.loc[0, 'abs_hardness'] == orca.abs_hardness(mol)
    assert df.loc[0, 'abs_softness'] == orca.abs_softness(mol)
    frontier_density = orca.frontier_electron_density(mol)
    if frontier_density:
        max_frontier = max(charge for _, charge in frontier_density)
        assert df.loc[0, 'frontier_electron_density'] == max_frontier
    else:
        assert df.loc[0, 'frontier_electron_density'] == 0.0
    assert df.loc[0, 'entropy'] == orca.entropy(mol)
    assert df.loc[0, 'enthalpy'] == orca.enthalpy(mol)
    
    # Verify that only selected descriptors are present
    for col in df.columns:
        if col not in ['smiles'] and col not in selected_descriptors:
            assert col not in df.columns, f"Unexpected descriptor column: {col}"


def test_pandas_compatibility_all_descriptors(orca: Orca):
    """Test pandas compatibility with all descriptors (default behavior)."""
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    df = DataFrame(["C1=CC=CC=C1"], columns=["smiles"])
    
    # Test with all descriptors (default, descriptors=None)
    df_all = orca.calculate_descriptors(smiles_column=df['smiles'])
    
    # Check that all expected descriptors are present
    expected_descriptors = [
        'homo_energy', 'lumo_energy', 'gap_energy', 'ch_potential',
        'electronegativity', 'abs_hardness', 'abs_softness',
        'frontier_electron_density', 'total_energy', 'dipole_moment',
        'polar_surface_area', 'gibbs_free_energy', 'entropy', 'enthalpy',
        'molecular_volume', 'num_rotatable_bonds', 'wiener_index',
        'solvent_accessible_surface_area', 'get_min_h_charge', 'xy_shadow',
        'meric', 'm_log_p', 'moran_autocorrelation', 'autocorrelation_hats'
    ]
    
    for desc in expected_descriptors:
        assert desc in df_all.columns, f"Descriptor '{desc}' not found in DataFrame"
    
    # Verify values match individual method calls
    assert df_all.loc[0, 'homo_energy'] == orca.homo_energy(mol)
    assert df_all.loc[0, 'total_energy'] == orca.total_energy(mol)
    assert df_all.loc[0, 'dipole_moment'] == orca.dipole_moment(mol)
    assert df_all.loc[0, 'num_rotatable_bonds'] == orca.num_rotatable_bonds(mol)


def test_pandas_compatibility_invalid_descriptor(orca: Orca):
    """Test that invalid descriptor names raise ValueError."""
    from pytest import raises
    
    df = DataFrame(["C1=CC=CC=C1"], columns=["smiles"])
    
    # Test with invalid descriptor name
    with raises(ValueError, match="Invalid descriptor names"):
        orca.calculate_descriptors(
            smiles_column=df['smiles'],
            descriptors=['homo_energy', 'invalid_descriptor', 'lumo_energy']
        )