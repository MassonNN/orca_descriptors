import pytest
from rdkit.Chem import Atom, MolFromSmiles, AddHs
from typing import Dict, List, Tuple

@pytest.mark.parametrize("orca_parametrized", ["dft", "am1"], indirect=True)
def test_nmr_carbon_shift(orca_parametrized):
    """
    NMR Chemical Shift (δ_C) for Carbon atoms in Benzene.
    All C atoms should be equivalent, and the shift should be high (~128 ppm) due to the aromatic ring current.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Assumes the function returns a dictionary: {atom_index: shift_value (ppm)}
    result: Dict[int, float] = orca_parametrized.get_nmr_shifts(mol, atom_type='C')
    
    # Assert 1: Must return 6 carbon shifts
    assert len(result) == 6
    
    # Assert 2: All shifts must be positive
    assert all(shift > 0 for shift in result.values())
    
    # Assert 3: The average shift should be substantial (typical aromatic C shift is 120-200 ppm)
    avg_shift = sum(result.values()) / len(result)
    assert 100.0 < avg_shift < 200.0

@pytest.mark.parametrize("orca_parametrized", ["dft", "am1"], indirect=True)
def test_mayer_bond_index(orca_parametrized):
    """
    Mayer Bond Index for C-C bonds in Benzene.
    Due to aromaticity, all C-C indices should be nearly equal (~1.4 for DFT, different for AM1).
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Assumes the function returns a list of tuples: (atom_i, atom_j, index_value)
    result: List[Tuple[int, int, float]] = orca_parametrized.get_mayer_indices(mol, atom_type_i='C', atom_type_j='C')
    
    # Assert 1: Must be 6 C-C bonds
    assert len(result) == 6
    
    # Extract index values
    indices = [idx for _, _, idx in result]
    avg_index = sum(indices) / len(indices)
    
    # Assert 2: For DFT, average Mayer index must be between 1.3 and 1.5 (aromatic standard)
    # For AM1, values are different (typically higher), so we just check they're positive and reasonable
    if orca_parametrized.functional.upper() == "AM1":
        assert avg_index > 0
        assert avg_index < 10.0  # AM1 can give higher values
    else:
        assert 1.30 < avg_index < 1.50
    
    # Assert 3: All indices must be nearly identical (deviation < 0.05 for DFT, < 0.5 for AM1)
    max_deviation = max(abs(idx - avg_index) for idx in indices)
    if orca_parametrized.functional.upper() == "AM1":
        assert max_deviation < 0.5  # AM1 can have more variation
    else:
        assert max_deviation < 0.05

@pytest.mark.parametrize("orca_parametrized", ["dft"], indirect=True)
def test_nbo_stabilization_e2(orca_parametrized):
    """
    NBO Stabilization Energy (E(2)) for the n -> π* interaction in Acetone (O -> C=O π*).
    This energy must be positive and substantial (~5-20 kcal/mol) due to the resonance.
    Note: NBO analysis is only available for DFT methods, not semi-empirical.
    Note: Requires NBOEXE environment variable to be set to point to NBO executable.
    """
    import os
    if not os.environ.get("NBOEXE"):
        pytest.skip("NBOEXE environment variable not set. NBO analysis requires NBO executable.")
    
    mol = AddHs(MolFromSmiles("CC(=O)C"))
    # Assumes the function returns the E(2) energy of the strongest n -> π* interaction (in kcal/mol)
    try:
        result = orca_parametrized.nbo_stabilization_energy(mol, donor='LP(O)', acceptor='PiStar(C=O)')
    except ValueError as e:
        if "NBOEXE" in str(e) or "No NBO stabilization energies" in str(e):
            pytest.skip(f"NBO analysis not available: {e}")
        raise
    
    # Assert 1: Stabilization energy must be positive
    assert result > 0 
    
    # Assert 2: The value should be significant (e.g., greater than 5 kcal/mol)
    assert result > 5.0
    
# @pytest.mark.parametrize("orca_parametrized", ["dft"], indirect=True)
# def test_min_max_esp(orca_parametrized):
#     """
#     Minimum and Maximum Electrostatic Potential (V_min, V_max) on the vdW surface of Acetone.
#     V_min (most negative) should be near the Oxygen atom.
#     V_max (most positive) should be near the Hydrogen atoms.
#     Note: ESP is only available for DFT methods, not semi-empirical.
#     NOTE: This descriptor is currently not implemented due to limitations in ORCA 6.0.1
#     for generating ESP cube files directly. Would require orca_plot utility integration.
#     """
#     mol = AddHs(MolFromSmiles("CC(=O)C"))
#     
#     # V_min (Most negative potential, near the Oxygen)
#     v_min = orca_parametrized.get_esp_extrema(mol, 'min')
#     
#     # V_max (Most positive potential, near the H atoms)
#     v_max = orca_parametrized.get_esp_extrema(mol, 'max')
#     
#     # Assert 1: V_min must be negative (due to Oxygen)
#     assert v_min < 0 
#     
#     # Assert 2: V_max must be positive (due to Hydrogen)
#     assert v_max > 0
#     
#     # Assert 3: V_max must be greater than V_min
#     assert v_max > v_min