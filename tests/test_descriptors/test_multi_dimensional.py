import pytest
from orca_descriptors import Molecule

# Structures and expected values for parametrization
# [SMILES, Name, Expected Dipole Moment (D), Expected PSA (Å²), Expected Nrot]

TEST_MOLECULES_DATA = [
    # 1. Benzene (Aromatic, Non-polar, Rigid)
    ("C1=CC=CC=C1", "Benzene", 0.0, 0.0, 0),
    # 2. Acetone (Polar, Rigid, Heteroatom O)
    ("CC(=O)C", "Acetone", 2.8, 17.07, 0), # Dipole ~2.8 D, PSA ~17.07 Å²
    # 3. Ethane (Non-polar, Flexible)
    ("CC", "Ethane", 0.0, 0.0, 1),
    # 4. Water (Highly polar, Heteroatom O)
    ("O", "Water", 1.85, 20.23, 0) # Dipole ~1.85 D, PSA ~20.23 Å²
]

# --- 1. Dipole Moment Test ---
@pytest.mark.parametrize(
    "smiles, name, expected_dipole, expected_psa, expected_nrot",
    TEST_MOLECULES_DATA
)
def test_dipole_moment(orca, smiles, name, expected_dipole, expected_psa, expected_nrot):
    """
    Dipole moment (mu_D) test for various molecules.
    Highly symmetric molecules (Benzene, Ethane) should have a near-zero dipole.
    Polar molecules (Acetone, Water) should have a significant dipole (> 1.0 D).
    """
    mol = Molecule.from_smiles(smiles)
    result = orca.dipole_moment(mol)
    
    # Assert: Checks that the result is close to the expected value (0.5 D tolerance)
    tolerance = 0.5 
    assert abs(result - expected_dipole) < tolerance, f"Dipole moment for {name} ({result:.2f} D) deviated significantly from {expected_dipole} D."
    
    # Additional check for polar/non-polar molecules
    if expected_dipole == 0.0:
        assert result < 0.1, f"Non-polar {name} has an unexpectedly high dipole: {result:.2f} D"
    else:
        assert result > 1.0, f"Polar {name} has an unexpectedly low dipole: {result:.2f} D"

# --- 2. Polar Surface Area (PSA) Test ---
@pytest.mark.parametrize(
    "smiles, name, expected_dipole, expected_psa, expected_nrot",
    TEST_MOLECULES_DATA
)
def test_polar_surface_area(orca, smiles, name, expected_dipole, expected_psa, expected_nrot):
    """
    Polar Surface Area (PSA) test. Molecules without N/O atoms (Benzene, Ethane) 
    should have a PSA of zero. Molecules with N/O (Acetone, Water) should have a non-zero PSA.
    """
    mol = Molecule.from_smiles(smiles)
    result = orca.polar_surface_area(mol)
    
    # Assert 1: PSA must be non-negative
    assert result >= 0.0, f"PSA for {name} is negative: {result:.2f} Å²"
    
    # Assert 2: Checks that the result is close to the expected value (2.0 Å² tolerance)
    tolerance = 2.0 
    assert abs(result - expected_psa) < tolerance, f"PSA for {name} ({result:.2f} Å²) deviated significantly from {expected_psa} Å²."
    
    # Additional check for zero values
    if expected_psa == 0.0:
        assert result < 0.01, f"Non-polar {name} has an unexpectedly high PSA: {result:.2f} Å²"
    else:
        assert result > 5.0, f"Polar {name} has an unexpectedly low PSA: {result:.2f} Å²"

# --- 3. Rotatable Bonds Test (Nrot) ---
@pytest.mark.parametrize(
    "smiles, name, expected_dipole, expected_psa, expected_nrot",
    TEST_MOLECULES_DATA
)
def test_rotatable_bonds(orca, smiles, name, expected_dipole, expected_psa, expected_nrot):
    """
    Number of rotatable bonds (Nrot) test. Ethane has one C-C rotatable bond (Nrot=1). 
    Benzene, Acetone, and Water are rigid (Nrot=0).
    """
    mol = Molecule.from_smiles(smiles)
    result = orca.num_rotatable_bonds(mol)
    
    # Assert 1: Nrot must be an integer
    assert isinstance(result, int), f"Nrot for {name} is not an integer: {result}"
    
    # Assert 2: Checks for exact match to the expected number
    assert result == expected_nrot, f"Nrot for {name} is {result}, but expected {expected_nrot}."