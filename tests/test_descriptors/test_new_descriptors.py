import pytest
from rdkit.Chem import MolFromSmiles, AddHs

# Parametrization for E_HOMO-1
# [SMILES, Name, Expected E_HOMO-1 (eV)]
HOMO_TEST_DATA = [
    ("C1=CC=CC=C1", "Benzene", -0.09), # E_HOMO-1 for benzene (second-to-last occupied orbital)
]

@pytest.mark.parametrize("smiles, name, expected_homo_minus_1", HOMO_TEST_DATA)
def test_homo_minus_one_energy(orca, smiles, name, expected_homo_minus_1):
    """
    HOMO-1 Energy (E_HOMO-1) test. Should be negative and slightly lower than E_HOMO.
    """
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.mo_energy(mol, index) returns the MO energy by index
    result = orca.mo_energy(mol, index=-2) # -2 usually corresponds to HOMO-1
    
    # Assert 1: Energy must be nega tive
    assert result < 0 
    
    # Assert 2: Checks that the value is within the expected range (0.1 eV tolerance)
    tolerance = 0.1
    assert abs(result - expected_homo_minus_1) < tolerance, f"E_HOMO-1 for {name} ({result:.2f} eV) deviated significantly from {expected_homo_minus_1} eV."

# Parametrization for ESP-MNACH
# [SMILES, Name, Expected ESP-MNACH]
ESP_TEST_DATA = [
    ("CC(=O)C", "Acetone", 0.056), # Minimum H-charge in Acetone (H in CH3)
]

@pytest.mark.parametrize("smiles, name, expected_min_h_charge", ESP_TEST_DATA)
def test_esp_mnach(orca, smiles, name, expected_min_h_charge):
    """
    ESP-MNACH (Minimum Net Atomic Charge for Hydrogen) test. 
    Hydrogen atoms typically have positive charge.
    """
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.get_min_h_charge(mol, method="ESP") returns the minimum H charge
    result = orca.get_min_h_charge(mol, method="ESP") 
    
    # Assert 1: Minimum H charge must be positive
    assert result > 0 
    
    # Assert 2: Checks that the value is within the expected range (0.02 e tolerance)
    tolerance = 0.02
    assert abs(result - expected_min_h_charge) < tolerance, f"ESP-MNACH for {name} ({result:.2f} e) deviated significantly from {expected_min_h_charge} e."

# --- 2. Geometric and Reactivity Descriptors ---

# Parametrization for XY Shadow and MERIC
# [SMILES, Name, Expected XY Shadow (Å²), Expected MERIC (eV)]
GEOM_REACTIVITY_DATA = [
    ("C1=CC=CC=C1", "Benzene", 18.0, None), # XY Shadow: ~18 Å²
    ("O(CCO)C", "Ethylene Glycol", 12.0, -0.82), # MERIC (for C-O)
]

@pytest.mark.parametrize("smiles, name, expected_xy_shadow, expected_meric", GEOM_REACTIVITY_DATA)
def test_xy_shadow(orca, smiles, name, expected_xy_shadow, expected_meric):
    """
    XY Shadow/XY Rectangle test. Measures the area of the smallest rectangle 
    enclosing the molecule's projection onto the XY plane (3D descriptor).
    """
    if expected_xy_shadow is None:
        pytest.skip("XY Shadow not defined for this molecule in test data.")
        
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.xy_shadow(mol) returns the area in Å²
    result = orca.xy_shadow(mol)
    
    # Assert 1: Area must be positive
    assert result > 0 
    
    # Assert 2: Checks that the value is within the expected range (5.0 Å² tolerance)
    tolerance = 5.0
    assert abs(result - expected_xy_shadow) < tolerance, f"XY Shadow for {name} ({result:.2f} Å²) deviated significantly from {expected_xy_shadow} Å²."

@pytest.mark.parametrize("smiles, name, expected_xy_shadow, expected_meric", GEOM_REACTIVITY_DATA)
def test_meric(orca, smiles, name, expected_xy_shadow, expected_meric):
    """
    MERIC (Electrophilicity Index for Carbon) test. 
    Predicts sites susceptible to nucleophilic attack (should be negative for carbonyl/carboxyl C).
    """
    if expected_meric is None:
        pytest.skip("MERIC not defined for this molecule in test data.")
        
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.meric(mol) returns the minimum MERIC for C atoms (e.g., C in C-O)
    result = orca.meric(mol) 
    
    # Assert 1: For carbon bonded to a heteroatom, MERIC must be negative
    assert result < 0 
    
    # Assert 2: Checks that the value is within the expected range (0.2 eV tolerance)
    tolerance = 0.2
    assert abs(result - expected_meric) < tolerance, f"MERIC for {name} ({result:.2f} eV) deviated significantly from {expected_meric} eV."

# --- 3. Topological and Physicochemical Descriptors ---

# Parametrization for T(O...O) and M log P
# [SMILES, Name, Expected T(O...O), Expected M log P]
TOPO_PHYSCHEM_DATA = [
    ("O(CCO)C", "Ethylene Glycol", 3, -0.37), # T(O...O)=3 (O-C-C-O), M log P (RDKit Crippen)
    ("CCCCCCCC", "Octane", None, 3.37), # Octane, M log P (RDKit Crippen)
]

@pytest.mark.parametrize("smiles, name, expected_t_oo, expected_m_logp", TOPO_PHYSCHEM_DATA)
def test_t_oo(orca, smiles, name, expected_t_oo, expected_m_logp):
    """
    T(O...O) test. Sum of topological distances between all pairs of oxygen atoms. 
    Ethylene glycol has two O atoms with a distance of 3 (O-C-C-O).
    """
    if expected_t_oo is None:
        pytest.skip("T(O...O) not defined for this molecule in test data.")

    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.topological_distance(mol, 'O', 'O') returns the sum of distances
    result = orca.topological_distance(mol, 'O', 'O')
    
    # Assert 1: Result must be an integer
    assert isinstance(result, int) 
    
    # Assert 2: Checks for exact match
    assert result == expected_t_oo, f"T(O...O) for {name} is {result}, but expected {expected_t_oo}."

@pytest.mark.parametrize("smiles, name, expected_t_oo, expected_m_logp", TOPO_PHYSCHEM_DATA)
def test_m_log_p(orca, smiles, name, expected_t_oo, expected_m_logp):
    """
    M log P test (Moriguchi Log P). Measures the octanol/water partition coefficient. 
    Should be negative for polar/hydrophilic molecules (Ethylene Glycol) and positive for non-polar/hydrophobic (Octane).
    """
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.m_log_p(mol) returns the partition coefficient
    result = orca.m_log_p(mol)
    
    # Assert 1: Checks the sign
    if expected_m_logp > 0:
        assert result > 0, f"M log P for hydrophobic {name} should be positive, but got {result:.2f}."
    else:
        assert result < 0, f"M log P for hydrophilic {name} should be negative, but got {result:.2f}."
        
    # Assert 2: Checks that the value is within the expected range (0.05 tolerance)
    tolerance = 0.05
    assert abs(result - expected_m_logp) < tolerance, f"M log P for {name} ({result:.2f}) deviated significantly from {expected_m_logp}."

# --- 4. Autocorrelation Tests (MATS2v, HATS4u) ---

# Parametrization for autocorrelations
# [SMILES, Name, Expected MATS2v (range), Expected HATS4u (range)]
AUTOCORR_DATA = [
    ("CCCCCCCC", "Octane", -0.01, 0.05), # Octane: MATS2v (variance=0 when all atoms have same vdW volume), HATS4u ~0.05
]

@pytest.mark.parametrize("smiles, name, expected_mats2v, expected_hats4u", AUTOCORR_DATA)
def test_mats2v(orca, smiles, name, expected_mats2v, expected_hats4u):
    """
    MATS2v (Moran Autocorrelation - Lag 2, weighted by vdW volume) test. 
    Autocorrelation descriptors are typically centered around zero.
    """
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.moran_autocorrelation(mol, lag=2, weight='vdw_volume') returns the descriptor
    result = orca.moran_autocorrelation(mol, lag=2, weight='vdw_volume')
    
    # Assert 1: Value must be in a reasonable range (e.g., -1.0 to 1.0)
    assert -1.0 <= result <= 1.0 
    
    # Assert 2: Checks that the value is within the expected range (0.01 tolerance)
    tolerance = 0.01
    assert abs(result - expected_mats2v) < tolerance, f"MATS2v for {name} ({result:.2f}) deviated significantly from {expected_mats2v}."

@pytest.mark.parametrize("smiles, name, expected_mats2v, expected_hats4u", AUTOCORR_DATA)
def test_hats4u(orca, smiles, name, expected_mats2v, expected_hats4u):
    """
    HATS4u (Autocorrelation H-bond Acceptor/Donor - Lag 4, unweighted) test.
    This descriptor is often close to zero.
    """
    mol = AddHs(MolFromSmiles(smiles))
    # Assumes orca.autocorrelation_hats(mol, lag=4, unweighted=True) returns the descriptor
    result = orca.autocorrelation_hats(mol, lag=4, unweighted=True)
    
    # Assert 1: Value must be in a reasonable range (e.g., -1.0 to 1.0)
    assert -1.0 <= result <= 1.0 
    
    # Assert 2: Checks that the value is within the expected range (0.1 tolerance)
    tolerance = 0.1
    assert abs(result - expected_hats4u) < tolerance, f"HATS4u for {name} ({result:.2f}) deviated significantly from {expected_hats4u}."