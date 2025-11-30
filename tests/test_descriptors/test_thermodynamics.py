from orca_descriptors.orca import Molecule

def test_gibbs_free_energy(orca):
    """
    Gibbs Free Energy (ΔG) is the total energy available for work and is the 
    criterion for spontaneous change. For a stable molecule like benzene, 
    the total Gibbs Free Energy in absolute units (Hartree) must be a large 
    negative number.
    """
    mol = Molecule.from_smiles("C1=CC=CC=C1")
    # Result is typically in Atomic Units (Hartree)
    result = orca.gibbs_free_energy(mol)
    
    # Assert 1: Gibbs Free Energy must be negative
    assert result < 0 
    
    # Assert 2: The absolute magnitude must be large (e.g., less than -200 Hartree for C6H6)
    assert result < -200.0

def test_entropy(orca):
    """
    Standard Entropy (ΔS) measures the molecular disorder. 
    It must always be a positive value.
    """
    mol = Molecule.from_smiles("C1=CC=CC=C1")
    # Result is typically in J/(mol*K) or cal/(mol*K)
    result = orca.entropy(mol)
    
    # Assert 1: Entropy must be positive
    assert result > 0 
    
    # Assert 2: Entropy for a gas-phase molecule like benzene should be substantial 
    # (e.g., greater than 100 J/(mol*K) or 20 cal/(mol*K))
    assert result > 10.0 # Using a low threshold for robustness

def test_molecular_volume(orca):
    """
    Molecular volume is a physical quantity reflecting the size of the molecule. 
    It must be a positive value.
    """
    mol = Molecule.from_smiles("C1=CC=CC=C1")
    # Result is typically in cubic Angstroms (Å^3)
    result = orca.molecular_volume(mol)
    
    # Assert 1: Volume must be positive
    assert result > 0 
    
    # Assert 2: Benzene is a medium-sized molecule; volume should be substantial (e.g., greater than 50 Å^3)
    assert result > 50.0

def test_enthalpy(orca):
    """
    Enthalpy (ΔH) includes the total electronic and nuclear energy plus thermal corrections 
    (translational, rotational, vibrational). For a stable molecule like benzene, 
    the total Enthalpy in absolute units (Hartree) must be a large negative number.
    """
    mol = Molecule.from_smiles("C1=CC=CC=C1")
    # Result is typically in Atomic Units (Hartree)
    result = orca.enthalpy(mol)
    
    # Assert 1: Enthalpy must be negative
    assert result < 0 
    
    # Assert 2: The absolute magnitude must be large (e.g., less than -200 Hartree for C6H6)
    assert result < -200.0

def test_ch_bond_lengths(orca):
    """
    In benzene, due to symmetry, all six C-H bond lengths must be 
    equal and close to the typical aromatic C-H bond length (approx. 1.08 - 1.10 Å).
    
    Assume orca.get_bond_lengths returns a list of tuples: (atom_index_1, atom_index_2, length)
    """
    mol = Molecule.from_smiles("C1=CC=CC=C1")
    # Get all C-H bond lengths (there are 6 of them)
    # The 'C' and 'H' arguments specify which atoms to search for bonds between.
    c_h_bond_lengths: list[tuple[int, int, float]] = orca.get_bond_lengths(mol, 'C', 'H')
    
    # Extract just the length values
    lengths = [length for _, _, length in c_h_bond_lengths]
    
    # Assert 1: There should be exactly 6 C-H bonds
    assert len(lengths) == 6
    
    # Assert 2: The average C-H bond length must be within the expected range (1.08 Å to 1.10 Å)
    avg_length = sum(lengths) / len(lengths)
    assert avg_length > 1.08 
    assert avg_length < 1.10 
    
    # Assert 3: All C-H bond lengths must be nearly identical (less than 0.005 Å difference from the mean)
    max_deviation = max(abs(l - avg_length) for l in lengths)
    assert max_deviation < 0.005

def test_cc_bond_lengths(orca):
    """
    In benzene, due to aromaticity, all six C-C bond lengths must be 
    approximately equal, intermediate between a single (1.54 Å) and double (1.34 Å) bond.
    
    Assume orca.get_bond_lengths returns a list of tuples: (atom_index_1, atom_index_2, length)
    """
    mol = Molecule.from_smiles("C1=CC=CC=C1")
    # Get all C-C bond lengths (there are 6 of them in the ring)
    c_c_bond_lengths: list[tuple[int, int, float]] = orca.get_bond_lengths(mol, 'C', 'C')
    
    # Extract just the length values
    lengths = [length for _, _, length in c_c_bond_lengths]
    
    # Assert 1: There should be exactly 6 C-C bonds in the ring
    assert len(lengths) == 6
    
    # Assert 2: The average C-C bond length must be close to the aromatic standard (approx. 1.39 Å)
    avg_length = sum(lengths) / len(lengths)
    assert avg_length > 1.36 
    assert avg_length < 1.42 
    
    # Assert 3: All C-C bond lengths must be nearly identical (less than 0.01 Å difference from the mean)
    max_deviation = max(abs(l - avg_length) for l in lengths)
    assert max_deviation < 0.01