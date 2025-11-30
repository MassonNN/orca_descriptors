from rdkit.Chem import Atom, MolFromSmiles, AddHs


def test_chem_potential(orca):
    """
    The chemical potential (mu) for benzene (a stable molecule) 
    should be negative, as it is defined as: mu = -electronegativity.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Call the function to perform the calculation and store the result
    result = orca.ch_potential(mol)
    
    # Assert 1: Chemical potential must be negative
    assert result < 0 
    
    # Assert 2: Check that the value is within a reasonable range (e.g., less than zero but greater than -10 eV)
    assert result > -10.0 


def test_electronegativity(orca):
    """
    Electronegativity (chi) is defined as: chi = -mu. 
    Therefore, it must be positive.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.electronegativity(mol)
    
    # Assert 1: Electronegativity must be positive
    assert result > 0 
    
    # Assert 2: The value should be in a typical range (e.g., greater than 1.0 eV)
    assert result > 1.0 


def test_abs_hardness(orca):
    """
    Absolute hardness (eta) is always positive. 
    For a stable molecule like benzene, it should be substantial.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.abs_hardness(mol)
    
    # Assert 1: Hardness is always positive
    assert result > 0 
    
    # Assert 2: Hardness for benzene should be greater than a minimum threshold (e.g., 1.0 eV)
    assert result > 1.0


def test_abs_softness(orca):
    """
    Absolute softness (S) is the reciprocal of hardness (S = 1/(2*eta)) 
    and is also always positive.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.abs_softness(mol)
    
    # Assert 1: Softness is always positive
    assert result > 0 
    
    # Assert 2: Softness for benzene (which is hard) should be less than 1.0
    assert result < 1.0 


def test_frontier_electron_density(orca):
    """
    Frontier electron density (HOMO/LUMO) indicates reaction centers.
    In benzene, the frontier orbitals are delocalized across all carbon atoms.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Assumes the function returns a list of tuples: (Atom, significant_frontier_electron_density)
    result: list[tuple[Atom, float]] = orca.frontier_electron_density(mol)  
    
    # Assert 1: A non-empty list of atoms should be returned (all carbon atoms participate)
    assert len(result) > 0
    
    # Assert 2: All atoms (the first element of the tuple) in the list must be carbon ('C') atoms
    assert all(atom_density_pair[0].GetSymbol() == 'C' for atom_density_pair in result)
    
    # Assert 3: In benzene, all 6 C atoms should have significant frontier electron density
    assert len(result) == 6

    assert all(atom_density_pair[1] > 0.001 for atom_density_pair in result)