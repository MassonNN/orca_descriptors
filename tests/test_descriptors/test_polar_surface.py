from rdkit.Chem import MolFromSmiles, AddHs


def test_dipole_moment(orca):
    """
    Dipole moment (mu_D) for highly symmetric, non-polar molecules like benzene 
    should be zero or extremely close to zero.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Result is typically in Debye (D)
    result = orca.dipole_moment(mol)
    
    # Assert 1: The dipole moment should be close to zero (e.g., less than 0.1 D)
    assert abs(result) < 0.1 
    
    # Assert 2: The result must be non-negative (for the magnitude of the dipole moment)
    assert result >= 0.0

def test_polar_surface_area(orca):
    """
    Polar Surface Area (PSA) is the surface area contributed by nitrogen and oxygen 
    atoms. Since benzene has none, PSA must be zero.
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Result is typically in square Angstroms (Å^2)
    result = orca.polar_surface_area(mol)
    
    # Assert 1: PSA should be exactly or very close to zero
    assert result < 0.01 
    
    # Assert 2: PSA cannot be negative
    assert result >= 0.0

def test_atom_charges(orca):
    """
    Atom charges (e.g., Mulliken) should be symmetrical across the molecule.
    C atoms should be slightly negative/neutral, H atoms slightly positive.
    The sum of all charges must be zero (for a neutral molecule).
    
    Assume orca.get_atom_charges returns a dictionary: {atom_index: charge_value}
    """
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    # Assume the function returns a mapping of atom index to charge
    charges: dict[int, float] = orca.get_atom_charges(mol)
    
    # Get the RDKit Atom objects (not strictly necessary for the asserts but good practice)
    atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
    
    # Assert 1: The sum of all charges must be zero (or close to zero)
    sum_of_charges = sum(charges.values())
    assert abs(sum_of_charges) < 0.01 # Tolerance for numerical errors
    
    # Assert 2: All C atoms (indexes 0-5) should have similar (slightly negative) charges
    # C atoms should be slightly negative (or very close to zero)
    c_charges = [charges[i] for i in range(6)]
    # Check if the average C charge is slightly negative (e.g., between -0.5 and 0.05)
    avg_c_charge = sum(c_charges) / len(c_charges)
    assert avg_c_charge < 0.05 
    assert avg_c_charge > -0.5 
    
    # Assert 3: All H atoms (indexes 6-11) should have similar positive charges
    h_charges = [charges[i] for i in range(6, 12)]
    # Check if the average H charge is positive (e.g., between 0.01 and 0.5)
    # Real ORCA calculations may give smaller charges (~0.04) for benzene
    avg_h_charge = sum(h_charges) / len(h_charges)
    assert avg_h_charge > 0.01  # Relaxed threshold for real calculations
    assert avg_h_charge < 0.5