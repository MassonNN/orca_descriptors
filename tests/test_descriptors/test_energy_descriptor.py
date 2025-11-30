from rdkit.Chem import MolFromSmiles, AddHs


def test_homo_energy(orca):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.homo_energy(mol)
    assert result < 0


def test_lumo_energy(orca):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.lumo_energy(mol)
    assert result < 5

def test_gap_energy(orca):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.gap_energy(mol)
    assert result > 1

def test_total_energy(orca):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca.total_energy(mol)
    assert result < -100


