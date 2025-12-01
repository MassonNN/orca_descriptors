import pytest
from rdkit.Chem import MolFromSmiles, AddHs


@pytest.mark.parametrize("orca_parametrized", ["dft", "am1"], indirect=True)
def test_homo_energy(orca_parametrized):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca_parametrized.homo_energy(mol)
    assert result < 0


@pytest.mark.parametrize("orca_parametrized", ["dft", "am1"], indirect=True)
def test_lumo_energy(orca_parametrized):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca_parametrized.lumo_energy(mol)
    assert result < 5


@pytest.mark.parametrize("orca_parametrized", ["dft", "am1"], indirect=True)
def test_gap_energy(orca_parametrized):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca_parametrized.gap_energy(mol)
    assert result > 1


@pytest.mark.parametrize("orca_parametrized", ["dft", "am1"], indirect=True)
def test_total_energy(orca_parametrized):
    mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))
    result = orca_parametrized.total_energy(mol)
    # For AM1, energy is less negative than DFT, so adjust threshold
    if orca_parametrized.functional == "AM1":
        assert result < -5  # AM1 energies are less negative
    else:
        assert result < -100  # DFT energies are more negative


