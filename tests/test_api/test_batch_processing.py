from orca_descriptors import Orca, ORCABatchProcessing
from rdkit.Chem import MolFromSmiles, AddHs


def test_batch_processing():
    orca = Orca(
        script_path="orca",
        working_dir=".",
        functional="PBE0",
        basis_set="def2-SVP",
    )
    batch_processing = ORCABatchProcessing(
        orca=orca,
        working_dir=".",
        functional="PBE0",
        basis_set="def2-SVP",
    )

    smiles_list = ["C1=CC=CC=C1", "CCO", "CC(=O)C"]
    x = batch_processing.x_molecule()  # mock molecule for define descriptors
    result = batch_processing.calculate_descriptors(
        smiles_list,
        descriptors=[
            orca.ch_potential(x),
            orca.electronegativity(x),
            orca.abs_hardness(x),
            orca.abs_softness(x),
            orca.frontier_electron_density(x),
            orca.topological_distance(x, 'O', 'O'),  # can use x molecule to define parameters
            orca.mo_energy(x, -3)
        ]
    )

    # result for each molecule in smiles_list with given descriptors parameters expected
    assert result is not None