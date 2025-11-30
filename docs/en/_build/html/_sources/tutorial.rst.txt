Tutorial
========

This tutorial will guide you through using the ORCA Descriptors library both as a Python library and as a command-line tool.

Using as a Python Library
--------------------------

Basic Usage
~~~~~~~~~~~

First, import the necessary classes::

   from orca_descriptors import Orca
   from rdkit.Chem import MolFromSmiles, AddHs

Initialize the ORCA calculator with your preferred settings::

   orca = Orca(
       script_path="orca",
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       dispersion_correction="D3BJ",
       solvation_model="COSMO(Water)",
       n_processors=8,
   )

Create a molecule from a SMILES string::

   mol = AddHs(MolFromSmiles("C1=CC=CC=C1"))  # Benzene

Calculate descriptors::

   # Energy descriptors
   homo = orca.homo_energy(mol)
   lumo = orca.lumo_energy(mol)
   gap = orca.gap_energy(mol)
   
   # DFT descriptors
   mu = orca.ch_potential(mol)
   chi = orca.electronegativity(mol)
   eta = orca.abs_hardness(mol)
   
   # Thermodynamic descriptors
   energy = orca.total_energy(mol)
   gibbs = orca.gibbs_free_energy(mol)

Caching
~~~~~~~

The library automatically caches calculation results. If you calculate descriptors for the same molecule with the same parameters, it will use the cached result instead of running ORCA again::

   # First calculation - runs ORCA
   homo1 = orca.homo_energy(mol)  # Takes time
   
   # Second calculation - uses cache
   homo2 = orca.homo_energy(mol)  # Instant

Using as a Command-Line Tool
-----------------------------

The library can also be used as a command-line utility after installation.

Run Benchmark
~~~~~~~~~~~~~

Before estimating calculation times, run a benchmark::

   orca_descriptors run_benchmark --working_dir ./calculations

Estimate Calculation Time
~~~~~~~~~~~~~~~~~~~~~~~~~

Estimate how long a calculation will take::

   orca_descriptors approximate_time --molecule CCO --method_type Opt

All ORCA parameters are available as CLI arguments. For example::

   orca_descriptors approximate_time \\
       --molecule CCO \\
       --functional PBE0 \\
       --basis_set def2-TZVP \\
       --n_processors 4 \\
       --method_type Opt

Example Workflow
----------------

Here's a complete example of calculating descriptors for multiple molecules::

   from orca_descriptors import Orca
   from rdkit.Chem import MolFromSmiles, AddHs
   
   # Initialize calculator
   orca = Orca(
       working_dir="./calculations",
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )
   
   # List of molecules to process
   smiles_list = [
       "C1=CC=CC=C1",  # Benzene
       "CCO",           # Ethanol
       "CC(=O)C",       # Acetone
   ]
   
   results = []
   for smiles in smiles_list:
       mol = AddHs(MolFromSmiles(smiles))
       results.append({
           "smiles": smiles,
           "homo": orca.homo_energy(mol),
           "lumo": orca.lumo_energy(mol),
           "gap": orca.gap_energy(mol),
           "dipole": orca.dipole_moment(mol),
       })
   
   # Process results
   for r in results:
       print(f"{r['smiles']}: Gap = {r['gap']:.2f} eV")

