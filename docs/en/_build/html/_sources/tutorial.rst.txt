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
   
   # Molecular orbital descriptors
   homo_minus_1 = orca.mo_energy(mol, index=-2)  # HOMO-1 energy
   
   # Charge descriptors
   min_h_charge = orca.get_min_h_charge(mol, method="ESP")  # Minimum H charge
   
   # Geometric descriptors
   xy_area = orca.xy_shadow(mol)  # XY projection area
   
   # Reactivity descriptors
   meric = orca.meric(mol)  # Electrophilicity index for carbon
   
   # Topological descriptors
   t_oo = orca.topological_distance(mol, 'O', 'O')  # Sum of O-O distances
   nrot = orca.num_rotatable_bonds(mol)  # Number of rotatable bonds
   wiener = orca.wiener_index(mol)  # Wiener index
   
   # Physicochemical descriptors
   logp = orca.m_log_p(mol)  # Octanol/water partition coefficient
   sasa = orca.solvent_accessible_surface_area(mol)  # SASA
   
   # Autocorrelation descriptors
   mats2v = orca.moran_autocorrelation(mol, lag=2, weight='vdw_volume')
   hats4u = orca.autocorrelation_hats(mol, lag=4, unweighted=True)

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

Available Descriptors
---------------------

The library provides a comprehensive set of descriptors for QSAR analysis:

Energy Descriptors
~~~~~~~~~~~~~~~~~~

- ``homo_energy(mol)`` - HOMO energy (eV)
- ``lumo_energy(mol)`` - LUMO energy (eV)
- ``gap_energy(mol)`` - HOMO-LUMO gap (eV)
- ``mo_energy(mol, index)`` - Molecular orbital energy by index (eV)
- ``total_energy(mol)`` - Total energy (Hartree)

DFT Descriptors
~~~~~~~~~~~~~~~

- ``ch_potential(mol)`` - Chemical potential (eV)
- ``electronegativity(mol)`` - Electronegativity (eV)
- ``abs_hardness(mol)`` - Absolute hardness (eV)
- ``abs_softness(mol)`` - Absolute softness (1/eV)
- ``frontier_electron_density(mol)`` - Frontier electron density

Charge Descriptors
~~~~~~~~~~~~~~~~~~

- ``get_atom_charges(mol)`` - Mulliken atomic charges
- ``get_min_h_charge(mol, method="ESP")`` - Minimum hydrogen charge

Geometric Descriptors
~~~~~~~~~~~~~~~~~~~~~

- ``xy_shadow(mol)`` - XY projection area (Å²)
- ``molecular_volume(mol)`` - Molecular volume (Å³)
- ``get_bond_lengths(mol, atom1, atom2)`` - Bond lengths (Å)

Reactivity Descriptors
~~~~~~~~~~~~~~~~~~~~~~

- ``meric(mol)`` - Minimum electrophilicity index for carbon (eV)
- ``dipole_moment(mol)`` - Dipole moment (Debye)

Topological Descriptors
~~~~~~~~~~~~~~~~~~~~~~~

- ``topological_distance(mol, atom1, atom2)`` - Sum of topological distances
- ``num_rotatable_bonds(mol)`` - Number of rotatable bonds
- ``wiener_index(mol)`` - Wiener index

Physicochemical Descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``m_log_p(mol)`` - Moriguchi Log P (octanol/water partition coefficient)
- ``polar_surface_area(mol)`` - Polar surface area (Å²)
- ``solvent_accessible_surface_area(mol)`` - SASA (Å²)

Thermodynamic Descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``gibbs_free_energy(mol)`` - Gibbs free energy (Hartree)
- ``entropy(mol)`` - Entropy (J/(mol·K))
- ``enthalpy(mol)`` - Enthalpy (Hartree)

Autocorrelation Descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``moran_autocorrelation(mol, lag, weight)`` - Moran autocorrelation
- ``autocorrelation_hats(mol, lag, unweighted)`` - HATS autocorrelation

