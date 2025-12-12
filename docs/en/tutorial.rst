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
       pre_optimize=True,  # Enable geometry pre-optimization with MMFF94
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

Remote Cache
~~~~~~~~~~~~

The library supports remote caching via API, allowing you to share calculation results across different machines and users. A public cache server is available at ``https://api.orca-descriptors.massonnn.ru``.

To use remote cache, you need to:

1. Register at `orca-descriptors.massonnn.ru <https://orca-descriptors.massonnn.ru>`_ and generate an API token
2. Provide the API token when initializing the Orca instance::

   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       cache_api_token="your_api_token_here",
       # Optional: specify custom server URL
       # cache_server_url="https://api.orca-descriptors.massonnn.ru",
   )

The remote cache works transparently:
- Local cache is checked first, then remote cache if available
- Results found in remote cache are automatically downloaded and stored locally
- New calculations are automatically uploaded to the remote server (if you have upload permissions)
- If remote cache is unavailable, the system falls back to local-only caching

Cache-Only Mode
~~~~~~~~~~~~~~~

You can enable cache-only mode to use only cached results without running ORCA calculations. This is useful for quickly retrieving results from cache::

   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       cache_api_token="your_api_token_here",
       cache_only=True,  # Only use cache, don't run calculations
   )
   
   # If result is in cache, returns value
   # If result is not in cache, raises FileNotFoundError
   homo = orca.homo_energy(mol)

In cache-only mode:
- Only cached results (local or remote) are used
- If a result is not found in cache, descriptors return ``None`` in batch processing
- No ORCA calculations are performed
- Useful for quickly retrieving results from cache without running expensive calculations

Choosing Functionals and Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The library supports both DFT methods and semi-empirical methods:

DFT Methods
^^^^^^^^^^^

For high-accuracy calculations, use DFT functionals with basis sets::

   # High-accuracy DFT calculation
   orca_dft = Orca(
       functional="PBE0",           # Hybrid functional
       basis_set="def2-TZVP",      # Triple-zeta basis set
       method_type="Opt",           # Geometry optimization
       dispersion_correction="D3BJ", # Dispersion correction
       n_processors=8,
   )

Common DFT functionals:

- ``PBE0`` - Hybrid GGA functional, good balance of accuracy and speed
- ``B3LYP`` - Popular hybrid functional
- ``M06-2X`` - Meta-GGA functional, good for thermochemistry
- ``ωB97X-D`` - Range-separated hybrid with dispersion

Common basis sets:

- ``def2-SVP`` - Small, fast (default)
- ``def2-TZVP`` - Triple-zeta, more accurate
- ``def2-QZVP`` - Quadruple-zeta, very accurate but slow

Semi-Empirical Methods
^^^^^^^^^^^^^^^^^^^^^^

For faster calculations on large molecules, use semi-empirical methods::

   # Fast semi-empirical calculation
   orca_semi = Orca(
       functional="AM1",            # Semi-empirical method
       method_type="SP",            # Single point (no optimization)
       n_processors=1,
       pre_optimize=True,           # Pre-optimize with MMFF94
   )

Supported semi-empirical methods:

- ``AM1`` - Austin Model 1, good for organic molecules
- ``PM3`` - Parametric Method 3, improved over AM1
- ``PM6`` - Parametric Method 6, better for transition metals
- ``PM7`` - Parametric Method 7, improved accuracy
- ``RM1`` - Recife Model 1, optimized for organic compounds

Note: For semi-empirical methods, ``basis_set`` and ``dispersion_correction`` parameters are automatically ignored.

Geometry Pre-Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the library performs geometry pre-optimization using RDKit's MMFF94 force field before sending the molecule to ORCA. This can significantly speed up ORCA calculations, especially for geometry optimizations::

   # With pre-optimization (default)
   orca = Orca(
       functional="PBE0",
       method_type="Opt",
       pre_optimize=True,  # Default: True
   )
   
   # Without pre-optimization
   orca = Orca(
       functional="PBE0",
       method_type="Opt",
       pre_optimize=False,
   )

Benefits of pre-optimization:

- Faster ORCA convergence (fewer optimization steps needed)
- More stable calculations (better starting geometry)
- Reduced computational cost

The pre-optimization uses MMFF94, which requires explicit hydrogens. The library automatically adds hydrogens if needed and generates 3D coordinates if not present.

Method Types
~~~~~~~~~~~~

Choose the appropriate calculation type based on your needs::

   # Single point energy calculation (fastest)
   orca_sp = Orca(
       functional="PBE0",
       method_type="SP",  # Single point
   )
   
   # Geometry optimization (slower, but provides optimized geometry)
   orca_opt = Orca(
       functional="PBE0",
       method_type="Opt",  # Optimization
   )

Note: Some descriptors (like ``molecular_volume``, ``polar_surface_area``, ``solvent_accessible_surface_area``) require optimized geometries and may not work correctly with ``method_type="SP"``.

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

