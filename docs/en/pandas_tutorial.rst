Pandas Integration Tutorial
============================

This tutorial demonstrates how to use the ORCA Descriptors library with pandas for efficient batch processing of molecular descriptors.

Installation
------------

To use pandas integration, install pandas as an optional dependency::

   pip install 'orca-descriptors[pandas]'

Or install pandas separately::

   pip install pandas

Basic Usage
-----------

The ``calculate_descriptors`` method provides seamless integration with pandas DataFrames and Series.

Calculating All Descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ``calculate_descriptors`` calculates all available descriptors for each molecule::

   from orca_descriptors import Orca
   import pandas as pd

   # Initialize ORCA calculator
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

   # Create a DataFrame with SMILES strings
   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'CCO', 'CC(=O)C'],
       'name': ['Benzene', 'Ethanol', 'Acetone']
   })

   # Calculate all descriptors
   df_with_descriptors = orca.calculate_descriptors(smiles_column=df['smiles'])

   # The DataFrame now contains all descriptor columns
   print(df_with_descriptors.columns)
   # Output: Index(['smiles', 'name', 'homo_energy', 'lumo_energy', 'gap_energy', ...])

   # Access descriptor values
   print(df_with_descriptors[['name', 'homo_energy', 'lumo_energy', 'gap_energy']])

Selecting Specific Descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can specify which descriptors to calculate using the ``descriptors`` parameter::

   # Calculate only specific descriptors
   selected_descriptors = [
       'homo_energy',
       'lumo_energy',
       'gap_energy',
       'dipole_moment',
       'molecular_volume'
   ]

   df_subset = orca.calculate_descriptors(
       smiles_column=df['smiles'],
       descriptors=selected_descriptors
   )

   # Only selected descriptors are calculated
   print(df_subset.columns)
   # Output: Index(['smiles', 'name', 'homo_energy', 'lumo_energy', 'gap_energy', 'dipole_moment', 'molecular_volume'])

Working with Series
-------------------

You can also pass a pandas Series directly::

   # Create a Series with SMILES
   smiles_series = pd.Series(['C1=CC=CC=C1', 'CCO', 'CC(=O)C'])

   # Calculate descriptors
   df_result = orca.calculate_descriptors(smiles_column=smiles_series)

   # Result is a DataFrame with 'smiles' and all descriptor columns
   print(df_result.head())

Progress Tracking
-----------------

The method provides progress information by default::

   # Process many molecules with progress output
   large_df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1'] * 10  # 10 molecules
   })

   df_result = orca.calculate_descriptors(
       smiles_column=large_df['smiles'],
       progress=True  # Default: shows progress
   )

   # Output:
   # INFO - Processing molecule 1/10: C1=CC=CC=C1
   # INFO - Processing molecule 2/10: C1=CC=CC=C1
   # ...

To disable progress output::

   df_result = orca.calculate_descriptors(
       smiles_column=df['smiles'],
       progress=False
   )

Example: QSAR Dataset Preparation
----------------------------------

Here's a complete example of preparing a QSAR dataset::

   import pandas as pd
   from orca_descriptors import Orca

   # Initialize calculator
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

   # Load your molecular dataset
   df = pd.read_csv('molecules.csv')  # Contains 'smiles' column

   # Calculate all descriptors
   df_descriptors = orca.calculate_descriptors(smiles_column=df['smiles'])

   # Filter molecules based on descriptors
   active_molecules = df_descriptors[df_descriptors['gap_energy'] < 3.0]

   # Save results
   df_descriptors.to_csv('molecules_with_descriptors.csv', index=False)

   # Statistical analysis
   print(df_descriptors[['homo_energy', 'lumo_energy', 'gap_energy']].describe())

Example: Custom Descriptor Selection
-------------------------------------

Calculate only energy-related descriptors for a large dataset::

   energy_descriptors = [
       'homo_energy',
       'lumo_energy',
       'gap_energy',
       'total_energy',
       'gibbs_free_energy',
       'entropy',
       'enthalpy'
   ]

   df_energy = orca.calculate_descriptors(
       smiles_column=df['smiles'],
       descriptors=energy_descriptors
   )

   # Analyze energy trends
   print(df_energy.groupby('smiles')[energy_descriptors].mean())

Error Handling
--------------

The method handles errors gracefully. If a descriptor calculation fails for a molecule, it sets that value to ``None`` and continues processing::

   # Some molecules might fail descriptor calculation
   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'INVALID_SMILES', 'CCO']
   })

   df_result = orca.calculate_descriptors(smiles_column=df['smiles'])

   # Failed calculations are marked as None
   print(df_result[df_result['homo_energy'].isna()])

Available Descriptors
---------------------

The following descriptors are available for calculation:

**Energy Descriptors:**
- ``homo_energy`` - HOMO energy (eV)
- ``lumo_energy`` - LUMO energy (eV)
- ``gap_energy`` - HOMO-LUMO gap (eV)
- ``total_energy`` - Total energy (Hartree)
- ``gibbs_free_energy`` - Gibbs free energy (Hartree)
- ``entropy`` - Entropy (J/(mol·K))
- ``enthalpy`` - Enthalpy (Hartree)

**DFT Descriptors:**
- ``ch_potential`` - Chemical potential (eV)
- ``electronegativity`` - Electronegativity (eV)
- ``abs_hardness`` - Absolute hardness (eV)
- ``abs_softness`` - Absolute softness (1/eV)
- ``frontier_electron_density`` - Maximum frontier electron density

**Polar Surface Descriptors:**
- ``dipole_moment`` - Dipole moment (Debye)
- ``polar_surface_area`` - Polar surface area (Å²)
- ``get_min_h_charge`` - Minimum hydrogen charge

**Thermodynamic Descriptors:**
- ``molecular_volume`` - Molecular volume (Å³)

**Multi-dimensional Descriptors:**
- ``num_rotatable_bonds`` - Number of rotatable bonds
- ``wiener_index`` - Wiener index
- ``solvent_accessible_surface_area`` - SASA (Å²)

**Geometric Descriptors:**
- ``xy_shadow`` - XY projection area (Å²)

**Reactivity Descriptors:**
- ``meric`` - Minimum electrophilicity index for carbon (eV)

**Physicochemical Descriptors:**
- ``m_log_p`` - Moriguchi Log P (octanol/water partition coefficient)

**Autocorrelation Descriptors:**
- ``moran_autocorrelation`` - Moran autocorrelation (lag=2, weight='vdw_volume')
- ``autocorrelation_hats`` - HATS autocorrelation (lag=4, unweighted=True)

Tips and Best Practices
------------------------

1. **Caching**: The library automatically caches calculation results. If you recalculate descriptors for the same molecules, it will use cached results.

2. **Batch Processing**: For large datasets, process molecules in batches to manage memory and computation time.

3. **Descriptor Selection**: Use the ``descriptors`` parameter to calculate only the descriptors you need, which can significantly reduce computation time.

4. **Progress Monitoring**: Keep ``progress=True`` (default) to monitor long-running calculations.

5. **Data Validation**: Check for ``None`` values in the result DataFrame to identify molecules that failed descriptor calculation.

Example: Batch Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

   # Process in batches
   batch_size = 10
   all_results = []

   for i in range(0, len(df), batch_size):
       batch = df.iloc[i:i+batch_size]
       batch_result = orca.calculate_descriptors(
           smiles_column=batch['smiles'],
           descriptors=['homo_energy', 'lumo_energy', 'gap_energy']
       )
       all_results.append(batch_result)

   # Combine results
   final_df = pd.concat(all_results, ignore_index=True)

