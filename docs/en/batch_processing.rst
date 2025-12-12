Batch Processing Tutorial
===========================

This tutorial demonstrates how to use the ORCA Descriptors library for efficient batch processing of molecular descriptors with support for pandas integration, parallelization, and advanced descriptor definition.

Installation
-------------

To use batch processing with pandas integration, install pandas as an optional dependency::

   pip install 'orca-descriptors[pandas]'

Or install pandas separately::

   pip install pandas

Basic Usage
-----------

The ``ORCABatchProcessing`` class provides efficient batch processing of molecular descriptors.

Creating a Batch Processor
~~~~~~~~~~~~~~~~~~~~~~~~~~

Initialize a batch processor with your ORCA configuration::

   from orca_descriptors import Orca, ORCABatchProcessing

   # Create an Orca instance
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

   # Create a batch processor
   batch_processing = ORCABatchProcessing(
       orca=orca,
       working_dir=".",
   )

You can also create a batch processor without an existing Orca instance::

   batch_processing = ORCABatchProcessing(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )

Calculating Descriptors
~~~~~~~~~~~~~~~~~~~~~~~

Calculate descriptors for a list of SMILES strings::

   smiles_list = ["C1=CC=CC=C1", "CCO", "CC(=O)C"]
   
   # Calculate all available descriptors
   result = batch_processing.calculate_descriptors(smiles_list)
   
   # Result is a DataFrame (if pandas available) or list of dictionaries
   print(result)

Working with Pandas
-------------------

The batch processor seamlessly integrates with pandas DataFrames and Series.

DataFrame Input
~~~~~~~~~~~~~~~

Pass a DataFrame with a 'smiles' column::

   import pandas as pd
   
   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'CCO', 'CC(=O)C'],
       'name': ['Benzene', 'Ethanol', 'Acetone']
   })
   
   # Calculate descriptors - original columns are preserved
   df_result = batch_processing.calculate_descriptors(df['smiles'])
   
   # DataFrame contains original columns + descriptor columns
   print(df_result.columns)
   # Output: Index(['smiles', 'name', 'homo_energy', 'lumo_energy', ...])

Series Input
~~~~~~~~~~~~

You can also pass a pandas Series directly::

   smiles_series = pd.Series(['C1=CC=CC=C1', 'CCO', 'CC(=O)C'])
   
   df_result = batch_processing.calculate_descriptors(smiles_series)
   
   # Result is a DataFrame with descriptor columns only
   print(df_result.head())

List Input
~~~~~~~~~~

Plain Python lists are also supported::

   smiles_list = ['C1=CC=CC=C1', 'CCO', 'CC(=O)C']
   
   result = batch_processing.calculate_descriptors(smiles_list)
   
   # Result is a DataFrame (if pandas available) or list of dictionaries

Defining Descriptors with XMolecule API
----------------------------------------

The XMolecule API allows you to define descriptors with their parameters using a special placeholder molecule. This is especially useful for descriptors that require parameters.

Creating an X Molecule
~~~~~~~~~~~~~~~~~~~~~~~

Create an X molecule using the ``x_molecule()`` method::

   x = batch_processing.x_molecule()

Using X Molecule to Define Descriptors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Call descriptor methods on the Orca instance with the X molecule to define descriptors with parameters::

   descriptors = [
       orca.ch_potential(x),
       orca.electronegativity(x),
       orca.abs_hardness(x),
       orca.topological_distance(x, 'O', 'O'),  # Distance between oxygen atoms
       orca.mo_energy(x, -3),  # HOMO-2 energy (index -3)
   ]
   
   result = batch_processing.calculate_descriptors(
       smiles_list,
       descriptors=descriptors
   )
   
   # Result columns include parameterized descriptor names:
   # 'ch_potential', 'electronegativity', 'abs_hardness', 
   # 'topological_distance_O_O', 'mo_energy_-3'

Complete Example
~~~~~~~~~~~~~~~~

Here's a complete example using the XMolecule API::

   from orca_descriptors import Orca, ORCABatchProcessing
   import pandas as pd
   
   # Initialize
   orca = Orca(functional="PBE0", basis_set="def2-SVP")
   batch_processing = ORCABatchProcessing(orca=orca)
   
   # Create X molecule
   x = batch_processing.x_molecule()
   
   # Define descriptors with parameters
   descriptors = [
       orca.homo_energy(x),
       orca.lumo_energy(x),
       orca.gap_energy(x),
       orca.mo_energy(x, -1),  # HOMO
       orca.mo_energy(x, -2),  # HOMO-1
       orca.topological_distance(x, 'C', 'C'),  # C-C distances
       orca.topological_distance(x, 'O', 'O'),  # O-O distances
   ]
   
   # Load your dataset
   df = pd.read_csv('molecules.csv')
   
   # Calculate descriptors
   df_result = batch_processing.calculate_descriptors(
       df['smiles'],
       descriptors=descriptors
   )
   
   # Save results
   df_result.to_csv('molecules_with_descriptors.csv', index=False)

Selecting Descriptors
---------------------

You can specify which descriptors to calculate in two ways:

Method 1: Using Descriptor Names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass a list of descriptor names as strings::

   selected_descriptors = [
       'homo_energy',
       'lumo_energy',
       'gap_energy',
       'dipole_moment',
       'molecular_volume'
   ]
   
   result = batch_processing.calculate_descriptors(
       smiles_list,
       descriptors=selected_descriptors
   )

Method 2: Using XMolecule API (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the XMolecule API for descriptors with parameters::

   x = batch_processing.x_molecule()
   
   descriptors = [
       orca.homo_energy(x),
       orca.mo_energy(x, -3),  # With parameter
       orca.topological_distance(x, 'O', 'O'),  # With parameters
   ]
   
   result = batch_processing.calculate_descriptors(
       smiles_list,
       descriptors=descriptors
   )

Parallelization
---------------

The batch processor supports multiple parallelization modes for efficient processing of large datasets.

Sequential Processing
~~~~~~~~~~~~~~~~~~~~~

Default mode - processes molecules one by one::

   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="sequential"
   )

Multiprocessing
~~~~~~~~~~~~~~~

Use Python multiprocessing to run multiple ORCA calculations in parallel::

   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="multiprocessing",
       n_workers=4  # Number of parallel workers
   )
   
   # Process molecules in parallel
   result = batch_processing.calculate_descriptors(smiles_list)

The multiprocessing mode automatically adjusts time estimates based on the number of workers and parallel efficiency.

MPI (mpirun)
~~~~~~~~~~~~

Use ORCA's built-in MPI parallelization::

   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="mpirun",
       use_mpirun=True,
       n_processors=4  # Processors per ORCA calculation
   )

Progress Tracking
-----------------

The batch processor provides detailed progress information.

Progress Output
~~~~~~~~~~~~~~~

By default, progress information is displayed::

   result = batch_processing.calculate_descriptors(
       smiles_list,
       progress=True  # Default
   )
   
   # Output:
   # INFO - Estimating calculation times...
   # INFO - Processing molecule 1/10 (remaining: 10, estimated time: ~5m)
   # INFO - Processing molecule 2/10 (remaining: 9, estimated time: ~4m 30s, avg: 30.5s/molecule)
   # INFO - Processing molecule 3/10 (remaining: 8, CACHED, avg: 0.1s/molecule)

Cached Molecules
~~~~~~~~~~~~~~~~

When a molecule is found in cache, the progress shows "CACHED" instead of time estimates::

   # First calculation
   result1 = batch_processing.calculate_descriptors(smiles_list)
   
   # Second calculation (from cache)
   result2 = batch_processing.calculate_descriptors(smiles_list)
   # Output: INFO - Processing molecule 1/10 (remaining: 10, CACHED)

Disable Progress
~~~~~~~~~~~~~~~~

To disable progress output::

   result = batch_processing.calculate_descriptors(
       smiles_list,
       progress=False
   )

Error Handling
--------------

The batch processor handles errors gracefully. If a descriptor calculation fails for a molecule, it sets that value to ``None`` and continues processing::

   df = pd.DataFrame({
       'smiles': ['C1=CC=CC=C1', 'INVALID_SMILES', 'CCO']
   })
   
   result = batch_processing.calculate_descriptors(df['smiles'])
   
   # Failed calculations are marked as None
   print(result[result['homo_energy'].isna()])

Error messages are logged at different levels:

- Brief error summary: ``logging.INFO``
- Detailed error information: ``logging.DEBUG``

Example: QSAR Dataset Preparation
----------------------------------

Here's a complete example of preparing a QSAR dataset::

   import pandas as pd
   from orca_descriptors import Orca, ORCABatchProcessing
   
   # Initialize calculator
   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       method_type="Opt",
       n_processors=4,
   )
   
   # Create batch processor with multiprocessing
   batch_processing = ORCABatchProcessing(
       orca=orca,
       parallel_mode="multiprocessing",
       n_workers=4,
   )
   
   # Create X molecule for descriptor definition
   x = batch_processing.x_molecule()
   
   # Define descriptors
   descriptors = [
       orca.homo_energy(x),
       orca.lumo_energy(x),
       orca.gap_energy(x),
       orca.ch_potential(x),
       orca.electronegativity(x),
       orca.abs_hardness(x),
       orca.dipole_moment(x),
       orca.molecular_volume(x),
   ]
   
   # Load your molecular dataset
   df = pd.read_csv('molecules.csv')  # Contains 'smiles' column
   
   # Calculate descriptors
   df_descriptors = batch_processing.calculate_descriptors(
       smiles_column=df['smiles'],
       descriptors=descriptors
   )
   
   # Filter molecules based on descriptors
   active_molecules = df_descriptors[df_descriptors['gap_energy'] < 3.0]
   
   # Save results
   df_descriptors.to_csv('molecules_with_descriptors.csv', index=False)
   
   # Statistical analysis
   print(df_descriptors[['homo_energy', 'lumo_energy', 'gap_energy']].describe())

Example: Batch Processing with Custom Parameters
-------------------------------------------------

Process molecules in batches with different configurations::

   # Process in batches
   batch_size = 10
   all_results = []
   
   for i in range(0, len(df), batch_size):
       batch = df.iloc[i:i+batch_size]
       
       x = batch_processing.x_molecule()
       descriptors = [
           orca.homo_energy(x),
           orca.lumo_energy(x),
           orca.gap_energy(x)
       ]
       
       batch_result = batch_processing.calculate_descriptors(
           smiles_column=batch['smiles'],
           descriptors=descriptors
       )
       all_results.append(batch_result)
   
   # Combine results
   final_df = pd.concat(all_results, ignore_index=True)

Remote Cache
-------------

The batch processor supports remote caching via API, allowing you to share calculation results across different machines and users. A public cache server is available at ``https://api.orca-descriptors.massonnn.ru``.

To use remote cache:

1. Register at `orca-descriptors.massonnn.ru <https://orca-descriptors.massonnn.ru>`_ and generate an API token
2. Provide the API token when initializing the Orca instance::

   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       cache_api_token="your_api_token_here",
   )
   
   batch_processing = ORCABatchProcessing(orca=orca)

The remote cache works transparently in batch processing:
- All molecules are checked for cache (both local and remote) before starting calculations
- Cached molecules are processed immediately and excluded from further calculations
- New calculations are automatically uploaded to the remote server (if you have upload permissions)
- Progress messages show when molecules are found in cache

Cache-Only Mode
---------------

You can enable cache-only mode to use only cached results without running ORCA calculations::

   orca = Orca(
       functional="PBE0",
       basis_set="def2-SVP",
       cache_api_token="your_api_token_here",
       cache_only=True,  # Only use cache, don't run calculations
   )
   
   batch_processing = ORCABatchProcessing(orca=orca)
   
   result = batch_processing.calculate_descriptors(smiles_list)

In cache-only mode:
- Only cached results (local or remote) are used
- Molecules not found in cache return ``None`` for descriptors
- No ORCA calculations are performed
- Progress messages show warnings for molecules not found in cache
- Useful for quickly retrieving results from cache without running expensive calculations

Tips and Best Practices
-----------------------

1. **Caching**: The library automatically caches calculation results. Recalculating descriptors for the same molecules uses cached results, significantly speeding up subsequent runs.

2. **Remote Cache**: Use remote cache to share results across machines and collaborate with others. Register at `orca-descriptors.massonnn.ru <https://orca-descriptors.massonnn.ru>`_ to get an API token.

3. **Multiprocessing**: For large datasets, use ``parallel_mode="multiprocessing"`` with an appropriate number of workers (typically equal to CPU cores).

4. **Descriptor Selection**: Use the XMolecule API to define descriptors with parameters, making your code more readable and maintainable.

5. **Progress Monitoring**: Keep ``progress=True`` (default) to monitor long-running calculations and identify cached molecules.

6. **Data Validation**: Check for ``None`` values in the result DataFrame to identify molecules that failed descriptor calculation or were not found in cache (when using cache-only mode).

7. **Memory Management**: For very large datasets, process molecules in batches to manage memory usage.

Available Descriptors
---------------------

See the :doc:`descriptors` documentation for a complete list of available descriptors and their parameters.

