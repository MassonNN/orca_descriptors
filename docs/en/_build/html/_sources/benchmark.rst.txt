Benchmark and Time Estimation
==============================

Overview
--------

The library includes a time estimation system that helps predict how long ORCA calculations will take. This is essential for planning large-scale QSAR descriptor calculations.

How It Works
------------

The time estimation system uses a benchmark calculation to calibrate performance on your machine:

1. **Benchmark Calculation**: A single-point calculation is run on benzene (C₁H₆) to measure:
   - Number of basis functions
   - Time per SCF cycle
   - Total calculation time

2. **Scaling Formula**: For new molecules, the system estimates time using:
   - Molecular size scaling (O(N².⁵) for DFT calculations)
   - Method type (SP, Opt, Freq)
   - Number of processors

3. **Parameter Scaling**: The system automatically adjusts for different:
   - Number of processors (accounts for parallel efficiency)
   - Basis sets (scales with size)
   - Functionals (accounts for computational cost)

Why It's Needed
---------------

* **Planning**: Estimate total time for large datasets
* **Resource Management**: Allocate computational resources efficiently
* **Progress Tracking**: Monitor calculation progress
* **Optimization**: Choose optimal parameters for your hardware

Running a Benchmark
-------------------

Using CLI
~~~~~~~~~

.. code-block:: bash

   orca_descriptors run_benchmark --working_dir ./calculations

The benchmark uses benzene as a standard test molecule and saves results to ``.orca_benchmark.json``.

Using Python
~~~~~~~~~~~~

.. code-block:: python

   from orca_descriptors import Orca
   
   orca = Orca(working_dir="./calculations")
   benchmark_data = orca.run_benchmark()
   
   print(f"SCF cycle time: {benchmark_data['scf_time']:.2f} seconds")
   print(f"Number of basis functions: {benchmark_data['n_basis']}")

Estimating Calculation Time
----------------------------

Using CLI
~~~~~~~~~

.. code-block:: bash

   orca_descriptors approximate_time --molecule CCO --method_type Opt

This estimates time without running the actual calculation.

Using Python
~~~~~~~~~~~~

.. code-block:: python

   from orca_descriptors import Orca
   from rdkit.Chem import MolFromSmiles, AddHs
   
   orca = Orca(working_dir="./calculations")
   mol = AddHs(MolFromSmiles("CCO"))
   
   estimated_time = orca.estimate_calculation_time(mol)
   print(f"Estimated time: {estimated_time:.2f} seconds")

Automatic Parameter Scaling
---------------------------

The system automatically scales benchmark data for different parameters. You don't need to re-run the benchmark if you change:

* **Number of processors**: Automatically accounts for parallel efficiency
* **Basis set**: Scales based on basis set size (O(N³.⁵))
* **Functional**: Adjusts for relative computational cost

Example: If your benchmark was run with 1 processor and def2-SVP, you can estimate time for 4 processors and def2-TZVP without re-running the benchmark.

Benchmark File Location
-----------------------

The benchmark data is saved to::

   <working_dir>/.orca_benchmark.json

This file contains:
* Functional and basis set used
* Number of processors
* Number of basis functions
* SCF cycle time
* Total calculation time

You can share this file between different working directories if using the same hardware and ORCA version.

