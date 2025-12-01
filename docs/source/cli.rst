Command-Line Interface
======================

The library provides a command-line interface (CLI) for running benchmarks and estimating calculation times.

Commands
--------

run_benchmark
~~~~~~~~~~~~~

Run a benchmark calculation to calibrate time estimation::

   orca_descriptors run_benchmark [OPTIONS]

Options:

* ``--working_dir``: Working directory for calculations (default: current directory)
* ``--functional``: DFT functional (default: PBE0)
* ``--basis_set``: Basis set (default: def2-SVP)
* ``--n_processors``: Number of processors (default: 1)
* All other ORCA parameters are also available

Example::

   orca_descriptors run_benchmark --working_dir ./calculations --n_processors 4

approximate_time
~~~~~~~~~~~~~~~~

Estimate calculation time for a molecule without running the calculation::

   orca_descriptors approximate_time --molecule SMILES [OPTIONS]

Required arguments:

* ``--molecule``: SMILES string of the molecule

Options:

* ``--method_type``: Calculation type - Opt, SP, or Freq (default: Opt)
* ``--n_opt_steps``: Expected number of optimization steps (for Opt method)
* All ORCA parameters are available

Example::

   orca_descriptors approximate_time \\
       --molecule CCO \\
       --method_type Opt \\
       --functional PBE0 \\
       --basis_set def2-TZVP \\
       --n_processors 8

Available Parameters
-------------------

All parameters from the ``Orca`` class are available as command-line arguments:

* ``--script_path``: Path to ORCA executable (default: 'orca')
* ``--working_dir``: Working directory for calculations
* ``--output_dir``: Directory for output files
* ``--functional``: DFT functional (default: PBE0)
* ``--basis_set``: Basis set (default: def2-SVP)
* ``--method_type``: Calculation type: Opt, SP, or Freq (default: Opt)
* ``--dispersion_correction``: Dispersion correction, e.g., D3BJ (default: D3BJ). Use 'None' to disable.
* ``--solvation_model``: Solvation model, e.g., 'COSMO(Water)' (default: None)
* ``--n_processors``: Number of processors (default: 1)
* ``--max_scf_cycles``: Maximum SCF cycles (default: 100)
* ``--scf_convergence``: SCF convergence threshold (default: 1e-6)
* ``--charge``: Molecular charge (default: 0)
* ``--multiplicity``: Spin multiplicity (default: 1)
* ``--cache_dir``: Directory for caching results
* ``--log_level``: Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)
* ``--max_wait``: Maximum time to wait for output file creation in seconds (default: 300)

Help
----

Get help for any command::

   orca_descriptors --help
   orca_descriptors run_benchmark --help
   orca_descriptors approximate_time --help

