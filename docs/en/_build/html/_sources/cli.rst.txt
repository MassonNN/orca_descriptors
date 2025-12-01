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
* ``--functional``: DFT functional or semi-empirical method (default: AM1)
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
       --functional AM1 \\
       --basis_set def2-TZVP \\
       --n_processors 8

clear
~~~~~

Remove all ORCA files in the working directory (useful for cleaning up files that weren't removed due to errors)::

   orca_descriptors clear [OPTIONS]

This command removes all ORCA-related files (`.inp`, `.out`, `.log`, `.gbw`, `.cube`, `.prop`, etc.) from the working directory.

Options:

* ``--working_dir``: Working directory to clean (default: current directory)
* All other ORCA parameters are available but not used

Example::

   orca_descriptors clear --working_dir ./calculations

purge_cache
~~~~~~~~~~~

Remove ORCA cache::

   orca_descriptors purge_cache [OPTIONS]

This command clears the ORCA cache directory, removing all cached calculation results.

Options:

* ``--cache_dir``: Cache directory to purge (default: output_dir/.orca_cache)
* ``--output_dir``: Output directory (used to determine cache location if cache_dir not specified)
* All other ORCA parameters are available but not used

Example::

   orca_descriptors purge_cache --output_dir ./calculations

Available Parameters
-------------------

All parameters from the ``Orca`` class are available as command-line arguments:

* ``--script_path``: Path to ORCA executable (default: 'orca')
* ``--working_dir``: Working directory for calculations
* ``--output_dir``: Directory for output files
* ``--functional``: DFT functional or semi-empirical method (default: AM1)
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
   orca_descriptors clear --help
   orca_descriptors purge_cache --help

