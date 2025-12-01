Changelog
=========


Version 0.3.0
-------------

Added
~~~~~

* Added new ``ORCABatchProcessing`` class for efficient batch processing of molecular descriptors with pandas compatibility
* Added support for semi-empirical methods (AM1, PM3, PM6, PM7, RM1, MNDO, MNDOD, OM1, OM2, OM3)
* Added ``pre_optimize`` parameter (default: ``True``) for pre-optimizing molecular geometry using MMFF94 force field before ORCA calculations
* Added multiprocessing support for parallel batch processing via ``parallel_mode="multiprocessing"`` parameter
* Added automatic cleanup of all ORCA temporary files (input, output, and all temporary files) after calculations since results are cached
* Added improved error parsing with brief summaries in ``logging.INFO`` and detailed information in ``logging.DEBUG``
* Added ``_pre_optimize_geometry()`` method for MMFF94 geometry optimization
* Added ``_is_semi_empirical()`` method to detect semi-empirical methods

Changed
~~~~~~~

* Refactored batch processing functionality from ``Orca.calculate_descriptors()`` into dedicated ``ORCABatchProcessing`` class
* ``Orca.calculate_descriptors()`` now uses ``ORCABatchProcessing`` internally for backward compatibility
* ``calculate_descriptors()`` now preserves original DataFrame columns (including 'smiles') instead of removing and re-adding them
* Improved molecule hash calculation to include ``pre_optimize`` parameter for proper caching
* Updated molecule hash calculation to exclude ``basis_set`` and ``dispersion_correction`` for semi-empirical methods
* Enhanced input file generation to support semi-empirical methods (no basis set or dispersion correction needed)
* Improved file cleanup to remove all ORCA files (including input and output files) since results are cached

Fixed
~~~~~

* Fixed DataFrame handling in batch processing to preserve all original columns
* Fixed error handling to provide concise error messages in INFO level and detailed information in DEBUG level

Technical Details
~~~~~~~~~~~~~~~~~

* ``ORCABatchProcessing`` supports three parallelization modes: "sequential", "multiprocessing", and "mpirun"
* Pre-optimization uses RDKit's MMFF94 force field for fast geometry optimization before quantum chemical calculations
* All ORCA files are automatically cleaned up after successful calculations, with results stored in cache
* Semi-empirical methods are automatically detected and handled differently from DFT methods
* Batch processing now includes time estimation based on benchmark machine performance


Version 0.2.2
--------------

Added
~~~~~

* Added ``numpy>=1.20.0`` to project dependencies (numpy was used but not declared)
* Added dynamic time estimation updates in batch processing - time estimates are now refined based on actual execution times of previous molecules
* Added ``_get_available_descriptors()`` method to dynamically discover available descriptor methods

Changed
~~~~~~~

* Updated dipole moment parser to prioritize gas-phase values when available (for calculations without solvation)
* Improved time estimation algorithm:
  * Changed scaling exponent from O(N^3.5) to O(N^2.5) for more realistic estimates
  * Uses ``total_time`` from benchmark instead of ``scf_time`` as base unit
  * More realistic optimization step estimation (15-35 steps instead of 10-50)
  * Removed artificial 24-hour time cap
* Refactored ``calculate_descriptors()`` method:
  * Removed redundant code duplication (replaced large if-elif chain with ``getattr``-based method calls)
  * Removed redundant ``all_descriptors`` list - descriptors are now discovered dynamically
  * Removed unnecessary comments
  * Improved code maintainability and readability

Fixed
~~~~~

* Fixed dipole moment parser to correctly extract gas-phase values from ORCA output when available
* Fixed time estimation showing unrealistic values (e.g., 47 hours for 2 molecules) - now provides accurate estimates based on actual benchmark data

Technical Details
~~~~~~~~~~~~~~~~~

* Time estimator now uses exponential moving average for better prediction accuracy
* Descriptor methods are called dynamically using ``getattr(self, desc_name)``
* Automatic descriptor discovery eliminates need to maintain manual descriptor lists

