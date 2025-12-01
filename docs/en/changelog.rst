Changelog
=========


Version 0.3.4
-------------

Added
~~~~~

* Added ``clear`` CLI command to remove all ORCA files in the working directory (useful for cleaning up files that weren't removed due to errors)
* Added ``purge_cache`` CLI command to remove ORCA cache

Changed
~~~~~~~

* Improved cache system to preserve original file extensions (`.out`, `.log`, `.smd.out`) when storing files
* Enhanced NBO stabilization energy parsing with better error messages indicating NBOEXE environment variable requirement
* Updated test for AM1 Mayer bond indices to account for different values compared to DFT methods

Removed
~~~~~~~

* Removed ESP extrema descriptor (``get_esp_extrema``) - not implemented due to limitations in ORCA 6.0.1 for generating ESP cube files directly (would require orca_plot utility integration)

Fixed
~~~~~

* Fixed cache system issue where files were saved with `.out` extension regardless of original extension, causing file not found errors
* Fixed NBO stabilization energy parsing to properly detect when NBOEXE environment variable is not set
* Fixed NMR chemical shifts parsing to work correctly with cached files
* Fixed test for AM1 Mayer bond indices to accept different value ranges compared to DFT

Technical Details
~~~~~~~~~~~~~~~~~

* Cache now preserves original file extensions to ensure correct file retrieval
* NBO analysis requires NBOEXE environment variable to be set to point to NBO executable (nbo6.exe or nbo5.exe)
* ESP extrema calculation would require integration with orca_plot utility, which is not currently implemented
* CLI commands ``clear`` and ``purge_cache`` help maintain clean working directories and cache management


Version 0.3.3
-------------

Changed
~~~~~~~

* Major code refactoring: split large ``orca.py`` file (1620 lines) into modular structure:
  * Created ``base.py`` with ``OrcaBase`` class containing common utility methods
  * Created ``calculation.py`` with ``CalculationMixin`` for calculation execution methods
  * Created ``decorators.py`` with ``handle_x_molecule`` decorator
  * Split descriptors into separate modules by category:
    * ``descriptors/electronic.py`` - Electronic property descriptors
    * ``descriptors/energy.py`` - Energy-related descriptors
    * ``descriptors/structural.py`` - Structural property descriptors
    * ``descriptors/topological.py`` - Topological descriptors
    * ``descriptors/misc.py`` - Miscellaneous descriptors
  * Main ``Orca`` class now uses multiple inheritance from mixins
  * Improved code organization and maintainability
  * Removed redundant comments throughout the codebase

Technical Details
~~~~~~~~~~~~~~~~~

* The new modular structure makes it easier to extend functionality and maintain code
* Descriptors are organized by category for better code navigation
* All functionality remains backward compatible


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

