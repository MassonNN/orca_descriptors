Changelog
=========


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

