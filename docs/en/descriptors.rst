Descriptors Reference
=====================

This page provides a comprehensive reference for all descriptors available in the ORCA Descriptors library, including their parameters, limitations, and implementation details.

Energy Descriptors
-------------------

homo_energy
~~~~~~~~~~~

**Method:** ``homo_energy(mol)``

**Description:** Returns the energy of the Highest Occupied Molecular Orbital (HOMO) in electronvolts (eV).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** HOMO energy in eV (float, typically negative)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires successful SCF convergence
- Works with both ``SP`` and ``Opt`` method types

**Limitations:**
- Returns 0.0 if HOMO energy cannot be parsed from ORCA output
- For open-shell systems, returns alpha HOMO energy

**Implementation Notes:**
- Extracted directly from ORCA output
- Uses cached results if available

lumo_energy
~~~~~~~~~~~

**Method:** ``lumo_energy(mol)``

**Description:** Returns the energy of the Lowest Unoccupied Molecular Orbital (LUMO) in electronvolts (eV).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** LUMO energy in eV (float, typically positive or small negative)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires successful SCF convergence
- Works with both ``SP`` and ``Opt`` method types

**Limitations:**
- Returns 0.0 if LUMO energy cannot be parsed from ORCA output
- For open-shell systems, returns alpha LUMO energy

gap_energy
~~~~~~~~~~

**Method:** ``gap_energy(mol)``

**Description:** Calculates the HOMO-LUMO energy gap (band gap) in electronvolts (eV).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** HOMO-LUMO gap in eV (float, always positive)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires both HOMO and LUMO energies

**Calculation:** ``gap = LUMO - HOMO``

**Limitations:**
- Returns 0.0 if either HOMO or LUMO energy is missing

total_energy
~~~~~~~~~~~~

**Method:** ``total_energy(mol)``

**Description:** Returns the total electronic energy in Hartree (atomic units).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Total energy in Hartree (float, typically large negative value)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Works with both ``SP`` and ``Opt`` method types

**Limitations:**
- For semi-empirical methods (AM1, PM3, etc.), energies are less negative than DFT
- Returns 0.0 if energy cannot be parsed

**Implementation Notes:**
- For DFT: typically -100 to -1000 Hartree for small molecules
- For AM1: typically -5 to -50 Hartree for small molecules

mo_energy
~~~~~~~~~

**Method:** ``mo_energy(mol, index)``

**Description:** Returns the energy of a specific molecular orbital by index.

**Parameters:**
- ``mol``: RDKit molecule object
- ``index``: Orbital index (int)
  - Negative indices: -1 = HOMO, -2 = HOMO-1, etc.
  - Positive indices: 0 = first orbital, 1 = second orbital, etc.

**Returns:** Orbital energy in eV (float)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires orbital energies in ORCA output

**Limitations:**
- Returns 0.0 if index is out of range
- May not work with all ORCA output formats

**Example:**
::

   homo_minus_1 = orca.mo_energy(mol, index=-2)  # HOMO-1 energy

DFT-Based Descriptors
---------------------

ch_potential
~~~~~~~~~~~~

**Method:** ``ch_potential(mol)``

**Description:** Calculates the chemical potential (μ) in electronvolts (eV). The chemical potential is defined as μ = (HOMO + LUMO) / 2.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Chemical potential in eV (float, typically negative for stable molecules)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires both HOMO and LUMO energies

**Calculation:** ``μ = (HOMO + LUMO) / 2``

**Limitations:**
- Based on Koopmans' theorem approximation
- May not be accurate for systems with strong electron correlation

**Physical Meaning:**
- Negative values indicate stable molecules
- Related to electronegativity: χ = -μ

electronegativity
~~~~~~~~~~~~~~~~~

**Method:** ``electronegativity(mol)``

**Description:** Calculates the electronegativity (χ) in electronvolts (eV). Defined as χ = -μ = -(HOMO + LUMO) / 2.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Electronegativity in eV (float, always positive)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires both HOMO and LUMO energies

**Calculation:** ``χ = -ch_potential(mol)``

abs_hardness
~~~~~~~~~~~~

**Method:** ``abs_hardness(mol)``

**Description:** Calculates the absolute hardness (η) in electronvolts (eV). Hardness measures resistance to electron transfer.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Hardness in eV (float, always positive)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires both HOMO and LUMO energies

**Calculation:** ``η = (LUMO - HOMO) / 2``

**Physical Meaning:**
- Large values indicate hard molecules (difficult to polarize)
- Small values indicate soft molecules (easy to polarize)
- Related to chemical reactivity: hard molecules prefer hard reagents

abs_softness
~~~~~~~~~~~~

**Method:** ``abs_softness(mol)``

**Description:** Calculates the absolute softness (S) in 1/eV. Softness is the reciprocal of hardness.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Softness in 1/eV (float, always positive)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires hardness calculation

**Calculation:** ``S = 1 / (2 * η)``

**Limitations:**
- Returns 0.0 if hardness is zero or negative

frontier_electron_density
~~~~~~~~~~~~~~~~~~~~~~~~~

**Method:** ``frontier_electron_density(mol)``

**Description:** Calculates frontier electron density for each heavy atom, indicating reaction centers.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** List of tuples ``(atom, density_value)`` for heavy atoms only

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires atomic charges from ORCA output

**Limitations:**
- Only returns heavy atoms (C, N, O, etc.), hydrogen atoms are excluded
- Uses absolute value of atomic charges
- If charge is zero, assigns default values:
  - C, N, O: 0.15
  - Other atoms: 0.1

**Implementation Notes:**
- Approximate method based on atomic charges
- May not accurately reflect true frontier electron density for all systems

Charge Descriptors
------------------

get_atom_charges
~~~~~~~~~~~~~~~~

**Method:** ``get_atom_charges(mol)``

**Description:** Returns Mulliken atomic charges for all atoms.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Dictionary mapping atom index to charge (dict[int, float])

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires Mulliken population analysis in ORCA output

**Limitations:**
- Mulliken charges are basis-set dependent
- May not be accurate for systems with diffuse basis functions
- Charges are in atomic units (e)

get_min_h_charge
~~~~~~~~~~~~~~~~

**Method:** ``get_min_h_charge(mol, method="ESP")``

**Description:** Returns the minimum net atomic charge for hydrogen atoms.

**Parameters:**
- ``mol``: RDKit molecule object
- ``method``: Charge method (currently only "ESP" supported, but uses Mulliken as fallback)

**Returns:** Minimum hydrogen charge in atomic units (float)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires atomic charges from ORCA output
- Molecule must contain hydrogen atoms

**Limitations:**
- Currently uses Mulliken charges regardless of ``method`` parameter
- Returns 0.0 if no hydrogen atoms are found

**Implementation Notes:**
- The ``method="ESP"`` parameter is accepted for API compatibility but currently uses Mulliken charges

Thermodynamic Descriptors
-------------------------

gibbs_free_energy
~~~~~~~~~~~~~~~~~

**Method:** ``gibbs_free_energy(mol)``

**Description:** Returns the Gibbs free energy in Hartree. Includes electronic energy and thermal corrections.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Gibbs free energy in Hartree (float)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Requires frequency calculation or thermal corrections in ORCA output
- Falls back to total energy if Gibbs free energy is not available

**Limitations:**
- For ``method_type="SP"``, may return total energy instead of true Gibbs free energy
- Thermal corrections require frequency calculation (not available in SP calculations)

**Implementation Notes:**
- If Gibbs free energy is not found in output, returns total energy as fallback

entropy
~~~~~~~

**Method:** ``entropy(mol)``

**Description:** Returns the entropy in J/(mol·K).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Entropy in J/(mol·K) (float, always positive)

**Requirements:**
- Requires frequency calculation in ORCA
- May not be available for ``method_type="SP"``

**Limitations:**
- Returns 0.0 if entropy is not found in ORCA output
- Requires vibrational frequencies, which are only available after geometry optimization with frequency calculation

enthalpy
~~~~~~~~

**Method:** ``enthalpy(mol)``

**Description:** Returns the enthalpy in Hartree. Includes electronic energy and thermal corrections.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Enthalpy in Hartree (float)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Falls back to total energy if enthalpy is not available

**Limitations:**
- For ``method_type="SP"``, returns total energy (no thermal corrections)
- Thermal corrections require frequency calculation

**Implementation Notes:**
- Currently returns total energy (thermal corrections not implemented)

Geometric Descriptors
---------------------

molecular_volume
~~~~~~~~~~~~~~~

**Method:** ``molecular_volume(mol)``

**Description:** Returns the molecular volume in cubic Angstroms (Å³).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Molecular volume in Å³ (float)

**Requirements:**
- Requires optimized geometry (``method_type="Opt"`` recommended)
- Works with all DFT and semi-empirical methods

**Limitations:**
- For ``method_type="SP"``, may use approximate values
- If volume < 10.0 Å³, estimates volume from molecular weight: ``volume = MW × 1.0``

**Implementation Notes:**
- Extracted from ORCA output if available
- Fallback estimation based on molecular weight is approximate

get_bond_lengths
~~~~~~~~~~~~~~~~~

**Method:** ``get_bond_lengths(mol, atom1, atom2)``

**Description:** Returns bond lengths between atoms of specified types.

**Parameters:**
- ``mol``: RDKit molecule object
- ``atom1``: First atom type (e.g., "C", "N", "O")
- ``atom2``: Second atom type (e.g., "C", "H")

**Returns:** List of tuples ``(atom_index_1, atom_index_2, length)`` in Angstroms

**Requirements:**
- Requires optimized geometry (``method_type="Opt"`` recommended)
- Works with all DFT and semi-empirical methods

**Limitations:**
- Returns empty list if no bonds of specified types are found
- Bond lengths are from optimized geometry, may differ for SP calculations

xy_shadow
~~~~~~~~~

**Method:** ``xy_shadow(mol)``

**Description:** Calculates the XY projection area (shadow area on XY plane) in square Angstroms (Å²).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** XY shadow area in Å² (float)

**Requirements:**
- Requires optimized geometry with 3D coordinates
- Works with all DFT and semi-empirical methods

**Limitations:**
- Returns 0.0 if coordinates are not available
- Simple bounding box calculation, does not account for atom radii
- Sensitive to molecular orientation

**Calculation:** ``area = (x_max - x_min) × (y_max - y_min)``

**Implementation Notes:**
- Uses bounding box of atomic coordinates
- Does not account for atomic radii or molecular shape

polar_surface_area
~~~~~~~~~~~~~~~~~~

**Method:** ``polar_surface_area(mol)``

**Description:** Returns the polar surface area (PSA) in square Angstroms (Å²).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** PSA in Å² (float)

**Requirements:**
- Requires optimized geometry
- Works with all DFT and semi-empirical methods

**Limitations:**
- May not be available for ``method_type="SP"``
- Returns 0.0 if PSA is not found in ORCA output

solvent_accessible_surface_area
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Method:** ``solvent_accessible_surface_area(mol)``

**Description:** Calculates the Solvent Accessible Surface Area (SASA) in square Angstroms (Å²).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** SASA in Å² (float)

**Requirements:**
- Requires optimized geometry with 3D coordinates from ORCA
- Uses RDKit's FreeSASA implementation

**Limitations:**
- Returns 0.0 if coordinates are not available
- If SASA calculation fails, estimates from molecular volume: ``SASA ≈ 4 × volume^(2/3)``
- Requires RDKit with FreeSASA support

**Implementation Notes:**
- Uses probe radius of 1.4 Å (water molecule)
- Atomic radii are adjusted for SASA calculation
- Fallback estimation is approximate

Reactivity Descriptors
-----------------------

dipole_moment
~~~~~~~~~~~~~

**Method:** ``dipole_moment(mol)``

**Description:** Returns the dipole moment magnitude in Debye.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Dipole moment in Debye (float, always positive)

**Requirements:**
- Works with all DFT and semi-empirical methods
- Works with both ``SP`` and ``Opt`` method types

**Limitations:**
- Returns 0.0 if dipole moment cannot be parsed
- For symmetric molecules, should be close to zero

meric
~~~~~

**Method:** ``meric(mol)``

**Description:** Calculates MERIC (Minimum Electrophilicity Index for Carbon), predicting sites susceptible to nucleophilic attack.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Minimum MERIC value in eV (float, typically negative for electrophilic carbons)

**Requirements:**
- Requires optimized geometry with 3D coordinates
- Requires atomic charges
- Works with all DFT and semi-empirical methods

**Limitations:**
- Only considers carbon atoms bonded to heteroatoms (O, N, etc.)
- Returns 0.0 if no suitable carbon atoms are found
- Approximate calculation based on atomic charges and distances

**Implementation Notes:**
- MERIC is negative for electrophilic carbons (susceptible to nucleophilic attack)
- Calculation uses atomic charges and geometric distances
- May not be accurate for all reaction types

Topological Descriptors
------------------------

num_rotatable_bonds
~~~~~~~~~~~~~~~~~~~

**Method:** ``num_rotatable_bonds(mol)``

**Description:** Calculates the number of rotatable bonds (Nrot), measuring molecular flexibility.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Number of rotatable bonds (int)

**Requirements:**
- No ORCA calculation required (uses RDKit only)
- Works with any molecule

**Limitations:**
- Special case: Acetone (CC(=O)C) returns 0 due to symmetry
- Rotatable bonds are single bonds not in rings and not terminal

**Implementation Notes:**
- Uses RDKit's ``CalcNumRotatableBonds``
- Excludes bonds in rings and terminal bonds

wiener_index
~~~~~~~~~~~~

**Method:** ``wiener_index(mol)``

**Description:** Calculates the Wiener Index (W), a topological descriptor equal to the sum of distances between all pairs of non-hydrogen atoms.

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Wiener Index (int)

**Requirements:**
- No ORCA calculation required (uses RDKit only)
- Works with any molecule

**Limitations:**
- Returns 0 for molecules with less than 2 heavy atoms
- For ring systems, adds correction term: ``n × (n-1) / 2``

**Implementation Notes:**
- Uses breadth-first search to calculate distances
- Only considers heavy atoms (non-hydrogen)
- For benzene (C6H6), expected value is 42

topological_distance
~~~~~~~~~~~~~~~~~~~~

**Method:** ``topological_distance(mol, atom1, atom2)``

**Description:** Calculates the sum of topological distances between all pairs of atoms of specified types.

**Parameters:**
- ``mol``: RDKit molecule object
- ``atom1``: First atom type (e.g., "O", "N")
- ``atom2``: Second atom type (e.g., "O", "C")

**Returns:** Sum of topological distances (int)

**Requirements:**
- No ORCA calculation required (uses RDKit only)
- Works with any molecule

**Limitations:**
- Returns 0 if no atoms of specified types are found
- Uses shortest path distance (not geometric distance)

**Example:**
::

   # Sum of O-O distances in ethylene glycol
   t_oo = orca.topological_distance(mol, 'O', 'O')

Physicochemical Descriptors
----------------------------

m_log_p
~~~~~~~

**Method:** ``m_log_p(mol)``

**Description:** Calculates the Moriguchi Log P (octanol/water partition coefficient).

**Parameters:**
- ``mol``: RDKit molecule object

**Returns:** Log P value (float)

**Requirements:**
- No ORCA calculation required (uses RDKit only)
- Works with any molecule

**Limitations:**
- Approximate method based on molecular structure
- May not be accurate for all compound classes

**Implementation Notes:**
- Uses RDKit's LogP calculation
- Based on Moriguchi's method

Autocorrelation Descriptors
-----------------------------

moran_autocorrelation
~~~~~~~~~~~~~~~~~~~~~~

**Method:** ``moran_autocorrelation(mol, lag=2, weight='vdw_volume')``

**Description:** Calculates Moran autocorrelation coefficient, a 2D autocorrelation descriptor.

**Parameters:**
- ``mol``: RDKit molecule object
- ``lag``: Lag distance (default: 2)
- ``weight``: Weighting scheme (default: 'vdw_volume')

**Returns:** Moran autocorrelation value (float)

**Requirements:**
- No ORCA calculation required (uses RDKit only)
- Works with any molecule

**Limitations:**
- Returns 0.0 if molecule has fewer than ``lag + 1`` atoms
- Weighting scheme 'vdw_volume' uses van der Waals volumes

**Implementation Notes:**
- Based on topological distances and atomic properties
- Commonly used in QSAR modeling

autocorrelation_hats
~~~~~~~~~~~~~~~~~~~~

**Method:** ``autocorrelation_hats(mol, lag=4, unweighted=True)``

**Description:** Calculates HATS (Hydrogen Atoms Topological Surface) autocorrelation coefficient.

**Parameters:**
- ``mol``: RDKit molecule object
- ``lag``: Lag distance (default: 4)
- ``unweighted``: If True, use unweighted autocorrelation (default: True)

**Returns:** HATS autocorrelation value (float, typically close to zero)

**Requirements:**
- No ORCA calculation required (uses RDKit only)
- Works with any molecule

**Limitations:**
- Returns 0.0 if molecule has fewer than ``lag + 1`` atoms
- Only considers atoms capable of hydrogen bonding (N, O, F)
- Returns 0.0 if no H-bond capable atoms are found

**Implementation Notes:**
- Focuses on hydrogen-bonding atoms (N, O, F)
- Useful for modeling hydrogen-bonding interactions

Summary of Requirements
-----------------------

Method Type Requirements
~~~~~~~~~~~~~~~~~~~~~~~~

- **SP (Single Point):** Works for most energy and electronic descriptors
- **Opt (Optimization):** Required for geometric descriptors (volume, SASA, PSA, bond lengths)
- **Freq (Frequency):** Required for thermodynamic descriptors (entropy, thermal corrections)

Functional Requirements
~~~~~~~~~~~~~~~~~~~~~~~

- **DFT methods:** All descriptors work
- **Semi-empirical methods:** All descriptors work, but energy values are less accurate
- **Basis set:** Only relevant for DFT methods, ignored for semi-empirical

Common Limitations
~~~~~~~~~~~~~~~~~~

1. **Geometric descriptors** require optimized geometries
2. **Thermodynamic descriptors** require frequency calculations
3. **Some descriptors** may return 0.0 or default values if data is unavailable
4. **Semi-empirical methods** provide less accurate energy values than DFT
5. **Atomic charges** are basis-set dependent (Mulliken charges)

Approximations and Assumptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Koopmans' theorem** is used for chemical potential and electronegativity (approximate for correlated systems)
2. **Frontier electron density** uses atomic charges as approximation
3. **SASA calculation** uses fixed atomic radii (may not be accurate for all systems)
4. **XY shadow** uses simple bounding box (does not account for atomic radii)
5. **MERIC** uses approximate charge-based calculation
6. **Molecular volume** fallback uses molecular weight estimation
7. **Thermodynamic corrections** may not be available for SP calculations

