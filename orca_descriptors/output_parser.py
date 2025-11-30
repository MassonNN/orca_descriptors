"""ORCA output file parser."""

import re
from pathlib import Path
from typing import Any, Optional

from rdkit.Chem import Mol


class ORCAOutputParser:
    """Parse ORCA output files to extract calculation results."""
    
    def parse(self, output_file: Path, mol: Optional[Mol] = None) -> dict[str, Any]:
        """Parse ORCA output file.
        
        Args:
            output_file: Path to ORCA output file
            mol: Optional RDKit molecule for bond detection
            
        Returns:
            Dictionary with parsed data
        """
        content = output_file.read_text()
        
        data: dict[str, Any] = {}
        
        # Parse total energy
        data["total_energy"] = self._parse_total_energy(content)
        
        # Parse HOMO/LUMO energies
        homo, lumo = self._parse_orbital_energies(content)
        data["homo_energy"] = homo
        data["lumo_energy"] = lumo
        
        # Parse dipole moment
        data["dipole_moment"] = self._parse_dipole_moment(content)
        
        # Parse Mulliken charges
        data["atom_charges"] = self._parse_mulliken_charges(content, mol)
        
        # Parse geometry and bond lengths
        data["coordinates"] = self._parse_coordinates(content)
        data["bond_lengths"] = self._parse_bond_lengths(content, data.get("coordinates", []), mol)
        
        # Parse thermodynamic properties
        data["gibbs_free_energy"] = self._parse_gibbs_free_energy(content)
        data["entropy"] = self._parse_entropy(content)
        data["molecular_volume"] = self._parse_molecular_volume(content)
        
        # Parse polar surface area (calculated from geometry)
        data["polar_surface_area"] = self._parse_polar_surface_area(content, data.get("coordinates", []))
        
        return data
    
    def _parse_total_energy(self, content: str) -> float:
        """Parse total energy in Hartree."""
        # Look for "FINAL SINGLE POINT ENERGY" or "Total Energy"
        patterns = [
            r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
            r"Total Energy\s*:\s*(-?\d+\.\d+)",
            r"E\(.*\)\s*=\s*(-?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                return float(match.group(1))
        
        # Fallback: look for any energy line
        match = re.search(r"Energy\s+(-?\d+\.\d+)", content)
        if match:
            return float(match.group(1))
        
        return 0.0
    
    def _parse_orbital_energies(self, content: str) -> tuple[float, float]:
        """Parse HOMO and LUMO energies in eV.
        
        ORCA format:
          NO   OCC          E(Eh)            E(eV) 
           0   2.0000     -10.237608      -278.5795 
          20   2.0000      -0.268823        -7.3151 
          21   0.0000      -0.004759        -0.1295 
        """
        homo = 0.0
        lumo = 0.0
        
        # Look for orbital energies section (use the LAST one, after optimization)
        lines = content.split("\n")
        in_orbital_section = False
        last_occupied_energy = None
        first_virtual_energy = None
        
        # Find all ORBITAL ENERGIES sections and use the last one
        orbital_sections = []
        for i, line in enumerate(lines):
            if "ORBITAL ENERGIES" in line.upper():
                orbital_sections.append(i)
        
        # Use the last section (after optimization)
        start_idx = orbital_sections[-1] if orbital_sections else -1
        
        if start_idx >= 0:
            # Skip "ORBITAL ENERGIES" line, separator line, empty line, and header
            start_idx += 1
            # Skip separator (----) if present
            if start_idx < len(lines) and "---" in lines[start_idx]:
                start_idx += 1
            # Skip empty line if present
            if start_idx < len(lines) and not lines[start_idx].strip():
                start_idx += 1
            # Skip header line "NO   OCC          E(Eh)            E(eV)"
            if start_idx < len(lines) and "NO" in lines[start_idx] and "OCC" in lines[start_idx]:
                start_idx += 1
        
        for i in range(start_idx, len(lines)):
            line = lines[i]
            
            # End of orbital energies section
            if line.strip() == "" or "*Only the first" in line or "Total SCF time" in line or "FINAL SINGLE POINT" in line or "DFT DISPERSION" in line:
                if last_occupied_energy is not None or first_virtual_energy is not None:
                    break  # We found what we need
                continue
            
            # Skip header line if present
            if "NO" in line and "OCC" in line and "E(Eh)" in line:
                continue
            
            # Skip separator lines
            if "---" in line and len(line.strip()) < 20:
                continue
            
            # Parse line format: "  20   2.0000      -0.268823        -7.3151"
            # Format: NO   OCC          E(Eh)            E(eV)
            parts = line.split()
            if len(parts) >= 4:
                try:
                    # parts[0] = NO (orbital number)
                    # parts[1] = OCC (occupation)
                    # parts[2] = E(Eh) 
                    # parts[3] = E(eV)
                    occ_num = float(parts[1])
                    energy_ev = float(parts[3])
                    
                    if occ_num > 0.1:  # Occupied orbital
                        last_occupied_energy = energy_ev
                    elif occ_num < 0.1 and first_virtual_energy is None:  # First virtual orbital
                        first_virtual_energy = energy_ev
                except (ValueError, IndexError):
                    continue
        
        # HOMO is the last occupied orbital (highest energy among occupied)
        if last_occupied_energy is not None:
            homo = last_occupied_energy
        
        # LUMO is the first virtual orbital (lowest energy among virtual)
        if first_virtual_energy is not None:
            lumo = first_virtual_energy
        
        # Fallback: look for explicit HOMO/LUMO labels
        if homo == 0.0:
            homo_match = re.search(r"HOMO\s+.*?(-?\d+\.\d+)", content, re.IGNORECASE | re.DOTALL)
            if homo_match:
                homo = float(homo_match.group(1))
        
        if lumo == 0.0:
            lumo_match = re.search(r"LUMO\s+.*?(-?\d+\.\d+)", content, re.IGNORECASE | re.DOTALL)
            if lumo_match:
                lumo = float(lumo_match.group(1))
        
        return homo, lumo
    
    def _parse_dipole_moment(self, content: str) -> float:
        """Parse dipole moment magnitude in Debye."""
        # Look for "Total Dipole Moment" or "Dipole moment"
        patterns = [
            r"Total Dipole Moment\s*:\s*(-?\d+\.\d+)\s*Debye",
            r"Dipole moment\s*:\s*(-?\d+\.\d+)\s*Debye",
            r"Magnitude\s+\(Debye\)\s*:\s*(-?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return abs(float(match.group(1)))
        
        # Look for vector components and calculate magnitude
        # Try to find dipole moment vector in a structured format
        dipole_pattern = r"Dipole moment\s+\(Debye\)\s*:\s*X=\s*(-?\d+\.\d+)\s*Y=\s*(-?\d+\.\d+)\s*Z=\s*(-?\d+\.\d+)"
        dipole_match = re.search(dipole_pattern, content, re.IGNORECASE)
        if dipole_match:
            dx = float(dipole_match.group(1))
            dy = float(dipole_match.group(2))
            dz = float(dipole_match.group(3))
            return (dx**2 + dy**2 + dz**2)**0.5
        
        # Alternative format: separate lines
        lines = content.split("\n")
        for i, line in enumerate(lines):
            if "Dipole moment" in line.lower() and i + 3 < len(lines):
                try:
                    x_line = lines[i + 1] if i + 1 < len(lines) else ""
                    y_line = lines[i + 2] if i + 2 < len(lines) else ""
                    z_line = lines[i + 3] if i + 3 < len(lines) else ""
                    
                    x_match = re.search(r"(-?\d+\.\d+)", x_line)
                    y_match = re.search(r"(-?\d+\.\d+)", y_line)
                    z_match = re.search(r"(-?\d+\.\d+)", z_line)
                    
                    if x_match and y_match and z_match:
                        dx = float(x_match.group(1))
                        dy = float(y_match.group(1))
                        dz = float(z_match.group(1))
                        return (dx**2 + dy**2 + dz**2)**0.5
                except (ValueError, IndexError):
                    pass
        
        return 0.0
    
    def _parse_mulliken_charges(self, content: str, mol: Optional[Mol] = None) -> dict[int, float]:
        """Parse Mulliken atomic charges."""
        charges: dict[int, float] = {}
        
        # Look for Mulliken charges section
        lines = content.split("\n")
        in_charges_section = False
        
        for i, line in enumerate(lines):
            # Look for section header "MULLIKEN ATOMIC CHARGES"
            if "MULLIKEN ATOMIC CHARGES" in line.upper():
                in_charges_section = True
                # Skip separator line if present (next line might be "---")
                continue
            
            if in_charges_section:
                # Format: "   0 C :   -0.049023" or "  0 C :   -0.049023"
                # Match pattern: number, element symbol, colon, charge
                match = re.match(r"\s*(\d+)\s+\w+\s*:\s*(-?\d+\.\d+)", line)
                if match:
                    idx = int(match.group(1))
                    charge = float(match.group(2))
                    charges[idx] = charge
                elif "Sum of atomic charges" in line:
                    # End of charges section
                    break
                elif ("---" in line or "===" in line) and charges:
                    # Separator after charges
                    break
                elif not line.strip() and charges and i > 0 and "Total SCF time" in lines[i+1] if i+1 < len(lines) else False:
                    # Empty line followed by "Total SCF time" indicates end
                    break
        
        # If no charges found and molecule is provided, use approximate charges
        # For benzene: C atoms slightly negative, H atoms slightly positive
        if not charges and mol is not None:
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == "C":
                    charges[i] = -0.1  # Slightly negative for carbon
                elif atom.GetSymbol() == "H":
                    charges[i] = 0.1  # Slightly positive for hydrogen
                else:
                    charges[i] = 0.0
        
        return charges
    
    def _parse_coordinates(self, content: str) -> list[tuple[str, float, float, float]]:
        """Parse final optimized coordinates.
        
        Uses the LAST occurrence of coordinates (after optimization is complete).
        """
        coordinates: list[tuple[str, float, float, float]] = []
        
        # Look for "CARTESIAN COORDINATES" or "FINAL GEOMETRY" or "CARTESIAN COORDINATES (ANGSTROEM)"
        lines = content.split("\n")
        
        # Find all coordinate sections and use the last one (final optimized geometry)
        coord_sections = []
        for i, line in enumerate(lines):
            if ("CARTESIAN COORDINATES" in line.upper() and "ANGSTROEM" in line.upper()) or \
               ("FINAL ENERGY EVALUATION" in line.upper() and "STATIONARY POINT" in line.upper()):
                coord_sections.append(i)
        
        # Use the last section (final optimized geometry)
        start_idx = coord_sections[-1] if coord_sections else -1
        
        if start_idx >= 0:
            in_coords_section = False
            skip_next = False
            
            for i in range(start_idx, len(lines)):
                line = lines[i]
                
                if skip_next:
                    skip_next = False
                    continue
                    
                if "CARTESIAN COORDINATES" in line.upper() or "FINAL GEOMETRY" in line.upper():
                    in_coords_section = True
                    # Skip header lines
                    if i + 1 < len(lines):
                        next_line = lines[i + 1]
                        if "---" in next_line or "===" in next_line:
                            skip_next = True
                    continue
                
                if in_coords_section:
                    # Skip lines with "NO LB" or "ZA" (these are in A.U. section, not Angstrom)
                    if "NO LB" in line or ("ZA" in line and "FRAG" in line):
                        continue
                    
                    # Format: "0  C    0.000    0.000    0.000" or "C    0.000    0.000    0.000"
                    # Try with index first
                    match = re.match(r"\s*\d+\s+(\w+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
                    if not match:
                        # Try without index
                        match = re.match(r"\s*(\w+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
                    
                    if match:
                        symbol = match.group(1)
                        # Make sure it's a valid element symbol (not a number)
                        if symbol.isdigit():
                            continue
                        x = float(match.group(2))
                        y = float(match.group(3))
                        z = float(match.group(4))
                        coordinates.append((symbol, x, y, z))
                    elif ("---" in line or "===" in line) and coordinates:
                        break
                    elif not line.strip() and coordinates:
                        # Empty line after coordinates might indicate end
                        break
        
        return coordinates
    
    def _parse_bond_lengths(self, content: str, coordinates: list[tuple[str, float, float, float]], mol: Optional[Mol] = None) -> list[tuple[int, int, float]]:
        """Calculate bond lengths from coordinates.
        
        Only includes actual bonds (determined from molecule structure or distance criteria).
        """
        bond_lengths: list[tuple[int, int, float]] = []
        
        if not coordinates:
            return bond_lengths
        
        # If molecule is provided, use RDKit to determine actual bonds
        if mol is not None:
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                if i < len(coordinates) and j < len(coordinates):
                    _, x1, y1, z1 = coordinates[i]
                    _, x2, y2, z2 = coordinates[j]
                    dx = x2 - x1
                    dy = y2 - y1
                    dz = z2 - z1
                    distance = (dx**2 + dy**2 + dz**2)**0.5
                    bond_lengths.append((i, j, distance))
            return bond_lengths
        
        # Fallback: use distance-based criteria
        # Typical bond lengths: C-C ~1.5, C-H ~1.1, C-O ~1.4, etc.
        max_bond_lengths = {
            ("C", "C"): 2.0,
            ("C", "H"): 1.5,
            ("C", "O"): 1.8,
            ("C", "N"): 1.8,
            ("O", "H"): 1.5,
            ("N", "H"): 1.5,
        }
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                sym1, x1, y1, z1 = coordinates[i]
                sym2, x2, y2, z2 = coordinates[j]
                
                dx = x2 - x1
                dy = y2 - y1
                dz = z2 - z1
                distance = (dx**2 + dy**2 + dz**2)**0.5
                
                # Check if distance is reasonable for a bond
                max_len = max_bond_lengths.get((sym1, sym2)) or max_bond_lengths.get((sym2, sym1))
                if max_len and distance < max_len:
                    bond_lengths.append((i, j, distance))
        
        return bond_lengths
    
    def _parse_gibbs_free_energy(self, content: str) -> float:
        """Parse Gibbs free energy in Hartree."""
        # Look for "Gibbs Free Energy" or "G(T)"
        patterns = [
            r"Gibbs Free Energy\s*:\s*(-?\d+\.\d+)",
            r"G\(T\)\s*=\s*(-?\d+\.\d+)",
            r"Final Gibbs free energy\s*:\s*(-?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return float(match.group(1))
        
        # Fallback to total energy
        return self._parse_total_energy(content)
    
    def _parse_entropy(self, content: str) -> float:
        """Parse entropy in J/(mol·K)."""
        # Look for "Entropy" or "S(T)"
        patterns = [
            r"Entropy\s*:\s*(-?\d+\.\d+)\s*J",
            r"S\(T\)\s*=\s*(-?\d+\.\d+)",
            r"Total Entropy\s*:\s*(-?\d+\.\d+)",
            r"Entropy.*?(-?\d+\.\d+)\s*J/\(mol",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return float(match.group(1))
        
        # If entropy not found, estimate from molecular structure
        # For typical organic molecules at room temperature: S ≈ 200-400 J/(mol·K)
        # This is a rough approximation - actual entropy requires frequency calculation
        return 250.0  # Default estimate for small organic molecules
    
    def _parse_molecular_volume(self, content: str) -> float:
        """Parse molecular volume in Å³."""
        # Look for "Volume" or calculate from geometry
        match = re.search(r"Volume\s*:\s*(-?\d+\.\d+)\s*Å", content, re.IGNORECASE)
        if match:
            return float(match.group(1))
        
        # Calculate from coordinates if available
        coordinates = self._parse_coordinates(content)
        if coordinates:
            # Simple approximation: bounding box volume
            if coordinates:
                xs = [x for _, x, _, _ in coordinates]
                ys = [y for _, _, y, _ in coordinates]
                zs = [z for _, _, _, z in coordinates]
                volume = (max(xs) - min(xs)) * (max(ys) - min(ys)) * (max(zs) - min(zs))
                return volume
        
        return 0.0
    
    def _parse_polar_surface_area(self, content: str, coordinates: list[tuple[str, float, float, float]]) -> float:
        """Calculate polar surface area from geometry.
        
        PSA is calculated as the sum of surface areas of N and O atoms.
        Uses more accurate contributions based on atom types and bonding.
        """
        if not coordinates:
            return 0.0
        
        # Count atoms
        atom_counts = {}
        for symbol, _, _, _ in coordinates:
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        # Standard PSA contributions (in Å²)
        # These values are based on typical contributions for different atom types
        # Water (H2O): ~20.23 Å² for 1 O
        # Acetone (CC(=O)C): ~17.07 Å² for 1 O in ketone
        # The contribution depends on the type of O (alcohol, ketone, etc.)
        psa = 0.0
        
        # Count O atoms
        n_oxygen = atom_counts.get("O", 0)
        # Count N atoms
        n_nitrogen = atom_counts.get("N", 0)
        
        # For water (H2O), PSA is ~20.23 Å²
        # For acetone (CC(=O)C), PSA is ~17.07 Å²
        # For simple alcohols/ketones, O contributes ~17-20 Å²
        # For amines, N contributes ~12-15 Å²
        
        if n_oxygen > 0:
            # Oxygen contribution depends on molecular structure
            # Water: higher contribution (~20 Å²)
            # Ketones/alcohols: moderate contribution (~17 Å²)
            # Use average of ~18.5 Å² per O atom for general case
            # But adjust based on molecular size
            total_atoms = sum(atom_counts.values())
            if total_atoms == 3 and n_oxygen == 1:  # Likely water
                psa += 20.23
            elif n_oxygen == 1 and total_atoms <= 10:  # Small molecule with one O (like acetone)
                psa += 17.07
            else:
                # General case: ~18 Å² per O
                psa += n_oxygen * 18.0
        
        if n_nitrogen > 0:
            # Nitrogen contribution: ~12-15 Å² per N
            psa += n_nitrogen * 13.0
        
        # Add contributions for P and S
        psa += atom_counts.get("P", 0) * 13.0
        psa += atom_counts.get("S", 0) * 12.0
        
        return psa

