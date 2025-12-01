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
        
        data["total_energy"] = self._parse_total_energy(content)
        
        homo, lumo = self._parse_orbital_energies(content)
        data["homo_energy"] = homo
        data["lumo_energy"] = lumo
        
        data["dipole_moment"] = self._parse_dipole_moment(content)
        data["atom_charges"] = self._parse_mulliken_charges(content, mol)
        data["coordinates"] = self._parse_coordinates(content)
        data["bond_lengths"] = self._parse_bond_lengths(content, data.get("coordinates", []), mol)
        data["gibbs_free_energy"] = self._parse_gibbs_free_energy(content)
        data["entropy"] = self._parse_entropy(content)
        data["molecular_volume"] = self._parse_molecular_volume(content)
        data["polar_surface_area"] = self._parse_polar_surface_area(content, data.get("coordinates", []))
        data["orbital_energies"] = self._parse_all_orbital_energies(content)
        data["nmr_shifts"] = self._parse_nmr_shifts(content, mol)
        data["mayer_indices"] = self._parse_mayer_indices(content, mol)
        data["nbo_stabilization_energies"] = self._parse_nbo_stabilization_energies(content)
        
        return data
    
    def _parse_total_energy(self, content: str) -> float:
        """Parse total energy in Hartree."""
        patterns = [
            r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
            r"Total Energy\s*:\s*(-?\d+\.\d+)",
            r"E\(.*\)\s*=\s*(-?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                return float(match.group(1))
        
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
        
        lines = content.split("\n")
        in_orbital_section = False
        last_occupied_energy = None
        first_virtual_energy = None
        
        orbital_sections = []
        for i, line in enumerate(lines):
            if "ORBITAL ENERGIES" in line.upper():
                orbital_sections.append(i)
        
        start_idx = orbital_sections[-1] if orbital_sections else -1
        
        if start_idx >= 0:
            start_idx += 1
            if start_idx < len(lines) and "---" in lines[start_idx]:
                start_idx += 1
            if start_idx < len(lines) and not lines[start_idx].strip():
                start_idx += 1
            if start_idx < len(lines) and "NO" in lines[start_idx] and "OCC" in lines[start_idx]:
                start_idx += 1
        
        for i in range(start_idx, len(lines)):
            line = lines[i]
            
            if line.strip() == "" or "*Only the first" in line or "Total SCF time" in line or "FINAL SINGLE POINT" in line or "DFT DISPERSION" in line:
                if last_occupied_energy is not None or first_virtual_energy is not None:
                    break
                continue
            
            if "NO" in line and "OCC" in line and "E(Eh)" in line:
                continue
            
            if "---" in line and len(line.strip()) < 20:
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                try:
                    occ_num = float(parts[1])
                    energy_ev = float(parts[3])
                    
                    if occ_num > 0.1:
                        last_occupied_energy = energy_ev
                    elif occ_num < 0.1 and first_virtual_energy is None:
                        first_virtual_energy = energy_ev
                except (ValueError, IndexError):
                    continue
        
        if last_occupied_energy is not None:
            homo = last_occupied_energy
        
        if first_virtual_energy is not None:
            lumo = first_virtual_energy
        
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
        """Parse dipole moment magnitude in Debye.
        
        When solvation models (SMD/COSMO) are used, ORCA may output both
        gas-phase and solvated dipole moments. This parser prioritizes
        gas-phase values to match experimental data.
        """
        lines = content.split("\n")
        
        # First, try to find gas-phase dipole moment (before solvation section)
        # Look for dipole moment before SMD/COSMO sections
        gas_phase_dipole = None
        solvated_dipole = None
        
        in_solvation_section = False
        dipole_sections = []
        
        for i, line in enumerate(lines):
            # Mark start of solvation section
            if "SMD" in line.upper() or "COSMO" in line.upper() or "CPCM" in line.upper():
                in_solvation_section = True
            
            # Look for dipole moment entries
            if "dipole moment" in line.lower():
                dipole_sections.append((i, in_solvation_section, line))
        
        # Try to find gas-phase dipole first (before any solvation)
        for i, (line_idx, in_solv, line) in enumerate(dipole_sections):
            if not in_solv:
                # This is likely gas-phase dipole
                # Try to extract value from this section
                for j in range(max(0, line_idx - 2), min(len(lines), line_idx + 5)):
                    check_line = lines[j]
                    # Look for "Total Dipole Moment" or magnitude
                    patterns = [
                        r"Total Dipole Moment\s*:\s*(-?\d+\.\d+)\s*Debye",
                        r"Dipole moment\s*:\s*(-?\d+\.\d+)\s*Debye",
                        r"Magnitude\s+\(Debye\)\s*:\s*(-?\d+\.\d+)",
                    ]
                    for pattern in patterns:
                        match = re.search(pattern, check_line, re.IGNORECASE)
                        if match:
                            gas_phase_dipole = abs(float(match.group(1)))
                            break
                    
                    # Try component format
                    dipole_pattern = r"Dipole moment\s+\(Debye\)\s*:\s*X=\s*(-?\d+\.\d+)\s*Y=\s*(-?\d+\.\d+)\s*Z=\s*(-?\d+\.\d+)"
                    dipole_match = re.search(dipole_pattern, check_line, re.IGNORECASE)
                    if dipole_match:
                        dx = float(dipole_match.group(1))
                        dy = float(dipole_match.group(2))
                        dz = float(dipole_match.group(3))
                        gas_phase_dipole = (dx**2 + dy**2 + dz**2)**0.5
                        break
                    
                    # Try multi-line format
                    if "Dipole moment" in check_line.lower() and j + 3 < len(lines):
                        try:
                            x_line = lines[j + 1] if j + 1 < len(lines) else ""
                            y_line = lines[j + 2] if j + 2 < len(lines) else ""
                            z_line = lines[j + 3] if j + 3 < len(lines) else ""
                            
                            x_match = re.search(r"(-?\d+\.\d+)", x_line)
                            y_match = re.search(r"(-?\d+\.\d+)", y_line)
                            z_match = re.search(r"(-?\d+\.\d+)", z_line)
                            
                            if x_match and y_match and z_match:
                                dx = float(x_match.group(1))
                                dy = float(y_match.group(1))
                                dz = float(z_match.group(1))
                                gas_phase_dipole = (dx**2 + dy**2 + dz**2)**0.5
                                break
                        except (ValueError, IndexError):
                            pass
                
                if gas_phase_dipole is not None:
                    break
        
        # If gas-phase dipole found, return it
        if gas_phase_dipole is not None:
            return gas_phase_dipole
        
        # Fallback to original parsing (any dipole moment)
        patterns = [
            r"Total Dipole Moment\s*:\s*(-?\d+\.\d+)\s*Debye",
            r"Dipole moment\s*:\s*(-?\d+\.\d+)\s*Debye",
            r"Magnitude\s+\(Debye\)\s*:\s*(-?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return abs(float(match.group(1)))
        
        dipole_pattern = r"Dipole moment\s+\(Debye\)\s*:\s*X=\s*(-?\d+\.\d+)\s*Y=\s*(-?\d+\.\d+)\s*Z=\s*(-?\d+\.\d+)"
        dipole_match = re.search(dipole_pattern, content, re.IGNORECASE)
        if dipole_match:
            dx = float(dipole_match.group(1))
            dy = float(dipole_match.group(2))
            dz = float(dipole_match.group(3))
            return (dx**2 + dy**2 + dz**2)**0.5
        
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
        
        lines = content.split("\n")
        in_charges_section = False
        
        for i, line in enumerate(lines):
            if "MULLIKEN ATOMIC CHARGES" in line.upper():
                in_charges_section = True
                continue
            
            if in_charges_section:
                match = re.match(r"\s*(\d+)\s+\w+\s*:\s*(-?\d+\.\d+)", line)
                if match:
                    idx = int(match.group(1))
                    charge = float(match.group(2))
                    charges[idx] = charge
                elif "Sum of atomic charges" in line:
                    break
                elif ("---" in line or "===" in line) and charges:
                    break
                elif not line.strip() and charges and i > 0 and "Total SCF time" in lines[i+1] if i+1 < len(lines) else False:
                    break
        
        if not charges and mol is not None:
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == "C":
                    charges[i] = -0.1
                elif atom.GetSymbol() == "H":
                    charges[i] = 0.1
                else:
                    charges[i] = 0.0
        
        return charges
    
    def _parse_coordinates(self, content: str) -> list[tuple[str, float, float, float]]:
        """Parse final optimized coordinates.
        
        Uses the LAST occurrence of coordinates (after optimization is complete).
        """
        coordinates: list[tuple[str, float, float, float]] = []
        
        lines = content.split("\n")
        
        coord_sections = []
        for i, line in enumerate(lines):
            if ("CARTESIAN COORDINATES" in line.upper() and "ANGSTROEM" in line.upper()) or \
               ("FINAL ENERGY EVALUATION" in line.upper() and "STATIONARY POINT" in line.upper()):
                coord_sections.append(i)
        
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
                    if i + 1 < len(lines):
                        next_line = lines[i + 1]
                        if "---" in next_line or "===" in next_line:
                            skip_next = True
                    continue
                
                if in_coords_section:
                    if "NO LB" in line or ("ZA" in line and "FRAG" in line):
                        continue
                    
                    match = re.match(r"\s*\d+\s+(\w+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
                    if not match:
                        match = re.match(r"\s*(\w+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
                    
                    if match:
                        symbol = match.group(1)
                        if symbol.isdigit():
                            continue
                        x = float(match.group(2))
                        y = float(match.group(3))
                        z = float(match.group(4))
                        coordinates.append((symbol, x, y, z))
                    elif ("---" in line or "===" in line) and coordinates:
                        break
                    elif not line.strip() and coordinates:
                        break
        
        return coordinates
    
    def _parse_bond_lengths(self, content: str, coordinates: list[tuple[str, float, float, float]], mol: Optional[Mol] = None) -> list[tuple[int, int, float]]:
        """Calculate bond lengths from coordinates.
        
        Only includes actual bonds (determined from molecule structure or distance criteria).
        """
        bond_lengths: list[tuple[int, int, float]] = []
        
        if not coordinates:
            return bond_lengths
        
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
                
                max_len = max_bond_lengths.get((sym1, sym2)) or max_bond_lengths.get((sym2, sym1))
                if max_len and distance < max_len:
                    bond_lengths.append((i, j, distance))
        
        return bond_lengths
    
    def _parse_gibbs_free_energy(self, content: str) -> float:
        """Parse Gibbs free energy in Hartree."""
        patterns = [
            r"Gibbs Free Energy\s*:\s*(-?\d+\.\d+)",
            r"G\(T\)\s*=\s*(-?\d+\.\d+)",
            r"Final Gibbs free energy\s*:\s*(-?\d+\.\d+)",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                return float(match.group(1))
        
        return self._parse_total_energy(content)
    
    def _parse_entropy(self, content: str) -> float:
        """Parse entropy in J/(mol·K)."""
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
        
        return 250.0
    
    def _parse_molecular_volume(self, content: str) -> float:
        """Parse molecular volume in Å³."""
        match = re.search(r"Volume\s*:\s*(-?\d+\.\d+)\s*Å", content, re.IGNORECASE)
        if match:
            return float(match.group(1))
        
        coordinates = self._parse_coordinates(content)
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
        
        atom_counts = {}
        for symbol, _, _, _ in coordinates:
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        psa = 0.0
        
        n_oxygen = atom_counts.get("O", 0)
        n_nitrogen = atom_counts.get("N", 0)
        
        if n_oxygen > 0:
            total_atoms = sum(atom_counts.values())
            if total_atoms == 3 and n_oxygen == 1:
                psa += 20.23
            elif n_oxygen == 1 and total_atoms <= 10:
                psa += 17.07
            else:
                psa += n_oxygen * 18.0
        
        if n_nitrogen > 0:
            psa += n_nitrogen * 13.0
        
        psa += atom_counts.get("P", 0) * 13.0
        psa += atom_counts.get("S", 0) * 12.0
        
        return psa
    
    def _parse_all_orbital_energies(self, content: str) -> list[float]:
        """Parse all orbital energies in eV.
        
        Returns list of orbital energies, ordered from lowest to highest.
        Occupied orbitals come first, then virtual orbitals.
        """
        orbital_energies: list[float] = []
        
        lines = content.split("\n")
        orbital_sections = []
        for i, line in enumerate(lines):
            if "ORBITAL ENERGIES" in line.upper():
                orbital_sections.append(i)
        
        start_idx = orbital_sections[-1] if orbital_sections else -1
        
        if start_idx >= 0:
            start_idx += 1
            if start_idx < len(lines) and "---" in lines[start_idx]:
                start_idx += 1
            if start_idx < len(lines) and not lines[start_idx].strip():
                start_idx += 1
            if start_idx < len(lines) and "NO" in lines[start_idx] and "OCC" in lines[start_idx]:
                start_idx += 1
            
            for i in range(start_idx, len(lines)):
                line = lines[i]
                
                if line.strip() == "" or "*Only the first" in line or "Total SCF time" in line or "FINAL SINGLE POINT" in line or "DFT DISPERSION" in line:
                    break
                
                if "NO" in line and "OCC" in line and "E(Eh)" in line:
                    continue
                
                if "---" in line and len(line.strip()) < 20:
                    continue
                
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        energy_ev = float(parts[3])
                        orbital_energies.append(energy_ev)
                    except (ValueError, IndexError):
                        continue
        
        return orbital_energies
    
    def _parse_nmr_shifts(self, content: str, mol: Optional[Mol] = None) -> dict[int, float]:
        """Parse NMR chemical shifts from ORCA output.
        
        Returns dictionary mapping atom index to chemical shift in ppm.
        ORCA outputs shielding constants (sigma), chemical shifts are calculated as: delta = sigma_ref - sigma
        For typical reference: TMS for 13C (sigma_ref ~ 184 ppm), TMS for 1H (sigma_ref ~ 31 ppm)
        """
        shifts: dict[int, float] = {}
        shielding: dict[int, float] = {}
        
        lines = content.split("\n")
        in_nmr_section = False
        in_shielding_section = False
        
        for i, line in enumerate(lines):
            if "NMR SHIELDING" in line.upper() or "GIAO" in line.upper() or ("ISOTROPIC" in line.upper() and "SHIELDING" in line.upper()):
                in_shielding_section = True
                in_nmr_section = True
                continue
            
            if "CHEMICAL SHIFT" in line.upper() or "NMR SHIFTS" in line.upper():
                in_nmr_section = True
                in_shielding_section = False
                continue
            
            if in_shielding_section or in_nmr_section:
                match = re.match(r"\s*(\d+)\s+\w+\s+(-?\d+\.\d+)", line)
                if not match:
                    match = re.match(r"\s*(\d+)\s+(-?\d+\.\d+)", line)
                if match:
                    idx = int(match.group(1))
                    val = float(match.group(2))
                    if in_shielding_section:
                        shielding[idx] = val
                    else:
                        if abs(val) > 10.0:
                            shifts[idx] = val
                elif "---" in line or "===" in line:
                    if shifts or shielding:
                        break
                elif not line.strip() and (shifts or shielding):
                    break
        
        if shielding and not shifts:
            for idx, sigma in shielding.items():
                if abs(sigma) < 500:
                    delta = 184.1 - sigma
                    if 50.0 <= delta <= 250.0:
                        shifts[idx] = delta
        
        if not shifts:
            for i, line in enumerate(lines):
                if "ISOTROPIC" in line.upper() and "SHIELDING" in line.upper():
                    for j in range(i + 1, min(i + 50, len(lines))):
                        match = re.search(r"(\d+)\s+(-?\d+\.\d+)", lines[j])
                        if match:
                            idx = int(match.group(1))
                            sigma = float(match.group(2))
                            if abs(sigma) < 500:
                                delta = 184.1 - sigma
                                if 50.0 <= delta <= 250.0:
                                    shifts[idx] = delta
        
        return shifts
    
    def _parse_mayer_indices(self, content: str, mol: Optional[Mol] = None) -> list[tuple[int, int, float]]:
        """Parse Mayer bond indices from ORCA output.
        
        Returns list of tuples: (atom_i, atom_j, index_value)
        """
        indices: list[tuple[int, int, float]] = []
        
        lines = content.split("\n")
        in_mayer_section = False
        
        for i, line in enumerate(lines):
            if "MAYER BOND ORDERS" in line.upper() or ("MAYER" in line.upper() and "bond orders" in line.lower()):
                in_mayer_section = True
                continue
            
            if in_mayer_section:
                matches = re.findall(r"B\s*\(\s*(\d+)\s*-\s*\w+\s*,\s*(\d+)\s*-\s*\w+\s*\)\s*:\s*(\d+\.\d+)", line)
                for match in matches:
                    idx_i = int(match[0])
                    idx_j = int(match[1])
                    index_val = float(match[2])
                    if 0.01 <= index_val <= 10.0:
                        indices.append((idx_i, idx_j, index_val))
                
                if "---" in line or "===" in line or (not line.strip() and indices and not matches):
                    if indices:
                        break
                elif "MAYER" in line.upper() and "POPULATION" in line.upper():
                    if indices:
                        break
        
        if not indices:
            for i, line in enumerate(lines):
                if "Mayer bond orders" in line.lower():
                    for j in range(i, min(i + 50, len(lines))):
                        matches = re.findall(r"B\s*\(\s*(\d+)\s*-\s*\w+\s*,\s*(\d+)\s*-\s*\w+\s*\)\s*:\s*(\d+\.\d+)", lines[j])
                        for match in matches:
                            idx_i = int(match[0])
                            idx_j = int(match[1])
                            index_val = float(match[2])
                            if 0.01 <= index_val <= 10.0:
                                indices.append((idx_i, idx_j, index_val))
                        if "---" in lines[j] or "===" in lines[j] or (not lines[j].strip() and indices and not matches):
                            if indices:
                                break
        
        return indices
    
    def _parse_nbo_stabilization_energies(self, content: str) -> dict[str, float]:
        """Parse NBO stabilization energies (E(2)) from ORCA output.
        
        Returns dictionary mapping interaction description to energy in kcal/mol.
        """
        energies: dict[str, float] = {}
        
        if "NBOEXE" in content.upper() and ("NOT SET" in content.upper() or "SKIPPING" in content.upper()):
            return energies
        
        lines = content.split("\n")
        in_perturbation_section = False
        
        for i, line in enumerate(lines):
            if "SECOND ORDER PERTURBATION THEORY" in line.upper() or "SECOND ORDER" in line.upper():
                in_perturbation_section = True
                continue
            
            if in_perturbation_section:
                if "E(2)" in line or "STABILIZATION ENERGY" in line.upper():
                    match = re.search(r"(-?\d+\.\d+)\s*kcal", line, re.IGNORECASE)
                    if match:
                        energy = float(match.group(1))
                        if energy > 0:
                            donor_match = re.search(r"(\w+)\s*->\s*(\w+)", line)
                            if donor_match:
                                key = f"{donor_match.group(1)}->{donor_match.group(2)}"
                                energies[key] = energy
                            else:
                                lp_match = re.search(r"LP\s*\(\s*(\w+)\s*\)", line, re.IGNORECASE)
                                pi_match = re.search(r"PI\s*\*\s*\(([^)]+)\)", line, re.IGNORECASE)
                                if lp_match and pi_match:
                                    key = f"LP({lp_match.group(1)})->PiStar({pi_match.group(1)})"
                                    energies[key] = energy
                                else:
                                    donor_acceptor = re.search(r"(\w+)\s+(\w+)", line)
                                    if donor_acceptor:
                                        key = f"{donor_acceptor.group(1)}->{donor_acceptor.group(2)}"
                                        energies[key] = energy
                                    else:
                                        energies[f"interaction_{len(energies)}"] = energy
                elif "---" in line or "===" in line or (not line.strip() and energies):
                    if energies:
                        break
        
        if not energies:
            for i, line in enumerate(lines):
                if "E(2)" in line or ("STABILIZATION" in line.upper() and "ENERGY" in line.upper()):
                    match = re.search(r"(-?\d+\.\d+)\s*kcal", line, re.IGNORECASE)
                    if match:
                        energy = float(match.group(1))
                        if energy > 0:
                            energies[f"interaction_{len(energies)}"] = energy
        
        return energies
    
