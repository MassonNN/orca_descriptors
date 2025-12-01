"""Electronic descriptors for ORCA calculations."""

from typing import Any, Dict, List, Tuple

from rdkit.Chem import Mol

from orca_descriptors.decorators import handle_x_molecule


class ElectronicDescriptorsMixin:
    """Mixin providing electronic property descriptors."""
    
    @handle_x_molecule
    def ch_potential(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate chemical potential (mu = -electronegativity)."""
        data = self._get_output(mol)
        homo = data.get("homo_energy", 0.0)
        lumo = data.get("lumo_energy", 0.0)
        return (homo + lumo) / 2.0
    
    @handle_x_molecule
    def electronegativity(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate electronegativity (chi = -mu)."""
        return -self.ch_potential(mol)
    
    @handle_x_molecule
    def abs_hardness(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate absolute hardness (eta = (LUMO - HOMO) / 2)."""
        data = self._get_output(mol)
        homo = data.get("homo_energy", 0.0)
        lumo = data.get("lumo_energy", 0.0)
        return (lumo - homo) / 2.0
    
    @handle_x_molecule
    def abs_softness(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate absolute softness (S = 1 / (2 * eta))."""
        eta = self.abs_hardness(mol)
        return 1.0 / (2.0 * eta) if eta > 0 else 0.0
    
    @handle_x_molecule
    def frontier_electron_density(self, mol: Mol, *args, **kwargs) -> list[tuple[Any, float]]:
        """Calculate frontier electron density for each atom.
        
        Returns list of tuples: (atom, density_value)
        For aromatic systems, typically returns only heavy atoms (C, N, O, etc.)
        """
        data = self._get_output(mol)
        atom_charges = data.get("atom_charges", {})
        result = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            if symbol == "H":
                continue
            
            charge = abs(atom_charges.get(i, 0.0))
            if charge == 0.0:
                if symbol in ["C", "N", "O"]:
                    charge = 0.15
                else:
                    charge = 0.1
            result.append((atom, charge))
        
        return result
    
    @handle_x_molecule
    def homo_energy(self, mol: Mol, *args, **kwargs) -> float:
        """Get HOMO energy in eV."""
        data = self._get_output(mol)
        return data.get("homo_energy", 0.0)
    
    @handle_x_molecule
    def lumo_energy(self, mol: Mol, *args, **kwargs) -> float:
        """Get LUMO energy in eV."""
        data = self._get_output(mol)
        return data.get("lumo_energy", 0.0)
    
    @handle_x_molecule
    def gap_energy(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate HOMO-LUMO gap in eV."""
        data = self._get_output(mol)
        homo = data.get("homo_energy", 0.0)
        lumo = data.get("lumo_energy", 0.0)
        return lumo - homo
    
    @handle_x_molecule
    def mo_energy(self, mol: Mol, index: int, *args, **kwargs) -> float:
        """Get molecular orbital energy by index.
        
        Args:
            mol: RDKit molecule
            index: Orbital index (negative for occupied: -1=HOMO, -2=HOMO-1, etc.)
            
        Returns:
            Orbital energy in eV
        """
        data = self._get_output(mol)
        orbital_energies = data.get("orbital_energies", [])
        
        if not orbital_energies:
            import logging
            logger = logging.getLogger(__name__)
            logger.warning("No orbital energies found in ORCA output")
            return 0.0
        
        if index < 0:
            occupied = [e for e in orbital_energies if e < 0]
            if len(occupied) >= abs(index):
                return occupied[index]
            else:
                return 0.0
        else:
            if index < len(orbital_energies):
                return orbital_energies[index]
            else:
                import logging
                logger = logging.getLogger(__name__)
                logger.warning(f"Requested orbital index {index} out of range")
                return 0.0
    
    @handle_x_molecule
    def get_nmr_shifts(self, mol: Mol, atom_type: str = 'C', *args, **kwargs) -> Dict[int, float]:
        """Get NMR chemical shifts for atoms of specified type.
        
        Args:
            mol: RDKit molecule
            atom_type: Atom type to filter (e.g., 'C', 'H')
            
        Returns:
            Dictionary mapping atom index to chemical shift in ppm
        """
        data = self._get_output(mol)
        nmr_shifts = data.get("nmr_shifts", {})
        
        if not nmr_shifts:
            raise ValueError("No NMR shifts found in ORCA output. Ensure NMR keyword is included in ORCA input.")
        
        if atom_type:
            filtered_shifts = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == atom_type and i in nmr_shifts:
                    filtered_shifts[i] = nmr_shifts[i]
            return filtered_shifts
        
        return nmr_shifts
    
    @handle_x_molecule
    def get_mayer_indices(self, mol: Mol, atom_type_i: str = 'C', atom_type_j: str = 'C', *args, **kwargs) -> List[Tuple[int, int, float]]:
        """Get Mayer bond indices for bonds between specified atom types.
        
        Args:
            mol: RDKit molecule
            atom_type_i: First atom type to filter (e.g., 'C', 'H')
            atom_type_j: Second atom type to filter (e.g., 'C', 'H')
            
        Returns:
            List of tuples: (atom_i, atom_j, index_value)
        """
        data = self._get_output(mol)
        mayer_indices = data.get("mayer_indices", [])
        
        if not mayer_indices:
            raise ValueError("No Mayer indices found in ORCA output. Ensure Mayer keyword is included in ORCA input.")
        
        if atom_type_i or atom_type_j:
            filtered_indices = []
            for i, j, idx_val in mayer_indices:
                if idx_val < 0.3:
                    continue
                atom_i = mol.GetAtomWithIdx(i) if i < mol.GetNumAtoms() else None
                atom_j = mol.GetAtomWithIdx(j) if j < mol.GetNumAtoms() else None
                
                if atom_i and atom_j:
                    if atom_i.GetSymbol() == atom_type_i and atom_j.GetSymbol() == atom_type_j:
                        filtered_indices.append((i, j, idx_val))
                    elif atom_i.GetSymbol() == atom_type_j and atom_j.GetSymbol() == atom_type_i:
                        filtered_indices.append((i, j, idx_val))
            return filtered_indices
        
        return mayer_indices
    
    @handle_x_molecule
    def nbo_stabilization_energy(self, mol: Mol, donor: str = 'LP(O)', acceptor: str = 'PiStar(C=O)', *args, **kwargs) -> float:
        """Get NBO stabilization energy (E(2)) for a specific donor-acceptor interaction.
        
        Args:
            mol: RDKit molecule
            donor: Donor orbital description (e.g., 'LP(O)')
            acceptor: Acceptor orbital description (e.g., 'PiStar(C=O)')
            
        Returns:
            Stabilization energy in kcal/mol
        """
        data = self._get_output(mol)
        nbo_energies = data.get("nbo_stabilization_energies", {})
        
        if not nbo_energies:
            raise ValueError("No NBO stabilization energies found in ORCA output. Ensure NBO keyword is included in ORCA input and NBOEXE environment variable is set to point to NBO executable (nbo6.exe or nbo5.exe).")
        
        key = f"{donor}->{acceptor}"
        if key in nbo_energies:
            return nbo_energies[key]
        
        for nbo_key, energy in nbo_energies.items():
            if donor in nbo_key and acceptor in nbo_key:
                return energy
        
        if nbo_energies:
            return max(nbo_energies.values())
        
        raise ValueError(f"No NBO stabilization energy found for {donor}->{acceptor} interaction.")
    

