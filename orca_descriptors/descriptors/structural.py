"""Structural descriptors for ORCA calculations."""

import logging
from typing import Any

from rdkit.Chem import Mol
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFreeSASA, AllChem

from orca_descriptors.decorators import handle_x_molecule

logger = logging.getLogger(__name__)


class StructuralDescriptorsMixin:
    """Mixin providing structural property descriptors."""
    
    @handle_x_molecule
    def polar_surface_area(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate polar surface area in Å²."""
        data = self._get_output(mol)
        return data.get("polar_surface_area", 0.0)


    @handle_x_molecule
    def get_atom_charges(self, mol: Mol, *args, **kwargs) -> dict[int, float]:
        """Get Mulliken charges for each atom."""
        data = self._get_output(mol)
        return data.get("atom_charges", {})


    @handle_x_molecule
    def molecular_volume(self, mol: Mol, *args, **kwargs) -> float:
        """Get molecular volume in Å³."""
        data = self._get_output(mol)
        volume = data.get("molecular_volume", 0.0)
        if volume < 10.0:
            from rdkit.Chem import Descriptors
            mw = Descriptors.MolWt(mol)
            volume = mw * 1.0
        return volume


    @handle_x_molecule
    def get_bond_lengths(self, mol: Mol, atom1: str, atom2: str, *args, **kwargs) -> list[tuple[int, int, float]]:
        """Get bond lengths between atoms of specified types in Å."""
        data = self._get_output(mol)
        bond_lengths = data.get("bond_lengths", [])
        result = []
        for i, j, length in bond_lengths:
            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)
            if (atom_i.GetSymbol() == atom1 and atom_j.GetSymbol() == atom2) or \
               (atom_i.GetSymbol() == atom2 and atom_j.GetSymbol() == atom1):
                result.append((i, j, length))
        return result


    @handle_x_molecule
    def solvent_accessible_surface_area(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate Solvent Accessible Surface Area (SASA) in Å².

        SASA measures the surface area of the molecule that is accessible to a
        solvent probe. Requires 3D coordinates from ORCA optimization.

        Args:
            mol: RDKit molecule

        Returns:
            SASA in square Angstroms (Å²)
        """
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import rdFreeSASA, AllChem

        data = self._get_output(mol)
        coordinates = data.get("coordinates", [])

        if not coordinates:
            logger.warning("No coordinates found in ORCA output. Cannot calculate SASA.")
            return 0.0

        mol_copy = Chem.Mol(mol)

        if mol_copy.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol_copy)

        if mol_copy.GetNumConformers() == 0:
            logger.warning("Failed to create conformer. Cannot calculate SASA.")
            return 0.0

        conf = mol_copy.GetConformer()

        if len(coordinates) != mol.GetNumAtoms():
            logger.warning(
                f"Coordinate count mismatch: {len(coordinates)} coordinates "
                f"for {mol.GetNumAtoms()} atoms. Cannot calculate SASA."
            )
            return 0.0

        for i, (symbol, x, y, z) in enumerate(coordinates):
            if i < mol.GetNumAtoms():
                conf.SetAtomPosition(i, (x, y, z))

        sasa_radii = {
            'H': 0.23,
            'C': 0.62,
            'N': 0.57,
            'O': 0.52,
            'F': 0.47,
            'P': 0.72,
            'S': 0.67,
            'Cl': 0.62,
            'Br': 0.72,
            'I': 0.82,
            'B': 0.62,
            'Si': 0.82,
            'Se': 0.72,
            'Te': 0.82
        }

        radii = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            radii.append(sasa_radii.get(symbol, 1.2))

        radii_array = np.array(radii, dtype=np.float64)

        try:
            sasa = rdFreeSASA.CalcSASA(mol_copy, radii_array, conf.GetId())
            return float(sasa)
        except Exception as e:
            logger.warning(f"Failed to calculate SASA: {e}")
            volume = data.get("molecular_volume", 0.0)
            if volume > 0:
                estimated_sasa = 4.0 * (volume ** (2.0/3.0))
                return estimated_sasa
            return 0.0


    @handle_x_molecule
    def xy_shadow(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate XY shadow area (projection area on XY plane).

        Args:
            mol: RDKit molecule

        Returns:
            XY shadow area in Å²
        """
        data = self._get_output(mol)
        coordinates = data.get("coordinates", [])

        if not coordinates:
            logger.warning("No coordinates found in ORCA output")
            return 0.0

        x_coords = [x for _, x, _, _ in coordinates]
        y_coords = [y for _, _, y, _ in coordinates]

        if not x_coords or not y_coords:
            return 0.0

        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)

        area = (x_max - x_min) * (y_max - y_min)
        return area


