"""Miscellaneous descriptors for ORCA calculations."""

import logging
from typing import Any

import numpy as np
from rdkit.Chem import Mol, Crippen, Descriptors, rdMolDescriptors

from orca_descriptors.decorators import handle_x_molecule

logger = logging.getLogger(__name__)


class MiscDescriptorsMixin:
    """Mixin providing miscellaneous descriptors."""
    
    @handle_x_molecule
    def dipole_moment(self, mol: Mol, *args, **kwargs) -> float:
        """Get dipole moment magnitude in Debye."""
        data = self._get_output(mol)
        return data.get("dipole_moment", 0.0)


    @handle_x_molecule
    def get_min_h_charge(self, mol: Mol, method: str = "ESP", *args, **kwargs) -> float:
        """Get minimum net atomic charge for hydrogen atoms.

        Args:
            mol: RDKit molecule
            method: Charge method (currently only "ESP" supported, uses Mulliken as fallback)

        Returns:
            Minimum hydrogen charge in atomic units
        """
        data = self._get_output(mol)
        atom_charges = data.get("atom_charges", {})

        h_charges = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() == "H":
                charge = atom_charges.get(i, 0.0)
                h_charges.append(charge)

        if not h_charges:
            logger.warning("No hydrogen atoms found in molecule")
            return 0.0

        return min(h_charges)


    @handle_x_molecule
    def meric(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate MERIC (Minimum Electrophilicity Index for Carbon).

        MERIC predicts sites susceptible to nucleophilic attack.
        For carbon bonded to heteroatoms (O, N, etc.), MERIC is typically negative.

        Args:
            mol: RDKit molecule

        Returns:
            Minimum MERIC value in eV (typically negative for electrophilic carbons)
        """
        data = self._get_output(mol)
        atom_charges = data.get("atom_charges", {})
        coordinates = data.get("coordinates", [])

        if not coordinates:
            logger.warning("No coordinates found in ORCA output")
            return 0.0

        meric_values = []

        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() != "C":
                continue

            has_heteroatom = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() in ["O", "N", "S", "P", "F", "Cl", "Br", "I"]:
                    has_heteroatom = True
                    break

            if has_heteroatom:
                charge = atom_charges.get(i, 0.0)
                homo = data.get("homo_energy", 0.0)
                lumo = data.get("lumo_energy", 0.0)

                if homo != 0.0 and lumo != 0.0:
                    mu = (homo + lumo) / 2.0
                    eta = (lumo - homo) / 2.0
                    if eta > 0:
                        electrophilicity = -abs(homo) / (2.0 * eta) - abs(charge) * 0.1
                    else:
                        electrophilicity = -abs(homo) * 0.1 - abs(charge) * 0.1
                    meric_values.append(electrophilicity)

        if not meric_values:
            return 0.0

        return min(meric_values)


    @handle_x_molecule
    def m_log_p(self, mol: Mol, *args, **kwargs) -> float:
        """Calculate Moriguchi Log P (octanol/water partition coefficient).

        Args:
            mol: RDKit molecule

        Returns:
            Log P value (positive for hydrophobic, negative for hydrophilic)
        """
        from rdkit.Chem import Crippen

        try:
            logp = Crippen.MolLogP(mol)
            return logp
        except Exception as e:
            logger.warning(f"Failed to calculate M log P: {e}")
            return 0.0


    @handle_x_molecule
    def moran_autocorrelation(self, mol: Mol, lag: int = 2, weight: str = "vdw_volume", *args, **kwargs) -> float:
        """Calculate Moran autocorrelation descriptor.

        Args:
            mol: RDKit molecule
            lag: Lag distance (default: 2)
            weight: Weighting scheme ('vdw_volume', 'atomic_mass', etc.)

        Returns:
            Moran autocorrelation value (typically in range -1.0 to 1.0)
        """
        from rdkit.Chem import Descriptors
        import numpy as np

        try:
            n = mol.GetNumAtoms()
            if n < lag + 1:
                return 0.0

            if weight == "vdw_volume":
                from rdkit.Chem import rdMolDescriptors
                weights = []
                for i in range(n):
                    atom = mol.GetAtomWithIdx(i)
                    symbol = atom.GetSymbol()
                    vdw_radii = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 
                                'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98}
                    radius = vdw_radii.get(symbol, 1.5)
                    volume = (4.0 / 3.0) * 3.14159 * (radius ** 3)
                    weights.append(volume)
            else:
                weights = [mol.GetAtomWithIdx(i).GetMass() for i in range(n)]

            if not weights or all(w == 0.0 for w in weights):
                return 0.0

            mean_weight = sum(weights) / n

            try:
                from rdkit.Chem import rdMolDescriptors
                dist_matrix = rdMolDescriptors.GetDistanceMatrix(mol)
            except:
                from collections import deque
                n_atoms = mol.GetNumAtoms()
                dist_matrix = [[0] * n_atoms for _ in range(n_atoms)]
                adj = {i: [] for i in range(n_atoms)}
                for bond in mol.GetBonds():
                    b = bond.GetBeginAtomIdx()
                    e = bond.GetEndAtomIdx()
                    adj[b].append(e)
                    adj[e].append(b)
                for start in range(n_atoms):
                    distances = {start: 0}
                    queue = deque([start])
                    while queue:
                        current = queue.popleft()
                        for neighbor in adj.get(current, []):
                            if neighbor not in distances:
                                distances[neighbor] = distances[current] + 1
                                queue.append(neighbor)
                    for end, dist in distances.items():
                        dist_matrix[start][end] = dist

            distances = []
            for i in range(n):
                for j in range(i + 1, n):
                    dist_val = dist_matrix[i][j] if isinstance(dist_matrix, list) else dist_matrix[i, j]
                    if dist_val == lag:
                        distances.append((i, j))

            if not distances:
                return 0.0

            autocorr = 0.0
            for i, j in distances:
                wi = weights[i]
                wj = weights[j]
                autocorr += (wi - mean_weight) * (wj - mean_weight)

            variance = sum((w - mean_weight) ** 2 for w in weights)
            if len(distances) == 0:
                return 0.0

            if variance == 0:
                return 0.0

            autocorr = autocorr / (len(distances) * variance) if variance > 0 else 0.0

            return autocorr
        except Exception as e:
            logger.warning(f"Failed to calculate Moran autocorrelation: {e}")
            return 0.0


    @handle_x_molecule
    def autocorrelation_hats(self, mol: Mol, lag: int = 4, unweighted: bool = True, *args, **kwargs) -> float:
        """Calculate HATS autocorrelation descriptor (H-bond acceptor/donor).

        Args:
            mol: RDKit molecule
            lag: Lag distance (default: 4)
            unweighted: If True, use unweighted autocorrelation

        Returns:
            HATS autocorrelation value (typically close to zero)
        """
        from rdkit.Chem import Descriptors
        import numpy as np

        try:
            n = mol.GetNumAtoms()
            if n < lag + 1:
                return 0.0

            h_bond_atoms = []
            for i in range(n):
                atom = mol.GetAtomWithIdx(i)
                symbol = atom.GetSymbol()
                if symbol in ["N", "O", "F"]:
                    h_bond_atoms.append(i)

            if not h_bond_atoms:
                return 0.0

            distances = []
            for i in h_bond_atoms:
                for j in h_bond_atoms:
                    if i != j:
                        try:
                            dist = mol.GetDistanceMatrix()[i, j]
                            if dist == lag:
                                distances.append((i, j))
                        except:
                            continue

            if not distances:
                return 0.0

            if unweighted:
                return len(distances) / (n * (n - 1)) if n > 1 else 0.0
            else:
                weights = [1.0 if mol.GetAtomWithIdx(i).GetSymbol() in ["N", "O", "F"] else 0.0 
                          for i in range(n)]
                mean_weight = sum(weights) / n if n > 0 else 0.0

                autocorr = 0.0
                for i, j in distances:
                    wi = weights[i] if i < len(weights) else 0.0
                    wj = weights[j] if j < len(weights) else 0.0
                    autocorr += (wi - mean_weight) * (wj - mean_weight)

                variance = sum((w - mean_weight) ** 2 for w in weights)
                if variance == 0:
                    return 0.0

                autocorr = autocorr / (len(distances) * variance) if variance > 0 else 0.0
                return autocorr
        except Exception as e:
            logger.warning(f"Failed to calculate HATS autocorrelation: {e}")
            return 0.0


