"""Energy-related descriptors for ORCA calculations."""

from rdkit.Chem import Mol

from orca_descriptors.decorators import handle_x_molecule


class EnergyDescriptorsMixin:
    """Mixin providing energy-related descriptors."""
    
    @handle_x_molecule
    def total_energy(self, mol: Mol, *args, **kwargs) -> float:
        """Get total energy in Hartree."""
        data = self._get_output(mol)
        return data.get("total_energy", 0.0)
    
    @handle_x_molecule
    def gibbs_free_energy(self, mol: Mol, *args, **kwargs) -> float:
        """Get Gibbs free energy in Hartree."""
        data = self._get_output(mol)
        return data.get("gibbs_free_energy", data.get("total_energy", 0.0))
    
    @handle_x_molecule
    def entropy(self, mol: Mol, *args, **kwargs) -> float:
        """Get entropy in J/(mol·K)."""
        data = self._get_output(mol)
        return data.get("entropy", 0.0)
    
    @handle_x_molecule
    def enthalpy(self, mol: Mol, *args, **kwargs) -> float:
        """Get enthalpy in Hartree."""
        data = self._get_output(mol)
        return data.get("total_energy", 0.0)

