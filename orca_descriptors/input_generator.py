"""ORCA input file generator for version 6.0.1."""

import re
from typing import Optional

from rdkit.Chem import AllChem, Mol, MolToXYZBlock


class ORCAInputGenerator:
    """Generate ORCA input files for version 6.0.1."""
    
    # List of supported semi-empirical methods
    SEMI_EMPIRICAL_METHODS = {
        "AM1", "PM3", "PM6", "PM7", "RM1", "MNDO", "MNDOD", "OM1", "OM2", "OM3"
    }
    
    def _is_semi_empirical(self, functional: str) -> bool:
        """Check if the functional is a semi-empirical method.
        
        Args:
            functional: Method name
            
        Returns:
            True if semi-empirical, False otherwise
        """
        return functional.upper() in self.SEMI_EMPIRICAL_METHODS
    
    def generate(
        self,
        mol: Mol,
        functional: str = "AM1",
        basis_set: str = "def2-SVP",
        method_type: str = "Opt",
        dispersion_correction: Optional[str] = None,
        solvation_model: Optional[str] = None,
        n_processors: int = 1,
        max_scf_cycles: int = 100,
        scf_convergence: float = 1e-6,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> str:
        """Generate ORCA input file content.
        
        Args:
            mol: RDKit molecule object
            functional: DFT functional or semi-empirical method (e.g., "PBE0", "AM1", "PM3", "PM6", "PM7")
            basis_set: Basis set (ignored for semi-empirical methods)
            method_type: Calculation type (Opt, SP, etc.)
            dispersion_correction: Dispersion correction (e.g., "D3BJ", ignored for semi-empirical)
            solvation_model: Solvation model (e.g., "COSMO(Water)")
            n_processors: Number of processors
            max_scf_cycles: Maximum SCF cycles
            scf_convergence: SCF convergence threshold
            charge: Molecular charge
            multiplicity: Spin multiplicity
            
        Returns:
            ORCA input file content as string
        """
        lines = []
        
        calc_line_parts = []
        is_semi_empirical = self._is_semi_empirical(functional)
        
        if method_type == "Opt":
            calc_line_parts.append("! Opt")
        elif method_type == "SP":
            calc_line_parts.append("! SP")
        else:
            calc_line_parts.append(f"! {method_type}")
        
        if is_semi_empirical:
            # For semi-empirical methods, just add the method name
            # Format: ! Opt AM1 or ! SP PM3
            calc_line_parts.append(functional.upper())
            calc_line_parts.append("Mayer")
        else:
            # For DFT methods, add functional and basis set
            calc_line_parts.append(f"{functional}")
            calc_line_parts.append(f"{basis_set}")
            
            if dispersion_correction:
                if dispersion_correction == "D3BJ":
                    calc_line_parts.append("D3BJ")
                else:
                    calc_line_parts.append(dispersion_correction)
            
            calc_line_parts.append("SlowConv")
            calc_line_parts.append("TightSCF")
            calc_line_parts.append("NMR")
            calc_line_parts.append("Mayer")
            calc_line_parts.append("NBO")
        
        lines.append(" ".join(calc_line_parts))
        lines.append("")
        
        lines.append("%pal")
        lines.append(f"  nprocs {n_processors}")
        lines.append("end")
        lines.append("")
        
        lines.append("%scf")
        lines.append(f"  MaxIter {max_scf_cycles}")
        lines.append(f"  ConvForced true")
        lines.append("end")
        lines.append("")
        
        if method_type == "Opt":
            lines.append("%geom")
            lines.append("  MaxIter 200")
            lines.append("  MaxStep 0.1")
            lines.append("end")
            lines.append("")
        
        if solvation_model:
            if "COSMO" in solvation_model.upper():
                solvent_match = re.search(r"\(([^)]+)\)", solvation_model)
                if solvent_match:
                    solvent = solvent_match.group(1)
                    lines.append("%cpcm")
                    lines.append("  smd true")
                    lines.append(f"  smdsolvent \"{solvent}\"")
                    lines.append("end")
                    lines.append("")
        
        lines.append("%output")
        lines.append("  PrintLevel 2")
        lines.append("  Print[P_Mulliken] 1")
        lines.append("  Print[P_AtCharges_M] 1")
        lines.append("  Print[P_Mayer] 1")
        if not is_semi_empirical:
            lines.append("  Print[P_NBO] 1")
        lines.append("end")
        lines.append("")
        lines.append("")
        
        lines.append(f"* xyz {charge} {multiplicity}")
        
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol)
        
        conf = mol.GetConformer()
        
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            lines.append(f"{symbol:2s}  {pos.x:15.10f}  {pos.y:15.10f}  {pos.z:15.10f}")
        
        lines.append("*")
        
        return "\n".join(lines)

