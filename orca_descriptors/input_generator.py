"""ORCA input file generator for version 6.0.1."""

import re
from typing import Optional

from rdkit.Chem import AllChem, Mol, MolToXYZBlock


class ORCAInputGenerator:
    """Generate ORCA input files for version 6.0.1."""
    
    def generate(
        self,
        mol: Mol,
        functional: str = "PBE0",
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
            functional: DFT functional
            basis_set: Basis set
            method_type: Calculation type (Opt, SP, etc.)
            dispersion_correction: Dispersion correction (e.g., "D3BJ")
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
        
        # Build calculation line
        calc_line_parts = []
        
        if method_type == "Opt":
            calc_line_parts.append("! Opt")
        elif method_type == "SP":
            calc_line_parts.append("! SP")
        else:
            calc_line_parts.append(f"! {method_type}")
        
        # Add functional and basis set
        calc_line_parts.append(f"{functional}")
        calc_line_parts.append(f"{basis_set}")
        
        # Add dispersion correction
        if dispersion_correction:
            if dispersion_correction == "D3BJ":
                calc_line_parts.append("D3BJ")
            else:
                calc_line_parts.append(dispersion_correction)
        
        # Add SCF settings
        calc_line_parts.append("SlowConv")
        
        # Add output options
        calc_line_parts.append("TightSCF")
        
        # Write calculation line
        lines.append(" ".join(calc_line_parts))
        lines.append("")
        
        # Write %pal block for parallelization
        lines.append("%pal")
        lines.append(f"  nprocs {n_processors}")
        lines.append("end")
        lines.append("")
        
        # Write %scf block
        # Note: In ORCA 6.0.1, Conv parameter is not valid in %scf block
        # Use ConvForced and convergence is controlled by TightSCF in command line
        lines.append("%scf")
        lines.append(f"  MaxIter {max_scf_cycles}")
        lines.append(f"  ConvForced true")
        # Conv parameter removed - not supported in ORCA 6.0.1
        # Convergence threshold is controlled by TightSCF keyword
        lines.append("end")
        lines.append("")
        
        # Write %geom block for geometry optimization settings
        # Limit the number of optimization cycles to prevent infinite loops
        if method_type == "Opt":
            lines.append("%geom")
            lines.append("  MaxIter 200")  # Maximum number of geometry optimization cycles
            lines.append("  MaxStep 0.1")   # Maximum step size (in Angstrom)
            lines.append("end")
            lines.append("")
        
        # Write %cpcm block for solvation (COSMO)
        if solvation_model:
            # Parse solvation model (e.g., "COSMO(Water)" -> solvent="Water")
            if "COSMO" in solvation_model.upper():
                # Extract solvent name from COSMO(Solvent)
                solvent_match = re.search(r"\(([^)]+)\)", solvation_model)
                if solvent_match:
                    solvent = solvent_match.group(1)
                    lines.append("%cpcm")
                    # In ORCA 6.0.1, CPCM uses smd true and smdsolvent
                    lines.append("  smd true")
                    lines.append(f"  smdsolvent \"{solvent}\"")
                    lines.append("end")
                    lines.append("")
        
        # Write %output block for detailed output
        # In ORCA 6.0.1, use minimal output settings to avoid errors
        lines.append("%output")
        lines.append("  PrintLevel 2")
        lines.append("  Print[P_Mulliken] 1")
        lines.append("  Print[P_AtCharges_M] 1")
        # Only include valid Print[] options for ORCA 6.0.1
        # P_Geo and P_EnGrad are not valid in ORCA 6.0.1
        lines.append("end")
        lines.append("")
        
        # Write coordinates block
        # For ORCA 6.0.1, use * xyzfile format or * xyz with charge and multiplicity on same line
        # Format: * xyzfile charge multiplicity filename
        # OR: * xyz charge multiplicity (inline coordinates)
        lines.append(f"* xyz {charge} {multiplicity}")
        
        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            # Add 3D coordinates using ETKDG
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
            if mol.GetNumConformers() == 0:
                # Fallback to basic embedding
                AllChem.EmbedMolecule(mol)
        
        # Get conformer
        conf = mol.GetConformer()
        
        # Generate XYZ coordinates manually for better control
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            # Format: "Symbol    x    y    z"
            lines.append(f"{symbol:2s}  {pos.x:15.10f}  {pos.y:15.10f}  {pos.z:15.10f}")
        
        lines.append("*")
        
        return "\n".join(lines)

