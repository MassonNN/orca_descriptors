"""Decorators for ORCA descriptor methods."""

import functools
from typing import Any

from rdkit.Chem import Mol


def handle_x_molecule(func):
    """Decorator to handle XMolecule for descriptor methods.
    
    If the first argument (mol) is an XMolecule instance, returns a DescriptorCall
    object instead of executing the function. Otherwise, calls the original function.
    
    The decorator correctly handles methods with required parameters after mol
    (e.g., mo_energy(mol, index), get_bond_lengths(mol, atom1, atom2)) by
    capturing all positional and keyword arguments.
    
    Args:
        func: The descriptor method to decorate
        
    Returns:
        Decorated function that handles XMolecule instances
    """
    @functools.wraps(func)
    def wrapper(self, mol: Mol, *args, **kwargs):
        try:
            from orca_descriptors.batch_processing import XMolecule, DescriptorCall
            if isinstance(mol, XMolecule):
                return DescriptorCall(func.__name__, args, kwargs)
        except ImportError:
            pass
        
        return func(self, mol, *args, **kwargs)
    
    return wrapper

