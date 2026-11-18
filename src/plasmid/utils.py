#!/usr/bin/env python
import importlib

def library_safeload(name: str):
    """
    Conditionally load a library by name.
    
    Args:
        name: Import name of the library (e.g., "gemmi", "biotite.structure").

    Returns:
        The imported module, or None if missing and required=False.
    """
    try:
        return importlib.import_module(name)
    except ImportError:
        print(f"failed to load {name}")
        return None
