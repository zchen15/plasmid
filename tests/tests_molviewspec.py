import pytest
import os
import tempfile
from pathlib import Path
from plasmid.pdb.record import Record, read_pdb
from plasmid.pdb.molview import Viewer

# Path to the target CIF file
CIF_FILE_PATH = Path("data/pdb/4ZT9.cif")

@pytest.fixture
def cif_file_path() -> Path:
    """Provides the path to the 4ZT9 CIF file."""
    # Ensure the file exists before yielding, although we rely on fixture setup
    if not CIF_FILE_PATH.exists():
        raise FileNotFoundError(f"CIF file not found at {CIF_FILE_PATH}")
    yield CIF_FILE_PATH
    # No cleanup needed as the file is static data

def test_viewer_creation_with_cif_file(cif_file_path):
    """Tests that Viewer can be created successfully using a local CIF file."""
    # 1. Read CIF data into a Record object
    # Assuming read_pdb handles mmCIF (which is common for PDB tools)
    rec = read_pdb(cif_file_path)
    
    # 2. Setup a temporary working directory for the Viewer to use
    with tempfile.TemporaryDirectory() as temp_dir:
        # Instantiate Viewer
        viewer = Viewer(rec=rec)
        
        # 3. Assert that the viewer was successfully created (no exceptions raised)
        assert viewer is not None
        assert viewer.builder is not None
        assert viewer.structure is not None
