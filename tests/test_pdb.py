import pytest
import pandas as pd
import os
from plasmid.pdb.record import read_pdb, Record

# Define the test file path
TEST_CIF_FILENAME = "data/pdb/9C3E.cif"

def test_read_cif_success():
    # Test reading the CIF file
    record = read_pdb(TEST_CIF_FILENAME)
    assert isinstance(record, Record)
    assert record.filename == TEST_CIF_FILENAME
    assert record.atom_array is not None
    assert len(record) > 0 # Check if atoms were successfully loaded

def test_record_to_dataframe():
    # Test conversion to DataFrame
    record = read_pdb(TEST_CIF_FILENAME)
    df = record.to_dataframe()
    
    assert isinstance(df, pd.DataFrame)
    # Expect at least 2 atoms (based on the dummy PDB structure, though CIF might be larger)
    assert len(df) > 0 
    # Check for expected columns
    expected_cols = ["atom_id", "chain_id", "res_id", "res_name", "hetero", "atom_name", "element", "x", "y", "z", "frame"]
    for col in expected_cols:
        assert col in df.columns

def test_read_cif_file_not_found():
    # Test file not found scenario
    non_existent_file = "test/non_existent.cif"
    with pytest.raises(FileNotFoundError):
        read_pdb(non_existent_file)