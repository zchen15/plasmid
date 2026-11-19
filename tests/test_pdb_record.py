import pytest
import numpy as np
import pandas as pd
import os
import biotite.structure.io.pdbqt as pdbqt
from plasmid.pdb.record import read_pdb, Record
from plasmid.pdb.io import atom_array_to_df, add_bonds

# Define the test file path
TEST_CIF_FILENAME = "data/pdb/9C3E.cif"

# Structure used for the split/merge round-trip test
TEST_COMPLEX_FILENAME = "data/pdb/4ZT9.cif"

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
    df = atom_array_to_df(record.atom_array)
    
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


def test_record_repr(capsys):
    # Test that __repr__ prints the sequence and ligand tables
    record = read_pdb(TEST_CIF_FILENAME)
    print(record)

    captured = capsys.readouterr()
    assert "List of Sequences:" in captured.out
    assert "List of Ligands:" in captured.out


def test_split_and_merge_structures(tmp_path):
    # Read the full dCas9/gRNA complex
    record = read_pdb(TEST_COMPLEX_FILENAME)

    # Split the complex into its gRNA (nucleotide chains) and
    # dCas9 protein (amino-acid chains)
    is_rna = np.isin(record.atom_array.chain_id, ["B", "D"])
    is_protein = np.isin(record.atom_array.chain_id, ["A", "C"])
    gRNA = record[is_rna]
    dCas9 = record[is_protein]

    # Write each component to the temporary directory
    gRNA_file = tmp_path / "gRNA.cif"
    dCas9_cif_file = tmp_path / "dCas9.cif"
    gRNA.save(str(gRNA_file))
    dCas9.save(str(dCas9_cif_file))

    # Write dCas9 to PDB format
    dCas9_pdb_file = tmp_path / "dCas9.pdb"
    dCas9.save(str(dCas9_pdb_file))

    # Write dCas9 to PDBQT format (requires a bond list)
    dCas9_pdbqt_file = tmp_path / "dCas9.pdbqt"
    dCas9_with_bonds = add_bonds(dCas9)
    pdbqt_file = pdbqt.PDBQTFile()
    pdbqt_file.set_structure(dCas9_with_bonds.atom_array)
    pdbqt_file.write(str(dCas9_pdbqt_file))

    for path in [gRNA_file, dCas9_cif_file, dCas9_pdb_file, dCas9_pdbqt_file]:
        assert path.exists()

    # Read the two CIF files back and recombine them
    recombined = (read_pdb(str(gRNA_file)) + read_pdb(str(dCas9_cif_file))).reset_index()
    expected = record.reset_index()

    pd.testing.assert_frame_equal(
        atom_array_to_df(recombined.atom_array).drop(columns=["atom_id"]),
        atom_array_to_df(expected.atom_array).drop(columns=["atom_id"]),
    )

