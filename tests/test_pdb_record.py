import pytest
import numpy as np
import pandas as pd
import biotite.structure.io.pdbqt as pdbqt

from plasmid.pdb.record import read_pdb
from plasmid.pdb.io import atom_array_to_df, add_bonds

import os
import glob

def test_read_pdb():
    # Test conversion to DataFrame
    files = glob.glob("data/pdb/*.cif")
    for fname in files:
        rec = read_pdb(fname)

        print(f"atom_array in {fname}: \n{rec.atom_array}")
        print(f"__repr__() of {fname}: \n{rec}")
        
        # check columns
        expected_cols = ["atom_id", "chain_id", "res_id", "res_name", 
                         "hetero", "atom_name", "element", "entity_type",
                         "x", "y", "z",
                         "b_factor", "occupancy", "charge",
                         "frame"]
        for col in expected_cols:
            assert col in rec.atom_array

def expected_locations(rec, chain_id):
    """Return the atom and residue positions for a chain."""
    c = rec.atom_array["chain_id"] == chain_id
    atoms = rec.atom_array.loc[c, "atom_id"].values
    res_ids = np.unique(rec.atom_array.loc[c, "res_id"].values)
    return atoms, res_ids

def test_annotate_gRNA():
    fname = "data/pdb/4ZT9.cif"
    rec = read_pdb(fname)

    # non inplace annotation of gRNA
    c = rec["chain_id"] == "B"
    rec1 = rec.annotate_by_atom_id(atom_id = rec[c, "atom_id"],
                            name = "gRNA B",
                            color = "blue",
                           )
    c = rec["chain_id"] == "D"
    rec2 = rec.annotate_by_atom_id(atom_id = rec[c, "atom_id"],
                            name = "gRNA D",
                            color = "red",
                           )    
    print(f"orig dCas9 contents:\n{rec}")
    print(f"rec1 dCas9 contents:\n{rec1}")
    print(f"rec2 dCas9 contents:\n{rec2}")
    
    # expected atom/residue locations for each gRNA chain
    expected_atoms_b, expected_res_b = expected_locations(rec, "B")
    expected_atoms_d, expected_res_d = expected_locations(rec, "D")
    
    assert len(rec.features) == 0
    assert rec1.features["name"].values[0] == "gRNA B"
    assert rec2.features["name"].values[0] == "gRNA D"
    
    # check the non-inplace annotations keep their locations
    feature_b = rec1.features.iloc[0]
    assert feature_b["chain_id"] == "B"
    assert feature_b["entity_type"] == "nucleic_acid"
    assert feature_b["color"] == "blue"
    assert feature_b["representation"] == "cartoon"
    assert np.array_equal(feature_b["atom_id"].positions, expected_atoms_b)
    assert np.array_equal(feature_b["res_id"].positions, expected_res_b)
    
    feature_d = rec2.features.iloc[0]
    assert feature_d["chain_id"] == "D"
    assert feature_d["entity_type"] == "nucleic_acid"
    assert feature_d["color"] == "red"
    assert feature_d["representation"] == "cartoon"
    assert np.array_equal(feature_d["atom_id"].positions, expected_atoms_d)
    assert np.array_equal(feature_d["res_id"].positions, expected_res_d)
    
    # annotate gRNAs inplace
    c = rec["chain_id"] == "B"
    rec.annotate_by_atom_id(atom_id = rec[c, "atom_id"],
                            name = "gRNA B",
                            color = "blue",
                            inplace = True,
                           )
    c = rec["chain_id"] == "D"
    rec.annotate_by_atom_id(atom_id = rec[c, "atom_id"],
                            name = "gRNA D",
                            color = "red",
                            inplace = True,
                           )
    
    print(f"dCas9 contents:\n{rec}")
    
    # check the inplace annotations keep their locations
    assert len(rec.features) == 2
    expected = {
        "gRNA B": ("B", "blue", expected_atoms_b, expected_res_b),
        "gRNA D": ("D", "red", expected_atoms_d, expected_res_d),
    }
    for _, feature in rec.features.iterrows():
        chain_id, color, atoms, res_ids = expected[feature["name"]]
        assert feature["chain_id"] == chain_id
        assert feature["entity_type"] == "nucleic_acid"
        assert feature["color"] == color
        assert feature["representation"] == "cartoon"
        assert np.array_equal(feature["atom_id"].positions, atoms)
        assert np.array_equal(feature["res_id"].positions, res_ids)

def test_annotate_add():
    fname = "data/pdb/4ZT9.cif"
    rec = read_pdb(fname)
    
    c1 = rec["chain_id"] == "B"
    rec.annotate_by_atom_id(atom_id = rec[c1, "atom_id"],
                            name = "gRNA",
                            color = "blue",
                            inplace = True,
                           )
    c2 = rec["chain_id"] == "A"
    rec.annotate_by_atom_id(atom_id = rec[c2, "atom_id"],
                            name = "dCas9",
                            color = "orange",
                            inplace = True,
                           )
    print(f"annotated record:\n{rec}")

    # filter for chain A and B
    rec = rec[c1 | c2]
    
    # add records
    x = rec.copy() + rec.copy()
    assert len(x.atom_array) > 0
    assert len(x.atom_array) == len(rec.atom_array)*2
    assert len(x.features) == 4

def test_reset_index():
    fname = "data/pdb/4ZT9.cif"
    rec = read_pdb(fname)
    rec.atom_array["atom_id"] = np.random.permutation(rec.atom_array["atom_id"].values)
    print(rec)
    x = rec.reset_index()
    print(x)

    # check at atom positions are contiguous
    df = x.get_amino_acids()
    assert all(len(loc.get_slices())==1 for loc in df["atom_id"])
    
    # check at atom positions are contiguous
    df = x.get_nucleic_acids()
    assert all(len(loc.get_slices())==1 for loc in df["atom_id"])
    
def test_split_and_merge_structures(tmp_path):
    # Structure used for the split/merge round-trip test
    TEST_COMPLEX_FILENAME = "data/pdb/4ZT9.cif"
    
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

    # atom_array is already a dataframe in the current API; compare the
    # structures after normalising row order
    sort_cols = ["chain_id", "res_id", "atom_name", "element", "res_name", "hetero"]
    recombined_df = (recombined.atom_array.drop(columns=["atom_id"])
                     .sort_values(sort_cols).reset_index(drop=True))
    expected_df = (expected.atom_array.drop(columns=["atom_id"])
                   .sort_values(sort_cols).reset_index(drop=True))
    pd.testing.assert_frame_equal(recombined_df, expected_df)

