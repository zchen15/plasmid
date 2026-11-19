from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import biotite.structure
import biotite.structure.io
import biotite.structure.io.pdbx

from .locations import Location

ENTITY_TYPE = {
    0:"ligand",
    1:"amino_acid",
    2:"nucleic_acid",
    3:"ion",
    4:"solvent"
}

def atom_array_to_df(
    atom_array: biotite.structure.AtomArray | biotite.structure.AtomArrayStack,
) -> pd.DataFrame:
    """
    Convert an AtomArray or AtomArrayStack to a pandas dataframe.
    
    Args:
        atom_array: biotite AtomArray or AtomArrayStack.
    
    Returns:
        A dataframe with one row per atom. For an AtomArrayStack the
        ``frame`` column identifies the model and ``atom_id`` is the index
        of the atom within its frame.
    """
    
    annotations = set(atom_array.get_annotation_categories())
    data: dict[str, Any] = {
        "atom_id": np.arange(len(atom_array)),
        "chain_id": atom_array.chain_id,
        "res_id": atom_array.res_id,
        "res_name": atom_array.res_name,
        "hetero": atom_array.hetero,
        "atom_name": atom_array.atom_name,
        "element": atom_array.element,
        "x": atom_array.coord[:, 0],
        "y": atom_array.coord[:, 1],
        "z": atom_array.coord[:, 2],
        "b_factor": (atom_array.get_annotation("b_factor")
                     if "b_factor" in annotations
                     else np.zeros(len(atom_array))),
        "occupancy": (atom_array.get_annotation("occupancy")
                      if "occupancy" in annotations
                      else np.ones(len(atom_array))),
        "charge": (atom_array.get_annotation("charge")
                   if "charge" in annotations
                   else np.zeros(len(atom_array))),
    }

    # preserve any other annotation attached to the AtomArray (e.g. ins_code)
    for key in annotations:
        if key not in data:
            data[key] = atom_array.get_annotation(key)

    # classify each atom for downstream filtering / display
    feature = np.zeros(len(atom_array))
    feature[biotite.structure.filter_solvent(atom_array)] = ENTITY_TYPE["solvent"]
    feature[biotite.structure.filter_monoatomic_ions(atom_array)] = ENTITY_TYPE["ions"]
    feature[biotite.structure.filter_amino_acids(atom_array)] = ENTITY_TYPE["amino_acid"]
    feature[biotite.structure.filter_nucleotides(atom_array)] = ENTITY_TYPE["nucleic_acid"]
    data["entity_type"] = feature
    
    df = pd.DataFrame(data)
    df["frame"] = 0
    return df

def add_bonds(record: "Record") -> "Record":
    """
    Add a bond list based on residue names and inter-atomic distances.

    Args:
        record: Record object.

    Returns:
        Copy of the record with ``AtomArray.bonds`` assigned.
    """
    df = atom_array_to_df(record.atom_array)
    df = df[df["frame"] == 0]
    c0 = biotite.structure.filter_amino_acids(record.atom_array)
    c0 |= biotite.structure.filter_nucleotides(record.atom_array)

    atom_list = []
    for chain_id in df["chain_id"].drop_duplicates():
        c1 = df["chain_id"] == chain_id
        c_poly = c0 & c1
        c_lig = ~c0 & c1
        if isinstance(record.atom_array, biotite.structure.AtomArrayStack):
            poly_atom_array = record.atom_array[:, c_poly]
            lig_atom_array = record.atom_array[:, c_lig]
        else:
            poly_atom_array = record.atom_array[c_poly]
            lig_atom_array = record.atom_array[c_lig]

        # bonds for residues
        bondlist = biotite.structure.connect_via_residue_names(poly_atom_array)
        poly_atom_array.bonds = bondlist
        atom_list.append(poly_atom_array)

        # bonds for ligands
        if isinstance(record.atom_array, biotite.structure.AtomArrayStack):
            bondlist = biotite.structure.connect_via_distances(lig_atom_array[0])
        else:
            bondlist = biotite.structure.connect_via_distances(lig_atom_array)
        lig_atom_array.bonds = bondlist
        atom_list.append(lig_atom_array)

    # concatenate atom_array
    for i in range(1, len(atom_list)):
        atom_list[0] = atom_list[0] + atom_list[i]

    # return copy
    out = record.copy()
    out.atom_array = atom_list[0]
    return out

def df_to_atom_stack(df):
    # build one AtomArray per frame and stack when there is more than one
    frames = [df_to_atom_array(frame_df) for _, frame_df in df.groupby("frame", sort=True)]
    return biotite.structure.stack(frames)

def df_to_atom_array(
    df: pd.DataFrame,
) -> biotite.structure.AtomArray | biotite.structure.AtomArrayStack:
    """
    Convert a dataframe produced by :func:`atom_array_to_df` to an AtomArray.

    Args:
        df: dataframe with the columns returned by :func:`atom_array_to_df`.

    Returns:
        A single AtomArray when the dataframe contains one frame, otherwise an
        AtomArrayStack with one model per ``frame`` value.
    """
    if "frame" in df.columns and df["frame"].nunique() > 1:
        return df_to_atom_stack(df)

    df = df.reset_index(drop=True)
    atom_array = biotite.structure.AtomArray(len(df))
    atom_array.coord = df[["x", "y", "z"]].to_numpy(dtype=float)
    atom_array.chain_id = df["chain_id"].to_numpy(dtype=str)
    atom_array.res_id = df["res_id"].to_numpy(dtype=int)
    atom_array.res_name = df["res_name"].to_numpy(dtype=str)
    atom_array.hetero = df["hetero"].to_numpy(dtype=bool)
    atom_array.atom_name = df["atom_name"].to_numpy(dtype=str)
    atom_array.element = df["element"].to_numpy(dtype=str)
    atom_array.set_annotation("b_factor", df["b_factor"].to_numpy(dtype=float))
    atom_array.set_annotation("occupancy", df["occupancy"].to_numpy(dtype=float))
    
    # check if charge is integer or float
    atom_array.set_annotation("charge", df["charge"].to_numpy(dtype=float))

    # restore every other non-derived column as an annotation
    skip = {"atom_id", "frame", "x", "y", "z", "feature",
            "chain_id", "res_id", "res_name", "hetero", "atom_name",
            "element", "b_factor", "occupancy", "charge"}
    for column in df.columns:
        if column in skip:
            continue
        values = df[column].to_numpy()
        if values.dtype == object:
            non_null = values[~pd.isna(values)]
            if all(isinstance(value, str) for value in non_null):
                values = values.astype(str)
        atom_array.set_annotation(column, values)
    return atom_array

def to_mmcif(filename:str, df_atom_array, df_pep=None, df_seq=None, df_nuc=None) -> None:
    """
    Save an AtomArray or AtomArrayStack to an mmCIF file.

    Args:
        atom_array: structure to write.
        filename: output mmCIF file name.
    """
    # accept either a dataframe or an AtomArray/AtomArrayStack
    if isinstance(df_atom_array, pd.DataFrame):
        df = df_atom_array
        df = df[df["frame"] == 0] if "frame" in df.columns else df
        atom_array = df_to_atom_array(df)
    else:
        atom_array = df_atom_array
        df = atom_array_to_df(atom_array)
        df = df[df["frame"] == 0]

    cols = ["chain_id", "res_id", "res_name"]

    # get list of peptides
    if df_pep is None:
        c0 = biotite.structure.filter_amino_acids(atom_array)
        df_pep = df.loc[c0, cols].drop_duplicates()
        df_pep = df_pep.sort_values(by=cols)

    # get list of nucleotides
    if df_nuc is None:
        c0 = biotite.structure.filter_nucleotides(atom_array)
        df_nuc = df.loc[c0, cols].drop_duplicates()

    # turn chain_id into integers
    peptide_list = df_pep["chain_id"].drop_duplicates().values
    for i in range(len(peptide_list)):
        c = df_pep["chain_id"] == peptide_list[i]
        df_pep.loc[c, "entity_id"] = i + 1
    df_pep["entity_type"] = "polymer"
    df_pep["poly_type"] = "polypeptide(L)"

    nucleotide_list = df_nuc["chain_id"].drop_duplicates().values
    for i in range(len(nucleotide_list)):
        c = df_nuc["chain_id"] == nucleotide_list[i]
        df_nuc.loc[c, "entity_id"] = i + len(peptide_list) + 1
    df_nuc["entity_type"] = "polymer"
    df_nuc["poly_type"] = "polyribonucleotide"
    x = pd.concat([df_pep, df_nuc])
    x["entity_id"] = x["entity_id"].astype(int)

    # write to entity_poly_seq fields
    entity_poly_seq = biotite.structure.io.pdbx.CIFCategory()
    entity_poly_seq["entity_id"] = list(x["entity_id"].values)
    entity_poly_seq["num"] = list(x["res_id"].values)
    entity_poly_seq["mon_id"] = list(x["res_name"].values)

    # write entity and entity_poly fields
    if df_seq is None:
        df_seq = sequence_df(atom_array)
    x = x[["chain_id", "entity_id", "entity_type", "poly_type"]].drop_duplicates()
    x = x.merge(df_seq[["chain_id", "sequence"]], on=["chain_id"], how="left")

    entity_poly = biotite.structure.io.pdbx.CIFCategory()
    entity_poly["entity_id"] = list(x["entity_id"].values)
    entity_poly["type"] = list(x["poly_type"].values)
    entity_poly["pdbx_seq_one_letter_code"] = list(x["sequence"].values)

    entity = biotite.structure.io.pdbx.CIFCategory()
    entity["id"] = list(x["entity_id"].values)
    entity["type"] = list(x["entity_type"].values)

    # write struct_asym which maps chain_id to entity_id
    struct_asym = biotite.structure.io.pdbx.CIFCategory()
    struct_asym["id"] = list(x["chain_id"].values)
    struct_asym["entity_id"] = list(x["entity_id"].values)

    block = {"entity": entity,
             "entity_poly": entity_poly,
             "struct_asym": struct_asym,
             "entity_poly_seq": entity_poly_seq}
    block = biotite.structure.io.pdbx.CIFBlock(block)
    file = biotite.structure.io.pdbx.CIFFile({"header": block})
    
    # add structure info
    biotite.structure.io.pdbx.set_structure(file, atom_array)
    file.write(str(filename))
