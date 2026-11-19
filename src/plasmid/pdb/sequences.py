"""
Sequence and secondary-structure annotation helpers for pdb/mmCIF records.

These functions operate on a :class:`plasmid.pdb.record.Record` and cache
their results on the record:

- ``get_sequence()`` populates ``record.sequence_annotations`` with
  ``atom_id``, ``chain_id``, ``res_id``, ``sequence`` and ``n_frames``.
- ``get_ligands()`` populates ``record.ligand_annotations`` with
  ``atom_id``, ``chain_id``, ``res_name``, ``res_id``, ``count`` and
  ``n_frames``.
- ``get_protein_sse()`` returns the sequence dataframe with an added
  ``dssp`` column.
- ``get_dssp_annotations()`` returns one row per secondary-structure segment
  with ``atom_id``, ``chain_id``, ``res_id``, ``dssp`` and ``length``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

import biotite.structure

from .io import atom_array_to_df, sequence_df
from .locations import Location

if TYPE_CHECKING:
    from .record import Record

# biotite 2.0 returns ``PSeaElement`` integer codes; older versions returned
# the DSSP letter codes directly. ``t`` (turn) is not produced by biotite's
# P-SEA implementation, but is accepted for backwards compatibility.
_SSE_CODES = {-1: "c", 0: "c", 1: "a", 2: "b"}


def get_sequence(record: "Record") -> pd.DataFrame:
    """
    Get the sequence of every amino-acid chain.

    The result is cached as ``record.sequence_annotations``.

    Args:
        record: Record object.

    Returns:
        Dataframe with the columns ``atom_id``, ``chain_id``, ``res_id``,
        ``sequence`` and ``n_frames``.
    """
    df = sequence_df(record.atom_array)

    # add info about the number of poses or frames
    if isinstance(record.atom_array, biotite.structure.AtomArrayStack):
        df["n_frames"] = len(record.atom_array)
    else:
        df["n_frames"] = 0

    record.sequence_annotations = df
    return df


def get_ligands(record: "Record") -> pd.DataFrame:
    """
    Get the list of non-polymer ligands in the structure.

    The result is cached as ``record.ligand_annotations``.

    Args:
        record: Record object.

    Returns:
        Dataframe with the columns ``atom_id``, ``chain_id``, ``res_name``,
        ``res_id``, ``count`` and ``n_frames``.
    """
    data = atom_array_to_df(record.atom_array)
    data = data[data["frame"] == 0]
    is_polymer = biotite.structure.filter_amino_acids(record.atom_array)
    is_polymer |= biotite.structure.filter_nucleotides(record.atom_array)
    data = data[~is_polymer]

    out = []
    for chain_id, res_name in data[["chain_id", "res_name"]].drop_duplicates().values:
        in_ligand = data["chain_id"] == chain_id
        in_ligand &= data["res_name"] == res_name
        atom_id = Location(data.loc[in_ligand, "atom_id"].values)
        res_id = Location(data.loc[in_ligand, "res_id"].values)
        out.append([atom_id, chain_id, res_name, res_id, len(res_id)])

    cols = ["atom_id", "chain_id", "res_name", "res_id", "count"]
    df = pd.DataFrame(out, columns=cols)

    # add info about number of poses or frames
    if isinstance(record.atom_array, biotite.structure.AtomArrayStack):
        df["n_frames"] = len(record.atom_array)
    else:
        df["n_frames"] = 0

    record.ligand_annotations = df
    return df


def get_protein_sse(record: "Record") -> pd.DataFrame:
    """
    Compute the secondary structure of each protein chain using P-SEA.

    Args:
        record: Record object.

    Returns:
        The sequence dataframe with an added ``dssp`` column containing one
        letter per residue.
    """
    df = get_sequence(record)
    out = []
    for atom_id in df["atom_id"].values:
        s1 = min(atom_id.positions)
        s2 = max(atom_id.positions)
        sse = biotite.structure.annotate_sse(record.atom_array[s1:s2])

        # biotite 2.0 uses integer PSeaElement codes, older versions letters
        try:
            sse = "".join(
                _SSE_CODES.get(int(code), "c") for code in sse
            )
        except (TypeError, ValueError):
            sse = "".join(sse)
        out.append(sse)
    df = df.copy()
    df["dssp"] = out
    return df


def get_dssp_annotations(record: "Record") -> pd.DataFrame:
    """
    Get secondary-structure annotations for the protein sequence.

    Args:
        record: Record object.

    Returns:
        Dataframe with one row per secondary-structure segment and the columns
        ``atom_id``, ``chain_id``, ``res_id``, ``dssp`` and ``length``.
    """
    data = atom_array_to_df(record.atom_array)
    data = data[data["frame"] == 0]
    df = get_protein_sse(record)

    out = []
    cols = ["chain_id", "atom_id", "dssp"]
    for chain_id, atom_id, sse in df[cols].values:
        # residues in this chain and within the annotated atom range
        in_chain = data["atom_id"] >= min(atom_id.positions)
        in_chain &= data["atom_id"] < max(atom_id.positions)
        in_chain &= data["chain_id"] == chain_id
        res_ids = data.loc[in_chain, "res_id"].drop_duplicates().values

        sse = np.array(list(sse))
        same_as_next = sse[:-1] == sse[1:]
        segment_start = np.concatenate(([True], ~same_as_next))
        segment_end = np.concatenate((~same_as_next, [True]))

        for start, end, code in zip(
            res_ids[segment_start],
            res_ids[segment_end],
            sse[segment_start],
        ):
            if code == "c":
                continue
            in_segment = data["chain_id"] == chain_id
            in_segment &= data["res_id"] >= start
            in_segment &= data["res_id"] <= end
            atom_ids = Location(data.loc[in_segment, "atom_id"].values)
            res_id = Location(data.loc[in_segment, "res_id"].values)
            out.append([atom_ids, chain_id, res_id, code, end - start])

    cols = ["atom_id", "chain_id", "res_id", "dssp", "length"]
    return pd.DataFrame(out, columns=cols)
