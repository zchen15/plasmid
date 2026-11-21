# Module for browsing and editing pdb and mmCIF files

# science libs
import json

import numpy as np
import pandas as pd
import biotite.structure.io

# custom libraries
from .locations import Location
from .io import atom_array_to_df, df_to_atom_array, df_to_atom_stack, ENTITY_TYPE

# canonical columns of io.atom_array_to_df(); ``frame`` is added by read_pdb()
ATOM_ARRAY_COLUMNS = [
    "atom_id", "chain_id", "res_id", "res_name", "hetero", "atom_name",
    "element", "x", "y", "z", "b_factor", "occupancy", "charge",
    "entity_type", "frame",
]
# canonical columns created by annotate_by_atom_id()
FEATURE_COLUMNS = [
    "atom_id", "chain_id", "res_id", "entity_type", "name",
    "representation", "color", "description",
]


def _trim_features(atom_array: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """Trim the feature table so that it only references atoms present in ``atom_array``.

    Features whose atoms were all removed are dropped; features which only partially
    overlap are trimmed to the surviving atoms.

    Args:
        atom_array: Atom dataframe to trim the features against.
        features: Feature table to trim.

    Returns:
        Copy of the feature table trimmed to ``atom_array``.
    """
    if len(features) == 0:
        return features

    atom_list = atom_array["atom_id"].drop_duplicates()
    features = features.copy()
    
    # keep only features that still have atoms
    features["atom_id"] = [loc.intersection(atom_list) for loc in features["atom_id"].values]
    c = [len(loc) > 0 for loc in features["atom_id"]]
    features = features[c].copy()
    
    # update res_id
    new_res_id = []
    for i in range(len(features)):
        loc = features.iloc[i]["atom_id"]
        c = atom_array["atom_id"].isin(loc.get_positions())
        new_res_id.append(Location(atom_array.loc[c, "res_id"]))
    features["res_id"] = new_res_id
    return features

def read_pdb(filename:str):
    """Reads a PDB or mmCIF file into a Record object.
    
    Args:
        filename: Path to the PDB or mmCIF file to read.
    
    Returns:
        A Record object containing the parsed structure.
    """
    atom_array = biotite.structure.io.load_structure(filename)
    if isinstance(atom_array, biotite.structure.AtomArrayStack):
        # todo vectorize better by extraccting all of biotite atom_array_stackz
        out = []
        for frame in range(len(atom_array)):
            df = atom_array_to_df(atom_array[frame])
            df["frame"] = frame
            out.append(df)
        df_atom_array = pd.concat(out)
    else:
        df_atom_array = atom_array_to_df(atom_array)
        df_atom_array["frame"] = 0
    
    rec = Record(filename = filename,
                 atom_array = df_atom_array,
                )
    return rec

class Record:
    """
    This class holds functions pertaining to modifying the atoms in a pdb file
    """
    filename = None
    atom_array = None
    features = None
    
    def __init__(self, **kwargs):
        """Initializes a Record object with default or provided attributes.
        
        Args:
            kwargs: Keyword arguments matchingRecord attributes to override
                defaults. Supported keys include ``filename``, ``atom_array``,
                and ``features``.
        
        Returns:
            None. The Record instance is modified in place.
        """
        self.atom_array = pd.DataFrame(columns=ATOM_ARRAY_COLUMNS)
        self.features = pd.DataFrame(columns=FEATURE_COLUMNS)
        
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        
    def copy(self):
        """Returns a deep copy of the Record instance.
        
        Returns:
            A new Record object with copied atom_array and features DataFrames.
        """
        rec = Record(filename = self.filename,
                     atom_array = self.atom_array.copy(),
                     features = self.features.copy(),
                    )
        return rec
        
    def reset_index(self):
        """Sorts the Record's atom array and feature table in place.
        
        Reassigns sequential atom IDs and re-sorts features by
        (start, length, feature, name) for consistent ordering.
        
        Returns:
            None. The Record is modified in place.
        """
        # sort and get new atom_id
        cols = ["atom_id","chain_id","res_id","res_name","hetero","atom_name","element"]
        new_id = self.atom_array[cols].drop_duplicates().sort_values(by=cols[1:])["atom_id"].values
        new_id = {new_id[i]:i for i in range(len(new_id))}
        
        atom_array = self.atom_array.copy()
        atom_array["atom_id"] = [new_id[i] for i in atom_array["atom_id"]]
        
        # remap atom_id of features
        features = self.features.copy()
        features["atom_id"] = [loc.remap(new_id) for loc in features["atom_id"]]
        stats = {"res_id":"sum", "atom_id":"sum"}
        cols = ["chain_id", "entity_type", "name", "representation", "color", "description"]
        features = features.groupby(cols).agg(stats).reset_index()

        # output
        out = Record(filename = self.filename,
                     atom_array = atom_array,
                     features = features,
                    )
        return out
        
    def save(self, filename:str, **kwargs):
        """Saves the Record's atom array to a PDB/mmCIF file.
        
        Args:
            filename: Path to the output file. Extension determines format
                (``.pdb`` or ``.cif``).
            **kwargs: Additional arguments passed to
                :func:`biotite.structure.io.save_structure`.
        
        Returns:
            None. The file is written to disk.
        """
        atom_array = df_to_atom_stack(self.atom_array)
        biotite.structure.io.save_structure(filename, atom_array, **kwargs)
        
    def get_ligands(self, frame:int=0):
        """Returns ligand and small-molecule atoms in the structure.
        
        Args:
            frame: Model index to select atoms from (default 0).
        
        Returns:
            DataFrame with columns ``[atom_id, chain_id, res_id, entity_type, res_name]``
            containing ligand/non-polymer atoms.
        """
        c1 = self.atom_array["entity_type"] == ENTITY_TYPE["nucleic_acid"]
        c1|= self.atom_array["entity_type"] == ENTITY_TYPE["amino_acid"]
        c2 = self.atom_array["frame"] == frame
        df = self.atom_array[~c1 & c2]
        
        REV_ENTITY_TYPE = {value:key for key, value in ENTITY_TYPE.items()}
        
        out = []
        cols = ["chain_id","res_name","entity_type"]
        for chain_id, res_name, etype in df[cols].drop_duplicates().values:
            c = df["chain_id"] == chain_id
            c&= df["res_name"] == res_name
            c&= df["entity_type"] == etype
            
            # get positions
            atom_id = Location(df.loc[c, "atom_id"].values)
            res_id = Location(df.loc[c, "res_id"].values)
            
            ename = REV_ENTITY_TYPE[etype]
            out.append([atom_id, chain_id, res_id, ename, res_name])
        out = pd.DataFrame(out, columns=["atom_id", "chain_id", "res_id", "entity_type", "res_name"])
        return out
        
    def get_amino_acids(self, frame:int=0):
        """Returns amino acid residue atoms in the structure.
        
        Args:
            frame: Model index to select atoms from (default 0).
        
        Returns:
            DataFrame with columns ``[atom_id, chain_id, sequence, res_id, entity_type]``
            containing per-residue amino acid information.
        """
        c1 = self.atom_array["entity_type"] == ENTITY_TYPE["amino_acid"]
        c2 = self.atom_array["frame"] == frame
        df = self.atom_array[c1 & c2]
        
        out = []
        for chain_id in df["chain_id"].drop_duplicates().values:
            c = df["chain_id"] == chain_id
            
            # get positions
            atom_id = Location(df.loc[c, "atom_id"].values)
            res_id = Location(df.loc[c, "res_id"].values)
            
            # convert three-letter residue names to the one-letter alphabet
            residues = df.loc[c, ["res_name", "res_id", "hetero"]].drop_duplicates()
            seqlist = ["-"]*(residues["res_id"].max() + 1)
            for ridx, res_name in residues[["res_id","res_name"]].values:
                try:
                    seqlist[ridx] = biotite.sequence.ProteinSequence.convert_letter_3to1(res_name)
                except Exception:
                    seqlist[ridx] = "-"

            out.append([atom_id, chain_id, "".join(seqlist), res_id, "amino_acid"])
        out = pd.DataFrame(out, columns=["atom_id", "chain_id", "sequence", "res_id", "entity_type"])
        return out
        
    def get_nucleic_acids(self, frame:int=0):
        """Returns nucleic acid residue atoms in the structure.
        
        Args:
            frame: Model index to select atoms from (default 0).
        
        Returns:
            DataFrame with columns ``[atom_id, chain_id, sequence, res_id, entity_type]``
            containing per-residue nucleic acid information.
        """
        c1= self.atom_array["entity_type"] == ENTITY_TYPE["nucleic_acid"]
        c2 = self.atom_array["frame"] == frame
        df = self.atom_array[c1 & c2]
        ename = "nucleic_acid"
        etype = ENTITY_TYPE[ename]
        
        out = []
        for chain_id in df["chain_id"].drop_duplicates().values:
            c = df["chain_id"] == chain_id
            
            # get positions
            atom_id = Location(df.loc[c, "atom_id"].values)
            res_id = Location(df.loc[c, "res_id"].values)
            
            # convert three-letter residue names to the one-letter alphabet
            residues = df.loc[c, ["res_name", "res_id", "hetero"]].drop_duplicates()
            seqlist = ["-"]*(residues["res_id"].max() + 1)
            for ridx, res_name in residues[["res_id","res_name"]].values:
                seqlist[ridx] = res_name
            
            out.append([atom_id, chain_id, "".join(seqlist), res_id, "nucleic_acid"])
        out = pd.DataFrame(out, columns=["atom_id", "chain_id", "sequence", "res_id", "entity_type"])
        return out
        
    def annotate_by_atom_id(self, atom_id:list, name:str, representation:str=None, color:str=None, description:dict=None, inplace:bool=False):
        """Adds annotation features to the PDB structure at specified atom positions.
        
        Args:
            atom_id: List of atom IDs where the feature should be added.
            name: Name/label of the new feature.
            representation: Visual representation style (e.g. ``'cartoon'``,
                ``'ball+stick'``). If ``None``, a default is chosen based on
                ``entity_type``.
            color: Color name or list for the feature. If ``None``, a default
                is chosen based on ``entity_type``.
            description: Optional dictionary of additional qualifiers.
            inplace: When ``True``, modifies ``self`` in place and returns it.
                When ``False`` (default), returns a copy with the new feature.
        
        Returns:
            Record object with the new feature annotation added.
        """
        # specify default representations
        DEFAULT_REPRESENTATIONS = {
            ENTITY_TYPE["solvent"]:"ball+stick",
            ENTITY_TYPE["ligand"]:"ball+stick",
            ENTITY_TYPE["ion"]:"ball+stick",
            ENTITY_TYPE["amino_acid"]:"cartoon",
            ENTITY_TYPE["nucleic_acid"]:"cartoon",
        }
        DEFAULT_COLORS = {
            ENTITY_TYPE["solvent"]:"ball+stick",
            ENTITY_TYPE["ligand"]:"ball+stick",
            ENTITY_TYPE["ion"]:"ball+stick",
            ENTITY_TYPE["amino_acid"]:"cartoon",
            ENTITY_TYPE["nucleic_acid"]:"cartoon",
        }
        REV_ENTITY_TYPE = {value:key for key, value in ENTITY_TYPE.items()}
        
        c = self.atom_array["atom_id"].isin(atom_id)
        df = self.atom_array[c]
        
        labels = [name, representation, color, description]
        out = []
        cols = ["chain_id","entity_type"]
        for chain_id, etype in df[cols].drop_duplicates().values:
            c = df["chain_id"] == chain_id
            c&= df["entity_type"] == etype
            atom_loc = Location(df.loc[c, "atom_id"].values)
            res_loc = Location(df.loc[c, "res_id"].values)
            
            # specify default representations
            if representation is None:
                labels[1] = DEFAULT_REPRESENTATIONS[etype]
            
            if color is None:
                labels[2] = DEFAULT_COLORS[etype]
            
            # build the dataframe
            out.append([atom_loc, chain_id, res_loc, REV_ENTITY_TYPE[etype]] + labels)
        cols = ["atom_id","chain_id","res_id","entity_type",
                "name","representation","color","description"]
        df = pd.DataFrame(out, columns=cols)

        if inplace:
            out = self
        else:
            out = self.copy()
        out.features = pd.concat([self.features, df])
        return out
        
    def superimpose_to(self, rec, atom_mask=None):
        """Superimposes another crystal structure onto this record.
        
        Aligns the atom arrays of ``self`` and ``rec`` using biotite's
        structural superposition, returning a new Record with the fitted
        coordinates.
        
        Args:
            rec: Record object whose structure will be superimposed onto
                this record.
            atom_mask: Optional boolean array selecting which atoms participate
                in the alignment. Atoms not in the mask are excluded from the
                superposition calculation.
        
        Returns:
            A new Record object with the coordinates of ``rec`` aligned to
            this record's coordinate frame.
        """
        out = self.copy()
        x1 = df_to_atom_array(rec.atom_array)
        x2 = df_to_atom_array(out.atom_array)
        result = biotite.structure.superimpose(
                    fixed=x1,
                    mobile=x2,
                    atom_mask=atom_mask
        )
        (fitted_arr, Aff) = result
        out.atom_array = atom_array_to_df(fitted_arr)
        return out
        
    def superimpose_homologs_to(self, rec, min_anchors:int=3):
        """Aligns this structure to another via homology-based superposition.
        
        Performs pairwise sequence alignment of residue backbones to identify
        matching anchors, then uses biotite's homologue superposition to fit
        ``rec`` onto this record's coordinate frame.
        
        Args:
            rec: Record object to align into this record's frame.
            min_anchors: Minimum number of anchor atom pairs required for a
                valid superposition. Lower values allow alignment with fewer
                matching residues.
        
        Returns:
            A new Record object with ``rec``'s coordinates aligned to this
            record's coordinate frame.
        """
        out = self.copy()
        x1 = df_to_atom_array(rec.atom_array)
        x2 = df_to_atom_array(out.atom_array)
        result = biotite.structure.superimpose_homologs(
                    fixed=x1, 
                    mobile=x2,
                    min_anchors=min_anchors
        )
        (fitted_arr, Aff, idx1, idx2) = result
        out.atom_array = atom_array_to_df(fitted_arr)
        return out

    def superimpose_by_sse_to(self, rec, min_anchors:int=3):
        """Aligns this structure to another using secondary structure element (SSE) matching.
        
        Performs pairwise alignment of secondary structure representations to
        identify matching anchor residues, then superimposes ``rec`` onto this
        record's coordinate frame while minimizing the TM score.
        
        Args:
            rec: Record object to align into this record's frame.
            min_anchors: Minimum number of anchor atom pairs required for a
                valid superposition.
        
        Returns:
            A new Record object with ``rec``'s coordinates aligned to this
            record's coordinate frame.
        """
        out = self.copy()
        x1 = df_to_atom_array(rec.atom_array)
        x2 = df_to_atom_array(out.atom_array)
        result = biotite.structure.superimpose_structural_homologs(
                    fixed=x1, 
                    mobile=x2,
        )
        (fitted_arr, Aff, idx1, idx2) = result
        out.atom_array = atom_array_to_df(fitted_arr)
        return out
        
    def __repr__(self):
        """Returns a human-readable representation of the Record.
        
        Displays summary information including filename, sequences,
        ligands, and features.
        
        Returns:
            Multi-line string showing the record's contents.
        """
        df_AA = self.get_amino_acids()
        df_nuc = self.get_nucleic_acids()
        df_seq = pd.concat([df_AA, df_nuc])
        df_lig = self.get_ligands()
        
        out = f"filename: {self.filename}\n\n"
        out+= f"List of Sequences:\n{df_seq}\n\n"
        out+= f"List of Ligands:\n{df_lig}\n\n"
        out+= f"List of Features:\n{self.features[FEATURE_COLUMNS]}\n\n"
        return out
        
    def __str__(self):
        return self.__repr__()
        
    def __delitem__(self, key):
        """Deletes features or atom slices from the Record.
        
        Supports boolean masks (to drop selected atoms) and slices (to remove
        a sequence range). Feature annotations are automatically trimmed to
        remain consistent with the remaining atoms.
        
        Args:
            key: Boolean array selecting atoms to keep (``~key`` are deleted),
                or a slice ``(start, stop)`` removing atoms in that range.
        
        Raises:
            ValueError: If ``key`` is not a boolean array or slice.
        """
        # deletes slices of atom_array by bool array
        if isinstance(key, np.ndarray) and key.dtype == bool:
            self.atom_array = self.atom_array.loc[~key]
            self.features = _trim_features(self.atom_array, self.features)
        
        # deletes slices of the atom array and associated features
        elif isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            self.atom_array = self.atom_array.loc[~c]
            self.features = _trim_features(self.atom_array, self.features)
        else:
            raise ValueError(f"Unsupported key type: {type(key)!r}")
        
    def __setitem__(self, key:tuple, value:"Record"):
        """Sets atom values and updates feature annotations accordingly.
        
        Supports replacing atom arrays with new Record objects via boolean mask
        or slice keys. Feature annotations are automatically trimmed to match
        the new atom configuration.
        
        Args:
            key: Tuple of (index/slice/boolean_mask, column_name) or just a
                boolean array/slice selecting which atoms to replace.
            value: New Record object whose atoms will replace the selected
                atoms in ``self``.
        
        Raises:
            ValueError: If ``key`` or ``value`` are of unsupported types.
        """
        if isinstance(value, Record) and isinstance(key, np.ndarray) and key.dtype == bool:
            self.atom_array = self.atom_array.loc[~key]
            features = _trim_features(self.atom_array, self.features)
            self = self + value
        
        elif isinstance(value, Record) and isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            self.atom_array = self.atom_array[~c]
            self.features = _trim_features(self.atom_array, self.features)
            self = self + value
        
        else:
            raise ValueError(f"Unsupported key type: {type(key)!r}")
        
    def __getitem__(self, key):
        """Retrieves a subset of the Record via indexing.
        
        Supports slice, boolean array, and integer array keys. When a boolean
        mask or index array is used, the returned Record's feature table is
        trimmed to match the selected atoms.
        
        Args:
            key: Selection key. Accepted forms:
                - ``slice`` ``(start, stop)`` — atoms in that range.
                - ``numpy.ndarray`` of bools — features at ``True`` positions.
                - ``pandas.Series`` of bools — same as numpy array.
                - ``str`` (column name) — returns that column from features.
                - ``(int, str)`` tuple — cell value at row ``index``, column ``name``.
        
        Returns:
            A new Record object containing the selected atoms and trimmed
            feature annotations, or a pandas Series/scalar depending on key type.
        
        Raises:
            ValueError: If ``key`` is of an unsupported type.
        """
        # slice through atom_array
        if isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            atom_array = self.atom_array[c].copy()
            features = _trim_features(atom_array, self.features)
        
        # slice via array indices or bool indices
        elif isinstance(key, np.ndarray) or isinstance(key, pd.Series):
            atom_array = self.atom_array.loc[key].copy()
            features = _trim_features(atom_array, self.features)

        else:
            raise ValueError(f"Unsupported key type: {type(key)!r}")
        
        rec = Record(filename = self.filename,
                     features = features,
                     atom_array = atom_array,
        )
        return rec
        
    def __add__(self, other):
        """
        Performs self + other
        
        Args:
            other: Record object
        
        Returns:
            modified copy of object
        """
        # sort and reassign atom_id of self
        out = self.reset_index()
        L = out.atom_array["atom_id"].max()
        
        # sort and reassign atom_id of other
        other = other.reset_index()
        other.atom_array["atom_id"]+= L
        other.features["atom_id"]+= L

        print('debug', other, out)
        
        # concatenate
        out.atom_array = pd.concat([out.atom_array, other.atom_array])
        out.features = pd.concat([out.features, other.features])
        return out

    def merge(self, other:"Record"):
        """
        Merges atoms from one Record into another

        Combines the atom arrays of both records and merges the feature
        annotations, combining annotation definitions of features which
        share a name.

        Args:
            other: Record object
        
        Returns:
            returns merged copy of self
        """
        atom_array = pd.concat([self.atom_array, other.atom_array])
        
        # populate joint atom_id slots with self
        cols = ["atom_id","chain_id","res_id"]
        atom_array = atom_array[cols].drop_duplicates()
        atom_array = atom_array.merge(self.atom_array, on=cols, how="left")
        self_features = _trim_features(atom_array, self.features)
        
        cols = ["atom_id","chain_id","res_id","res_name","atom_name","element"]
        x = atom_array[cols].drop_duplicates()
        if len(atom_array["atom_id"].drop_duplicates()) < len(x):
            raise ValueError(f"merge(): duplicates of atom_id entries found in self.atom_array")
        
        # populate joint atom_id slots with other
        c = atom_array["element"].notna()
        x1 = atom_array[c]
        
        cols = ["atom_id","chain_id","res_id"]
        atom_array = atom_array[~c, cols].drop_duplicates()
        x2 = atom_array.merge(other.atom_array, on=cols, how="left")
        other_features = _trim_features(x2, other.features)
        
        atom_array = pd.concat([x1, x2])
        
        cols = ["atom_id","chain_id","res_id","res_name","atom_name","element"]
        x = atom_array[cols].drop_duplicates()
        if len(atom_array["atom_id"].drop_duplicates()) < len(x):
            raise ValueError(f"merge(): duplicates of atom_id entries found in other.atom_array")
        
        # todo generate tests
        rec = Record(filename = self.filename,
                     features = pd.concat([self_features, other_features]),
                     atom_array = atom_array,
        )
        return rec
        
    def duplicated(self):
        """
        Find duplicated atom_array entries
        
        Returns:
            Boolean Series marking the duplicate features; the first
            occurrence of each set of annotations is ``False``.
        """
        # check unique atom entries are correct
        cols = ["atom_id","chain_id","res_id","res_name",
                "atom_name","hetero","element","entity_type"]
        x = self.atom_array[cols].drop_duplicates()        
        c = x[cols[:3]].duplicated()
        if sum(c) > 0:
            print(f"{sum(c)} duplicated atom_ids found")
            return self.atom_array["atom_id"].isin(x.loc[c, "atom_id"])
        
        # check for unique positions across all
        cols = ["atom_id", "chain_id", "res_id", "frame",
                "x", "y", "z", "b_factor", "occupancy", "charge"]
        x = self.atom_array[cols].drop_duplicates()
        c = x[cols[:4]].duplicated()
        if sum(c) > 0:
            print(f"{sum(c)} duplicated positions found")
            return self.atom_array["atom_id"].isin(x.loc[c, "atom_id"])
        return self.atom_array.duplicated()
