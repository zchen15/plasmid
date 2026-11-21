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

def read_pdb(filename:str):
    """
    Reads a pdb or mmCIF file
    Args:
        filename: filename to read
    
    Returns:
        Record object

    todo fixes, if each frame does not have the same number of total atoms,
    assume chain_id, res_id, res_name, element, atom_name uniquely identify the atom, then assign atom_id based on that
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
        """
        Initializes the Record object
        
        Args:
            kwargs: settings that override the class defaults.
        
        Returns:
            Record object
        """
        self.atom_array = pd.DataFrame(columns=ATOM_ARRAY_COLUMNS)
        self.features = pd.DataFrame(columns=FEATURE_COLUMNS)
        
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        
    def copy(self):
        """
        Returns a copy of self
        """
        rec = Record(filename = self.filename,
                     atom_array = self.atom_array.copy(),
                     features = self.features.copy(),
                    )
        return rec
        
    def reset_index(self):
        """
        Function for sorting record info
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
        out = Record(filename = self.filename,
                     atom_array = atom_array,
                     features = features,
                    )
        return out
        
    def save(self, filename:str, **kwargs):
        """
        Saves AtomArray to files
        
        Args:
            filename: name of file to save
        """
        atom_array = df_to_atom_stack(self.atom_array)
        biotite.structure.io.save_structure(filename, atom_array, **kwargs)
        
    def get_ligands(self, frame:int=0):
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
        """
        Add annotation to pdb structure
        
        Args:
            name: name of feature
        
        Returns:
            copy of pdb record object
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
        """
        Superimpose one crystal structure to another
        
        Args:
            rec: Record object
            atom_mask: bool array which selects which atoms will be used for alignment
        
        Returns:
            returns aligned copy of self
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
        """
        Align one crystal structure into another by first doing
        pairwise alignment of residue sequences to identify matching
        backbone residues. Then, perform superimpose with outliers on
        matching residues.
        
        Args:
            rec: Record object
            min_anchors: min number of anchor atoms
        
        Returns:
            returns aligned copy of self
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
        """
        Align one crystal structure into another by first doing
        pairwise alignment of sequence representation of secondary
        structures to identify matching anchor backbone residues. 
        Next, the structures are superimposed using these anchor
        residues while minimizing the Tm score.
        
        Args:
            rec: Record object
            min_anchors: min number of anchor atoms
        
        Returns:
            returns aligned copy of self
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
        """
        Gets representation of underlying data
        
        Returns:
            string representation of the underlying data
        """
        df_AA = self.get_amino_acids()
        df_nuc = self.get_nucleic_acids()
        df_seq = pd.concat([df_AA, df_nuc])
        df_lig = self.get_ligands()
        
        out = f"filename: {self.filename}\n\n"
        out+= f"List of Sequences:\n{df_seq}\n\n"
        out+= f"List of Ligands:\n{df_lig}\n\n"
        out+= f"List of Features:\n{self.features}\n\n"
        return out
        
    def __str__(self):
        return self.__repr__()
        
    def _trim_features(self, atom_array):
        """
        Trim the feature table so that it only references atoms which are
        present in ``atom_array``.
        
        Features whose atoms were all removed are dropped; features which
        only partially overlap are trimmed to the surviving atoms.
        
        Args:
            atom_array: Atom dataframe to trim the features against.
        
        Returns:
            Copy of the feature table trimmed to ``atom_array``.
        """
        features = self.features.copy()
        if len(features) == 0:
            return features
        
        kept = sorted({int(p) for p in atom_array["atom_id"].tolist()}) if len(atom_array) else []
        keep_mask = []
        for i, loc in enumerate(features["atom_id"].values):
            positions = np.asarray(loc.positions, dtype=int) if isinstance(loc, Location) else np.asarray([loc], dtype=int)
            surviving = [int(p) for p in positions if int(p) in kept]
            if len(surviving) == 0:
                keep_mask.append(False)
            else:
                keep_mask.append(True)
                features.iloc[i, features.columns.get_loc("atom_id")] = Location(surviving)
                features.iloc[i, features.columns.get_loc("res_id")] = Location(surviving)
        
        return features[np.asarray(keep_mask, dtype=bool)].reset_index(drop=True)

    def __delitem__(self, key):
        """
        Deletes parts of the atom_array
        """
        # deletes slices of atom_array by bool array
        if isinstance(key, np.ndarray) and key.dtype == bool:
            self.atom_array = self.atom_array.loc[~key]
            self.features = self._trim_features(self.atom_array)
        
        # deletes slices of the atom array and associated features
        elif isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            self.atom_array = self.atom_array.loc[~c]
            self.features = self._trim_features(self.atom_array)
        else:
            raise ValueError(f"__setitem__(): Unsupported: {type(key)} {type(value)}")
        
    def __setitem__(self, key:tuple, value:Record):
        """
        Replaces a set of atoms with new atoms
        """
        if isinstance(value, Record) and isinstance(key, np.ndarray) and key.dtype == bool:
            self.atom_array = self.atom_array.loc[~key]
            features = self._trim_features(self.atom_array)
            self = self + value
        
        elif isinstance(value, Record) and isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            self.atom_array = self.atom_array[~c]
            self.features = self._trim_features(self.atom_array)
            self = self + value
        
        else:
            raise ValueError(f"__setitem__(): Unsupported: {type(key)} {type(value)}")
        
    def __getitem__(self, key):
        """
        Get slices of an AtomArray
        
        Args:
            key: slice or list of bools or strings
        
        Returns:
            copy of the values of the underlying data
        """
        # slice through atom_array
        if isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            atom_array = self.atom_array[c].copy()
            features = self._trim_features(atom_array)
        
        # slice via array indices or bool indices
        elif isinstance(key, np.ndarray) or isinstance(key, pd.Series):
            atom_array = self.atom_array.loc[key].copy()
            features = self._trim_features(atom_array)

        else:
            raise ValueError(f"__setitem__(): Unsupported: {type(key)} {type(value)}")
        
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
        
        # concatenate
        out.atom_array = pd.concat([out.atom_array, other.atom_array])
        out.features = pd.concat([out.features, other.features])
        return out

    def merge(self, other:Record):
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
        self_features = self._trim_features(atom_array)
        
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
        other_features = other._trim_features(x2)
        
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
