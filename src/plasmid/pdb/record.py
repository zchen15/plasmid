# Module for browsing and editing pdb and mmCIF files

# science libs
import numpy as np
import pandas as pd
import biotite.structure.io

# custom libraries
from .locations import Location
from .io import atom_array_to_df, df_to_atom_array, read_pdb, ENTITY_TYPE
from .sequences import get_sequence, get_ligands

import copy

def read_pdb(filename:str):
    """
    Reads a pdb or mmCIF file
    Args:
        filename: filename to read
        
    Returns:
        Record object
    """
    
    atom_array = biotite.structure.io.load_structure(filename)
    if isinstance(atom_array, biotite.structure.AtomArrayStack):
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
                 atom_array = atom_array,
                )
    return rec

class Record:
    """
    This class holds functions pertaining to modifying the atoms in a pdb file
    """
    filename = None
    atom_array = None
    features = None
    
    def __init__(self, filename: str, **kwargs):
        """
        Initializes the Record object

        Args:
            filename: name of pdb or mmCIF file to read.
            kwargs: settings that override the class defaults.

        Returns:
            Record object
        """
        self.filename = filename
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        
    def copy(self):
        """
        Returns a copy of self
        """
        return copy.deepcopy(self)
        
    def reset_index(self):
        """
        Function for sorting record info
        """
        out = self.copy()
        
        out.atom_array = out.atom_array.sort_values(by=["chain_id","res_id","res_name","atom_name","element"]).reset_index(drop=True)
        out.atom_array["atom_id"] = np.arange(len(atom_array))
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
        
        out = []
        cols = ["chain_id","res_name","res_id"]
        for chain_id, res_name in df[cols].drop_duplicates().values:
            c = df["chain_id"] == chain_id
            c&= df["res_name"] == res_name
            
            # get positions
            atom_id = Location(df.loc[c, "atom_id"].values)
            res_id = Location(df.loc[c, "res_id"].values)
            
            out.append([atom_id, chain_id, res_id, res_name])
        out = pd.DataFrame(out, columns=["atom_id", "chain_id", "res_id", "res_name"])
        return out
        
    def get_amino_acids(self, frame:int=0):
        c1= self.atom_array["entity_type"] == ENTITY_TYPE["amino_acid"]
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
            seqlist = ["-"]*residues["res_id"].max()
            for ridx, res_name in residues[["res_id","res_name"]].values:
                try:
                    seqlist[i] = biotite.sequence.ProteinSequence.convert_letter_3to1(res_name)
                except Exception:
                    seqlist[i] = "-"
            
            out.append([atom_id, chain_id, res_id, ''.join(seqlist)])
        out = pd.DataFrame(out, columns=["atom_id", "chain_id", "res_id", "sequence"])
        return out
        
    def get_nucleic_acids(self, frame:int=0):
        c1= self.atom_array["entity_type"] == ENTITY_TYPE["nucleic_acid"]
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
            seqlist = ["-"]*residues["res_id"].max()
            for ridx, res_name in residues[["res_id","res_name"]].values:
                seqlist[i] = res_name
            
            out.append([atom_id, chain_id, res_id, ''.join(seqlist)])
        out = pd.DataFrame(out, columns=["atom_id", "chain_id", "res_id", "sequence"])
        return out
        
    def annotate_by_atom_id(self, name: str, atom_id:list):
        """
        Add annotation to pdb structure
        
        Args:
            name: name of feature
        
        Returns:
            copy of pdb record object
        """
        c = self.atom_array["atom_id"].isin(atom_id)
        df = self.atom_array[c]
        
        out = []
        cols = ["chain_id","res_id","entity_type"]
        for chain_id, etype in df[cols].drop_duplicates().values:
            c = df["chain_id"] == chain_id
            c&= df["entity_type"] == etype
            atom_loc = Location(df.loc[c, "atom_id"].values)
            res_loc = Location(df.loc[c, "res_id"].values)
            ename = ENTITY_TYPE[etype]
            out.append([atom_loc, chain_id, res_loc, f"{etype}:{ename}", name])
        cols = ["atom_id","chain_id","res_id","entity_type","name"]
        out = pd.DataFrame(out, columns=cols)
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
        
    def __len__(self):
        """
        Number of atoms
        
        Returns:
            length of atom_array
        """
        return self.atom_array["atom_id"].max()
        
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
        
        out = f"List of Sequences:\n{df_seq}\n\n"
        out+= f"List of Ligands:\n{df_lig}\n\n"
        return out
        
    def __delitem__(self, key):
        """
        Deletes a list of atoms in the AtomArray
        
        Args:
            key: list of bools
        
        Returns:
            modified copy of self
        """
        self.atom_array = self.atom_array.loc[~key]
        # todo clean up features

    def __setitem__(self, key:tuple, value:str):
        if isinstance(key, tuple) and isinstance(key[1], str):
            self.atom_array.loc[key] = value
        # todo clean up ceatures
    
    def __getitem__(self, key):
        """
        Get slices of an AtomArray
        
        Args:
            key: slice or list of bools or strings
        
        Returns:
            copy of the values of the underlying data
        """
        # return annotation values
        if isinstance(key, str):
            return self.atom_array[key]
        
        # slice through atom_ids
        elif isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            # todo clean up features
            rec = Record(filename = self.filename,
                         feature = self.features.copy(),
                         atom_array = self.atom_array.loc[c].copy(),
            )
        
        # slice via array indices or bool indices
        elif isinstance(key, np.ndarray):
            # todo clean up features
            rec = Record(filename = self.filename,
                         feature = self.features.copy(),
                         atom_array = self.atom_array.loc[key].copy(),
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
        # todo reset atom_id of other
        L = self.atom_array["atom_id"].max()
        other_atom_array = other.atom_array.copy()
        other_atom_array["atom_id"]+= L
        other_features = other.features.copy()
        other_features["atom_id"]+=L
        
        atom_array = pd.concat([self.atom_array, other_atom_array])
        features = pd.concat([self.features, other_features])
        rec = Record(filename = self.filename,
                     atom_array = atom_array,
                     features = features,
                    )
        return rec

    def merge(self, other):
        """
        Merges atoms from one Record into another
        
        Args:
            other: Record object
        
        Returns:
            returns merged copy of self
        """
        out = self.copy()
        df = pd.concat([out.atom_array, other.atom_array])
        
        # populate joint atom_id slots with self
        cols = ["atom_id","chain_id","res_id","res_name","frame"]
        df = df[cols].drop_duplicates().merge(out.atom_array, on=cols, how="left")
        
        # populate joint atom_id slots with other
        c = df["element"].notna()
        x1 = df[c]
        x2 = df.loc[~c, cols].drop_duplicates().merge(other.atom_array, on=cols, how="left")
        
        out.atom_array = pd.concat([x1, x2])
        return out
        
    def drop_duplicates(self):
        out = self.copy()
        out.atom_array = out.atom_array.drop_duplicates()
        return out
