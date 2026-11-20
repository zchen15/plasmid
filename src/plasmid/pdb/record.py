# Module for browsing and editing pdb and mmCIF files

# science libs
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
    """
    
    atom_array = biotite.structure.io.load_structure(filename)
    if isinstance(atom_array, biotite.structure.AtomArrayStack):
        # todo vectorize better
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

        out = f"filename: {self.filename}\n\n"
        out+= f"List of Sequences:\n{df_seq}\n\n"
        out+= f"List of Ligands:\n{df_lig}\n\n"
        out+= f"List of Features:\n{self.features}\n\n"
        return out
        
    def __str__(self):
        return self.__repr__()
        
    def __delitem__(self, key):
        """
        Deletes features or lists of features
        
        Args:
            key: list of bools
        
        Returns:
            modified copy of self
        """
        self.atom_array = self.atom_array.loc[~key]
        # todo update feature array
        
    def __setitem__(self, key:tuple, value:str):
        """
        Set items in atom_array
        """
        if isinstance(key, tuple) and isinstance(key[1], str):
            self.atom_array.loc[key] = value
        # todo update feature array
        
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
        
        # return position value
        elif isinstance(key, tuple) and isinstance(key[1], str):
            return self.atom_array.loc[key]
        
        # slice through atom_ids
        elif isinstance(key, slice):
            c = self.atom_array["atom_id"] >= key.start
            c&= self.atom_array["atom_id"] < key.stop
            keep = [loc.contains(list(range(start, stop))) for loc in self.features["atom_id"]]
            rec = Record(filename = self.filename,
                         features = self.features[keep].copy(),
                         atom_array = self.atom_array.loc[c].copy(),
            )
            return rec
        # slice via array indices or bool indices
        elif isinstance(key, np.ndarray) or isinstance(key, pd.Series):
            # todo clean up features
            rec = Record(filename = self.filename,
                         features = self.features.copy(),
                         atom_array = self.atom_array.loc[key].copy(),
            )
            return rec
        else:
            raise ValueError(f"__getitem__(self, key) has invalid key = {key}")
        
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
        
    def merge(self, other):
        """
        Merges atoms from one Record into another
        
        todo fixes
        
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
        
        # todo merge annotation definitions
        out.atom_array = pd.concat([x1, x2])
        out.features = pd.concat([out.features, other.features])
        return out.reset_index()

    def duplicated(self):
        """
        Todo, find duplicated atom positions
        """
    
    def drop_duplicates(self):
        """
        Drop duplicates annotations
        
        Todo
        """
        out = self.copy()
        out.atom_array = out.atom_array.drop_duplicates()
        return out
