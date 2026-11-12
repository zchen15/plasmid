# Module for browsing and editing pdb and mmCIF files

# science libs
import numpy as np
import pandas as pd

# custom libraries
from .locations import Location
from ..utils import load_libraries

import biotite.structure.io
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

import copy
import os
import tempfile


def _first(value):
    """
    Return the first value of a sequence-like value, otherwise the value itself.

    Args:
        value: scalar or sequence-like value.

    Returns:
        A scalar value or None.
    """
    if isinstance(value, (list, tuple, np.ndarray)):
        if len(value) == 0:
            return None
        return value[0]
    return value


def _text(value):
    """
    Convert a structure file value to a stripped string.

    Args:
        value: scalar or sequence-like value.

    Returns:
        A string or None.
    """
    # unwrap biotite CIF column objects
    if hasattr(value, 'as_array'):
        try:
            value = value.as_array()
        except Exception:
            try:
                value = value.as_item()
            except Exception:
                pass
    value = _first(value)
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode('utf-8')
    return str(value).strip()


def _column(category, key: str) -> list:
    """
    Extract a column from a biotite CIF category.

    Args:
        category: biotite CIFCategory object.
        key: column name.

    Returns:
        A list of column values, or an empty list if the column is absent.
    """
    try:
        column = category[key]
    except Exception:
        return []
    try:
        return column.as_array().tolist()
    except Exception:
        try:
            return [column.as_item()]
        except Exception:
            return []

def read_pdb(filename:str):
    """
    Reads a pdb or mmCIF file
    Args:
        filename: filename to read
        
    Returns:
        Record object
    """
    return Record(filename)
    
class Record:
    """
    This class holds functions pertaining to modifying the atoms in a pdb file
    
    Attributes:
        atom_array: biotite.structure.AtomArray object
        filename: name of the pdb or mmCIF file that was read
        annotations: pandas dataframe of annotations associated with the AtomArray
        sequence_annotations: pandas dataframe of sequence information from get_sequence()
        ligand_annotations: pandas dataframe of ligand information from get_ligands()
        chain_colors: dictionary mapping each chain_id to a default color
    """

    # default colors assigned to chains in the order they are encountered
    default_chain_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                            '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                            '#bcbd22', '#17becf']

    # annotations are stored in a dataframe with columns
    # [name, feature, location, color, description] similar to genbank features
    annotations: pd.DataFrame = pd.DataFrame(
        columns=["name", "feature", "location", "color", "description"]
    )
    sequence_annotations: pd.DataFrame = pd.DataFrame(
        columns=["atom_id", "chain_id", "res_id", "sequence", "n_frames", "color"]
    )
    ligand_annotations: pd.DataFrame = pd.DataFrame(
        columns=["atom_id", "chain_id", "res_name", "res_id", "count",
                 "n_frames", "color"]
    )
    chain_colors: dict = {}
    atom_array = None
    filename = None

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
        self.atom_array = biotite.structure.io.load_structure(filename)
        self.chain_colors = {}
        self.annotations = pd.DataFrame(
            columns=["name", "feature", "location", "color", "description"]
        )
        self.sequence_annotations = pd.DataFrame(
            columns=["atom_id", "chain_id", "res_id", "sequence", "n_frames", "color"]
        )
        self.ligand_annotations = pd.DataFrame(
            columns=["atom_id", "chain_id", "res_name", "res_id", "count",
                     "n_frames", "color"]
        )

        # parse annotations that are stored in the structure file
        self.annotations = self.parse_structure_annotations(filename)

        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)

        # populate the cached sequence and ligand annotation tables
        self.get_sequence()
        self.get_ligands()
        self.get_chain_colors()
        
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
        df = out.to_dataframe()
        c = df["frame"] == 0
        x = df[c].sort_values(by=["chain_id","res_id","res_name","atom_name","element"]).reset_index(drop=True)
        
        if len(df["frame"].drop_duplicates()) > 1:
            atom_stack = []
            for frame in df["frame"].drop_duplicates().values:
                atom_stack+= [out.atom_array[frame][x["atom_id"].values]]
            out.atom_array = biotite.structure.stack(atom_stack)
        else:
            out.atom_array = out.atom_array[x["atom_id"].values]
        return out
        
    def save(self, filename:str, **kwargs):
        """
        Saves AtomArray to files
        
        Args:
            filename: name of file to save
        """
        biotite.structure.io.save_structure(filename, self.atom_array, **kwargs)

    def to_mmcif(self, filename:str):
        """
        Saves AtomArray to mmcif file

        Args:
            filename: name of file to save
        """
        # get polypeptide info
        df = self.to_dataframe()
        c = df["frame"] == 0
        df = df[c]
        
        # get list of peptides
        c0 = biotite.structure.filter_amino_acids(self.atom_array)
        cols = ["chain_id","res_id","res_name"]
        df_pep = df.loc[c0, cols].drop_duplicates()
        df_pep = df_pep.sort_values(by=cols)
        
        # get list of nucleotides
        c0 = biotite.structure.filter_nucleotides(self.atom_array)
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
        df = self.get_sequence()
        x = x[["chain_id","entity_id","entity_type","poly_type"]].drop_duplicates()
        x = x.merge(df[["chain_id","sequence"]], on=["chain_id"], how="left")
        
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
        
        block = {"entity":entity,
                 "entity_poly":entity_poly,
                 "struct_asym":struct_asym,
                 "entity_poly_seq":entity_poly_seq}
        block = biotite.structure.io.pdbx.CIFBlock(block)
        file = biotite.structure.io.pdbx.CIFFile({"header":block})
        
        # add force add extra fields to cif file so Boltz doesn"t complain
        required_extra_fields = ["b_factor","occupancy","charge"]
        for key in required_extra_fields:
            if key not in self.atom_array.get_annotation_categories():
                if key=="occupancy":
                    L = len(self.atom_array.get_annotation("chain_id"))
                    self.atom_array.set_annotation(key, [1]*L)
                else:
                    L = len(self.atom_array.get_annotation("chain_id"))
                    self.atom_array.set_annotation(key, [0]*L)
        
        # add structure info
        biotite.structure.io.pdbx.set_structure(file, self.atom_array)
        file.write(filename)
        
    def atom_array_to_dataframe(self, atom_array):
        """
        Convert atom_array to pandas dataframe
        """
        df = pd.DataFrame(np.arange(len(atom_array)), columns=["atom_id"])
        df["chain_id"] = atom_array.chain_id
        df["res_id"] = atom_array.res_id
        df["res_name"] = atom_array.res_name
        df["hetero"] = atom_array.hetero
        df["atom_name"] = atom_array.atom_name
        df["element"] = atom_array.element
        df["x"] = atom_array.coord[:,0]
        df["y"] = atom_array.coord[:,1]
        df["z"] = atom_array.coord[:,2]
        return df
        
    def to_dataframe(self):
        """
        Convert atom array information to a dataframe
        """
        # process as stack of items
        if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
            out = []
            for frame in range(len(self.atom_array)):
                df = self.atom_array_to_dataframe(self.atom_array[frame])
                df["frame"] = frame
                out.append(df)
            df = pd.concat(out)
        else:
            df = self.atom_array_to_dataframe(self.atom_array)
            df["frame"] = 0
        #df["frame"] = df["frame"].astype(int)
        return df
        
    def get_chain_color(self, chain_id) -> str:
        """
        Get the default color assigned to a chain, assigning one if needed.

        Args:
            chain_id: chain identifier.

        Returns:
            A color string.
        """
        if chain_id not in self.chain_colors:
            i = len(self.chain_colors) % len(self.default_chain_colors)
            self.chain_colors[chain_id] = self.default_chain_colors[i]
        return self.chain_colors[chain_id]

    def get_chain_colors(self) -> dict:
        """
        Get the default color assigned to each chain.

        Chains that have not been seen by get_sequence() or get_ligands()
        (for example nucleotide-only chains) are assigned a color here.

        Returns:
            Dictionary mapping chain_id to color.
        """
        if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
            chain_ids = self.atom_array[0].chain_id
        else:
            chain_ids = self.atom_array.chain_id
        for chain_id in np.unique(chain_ids):
            self.get_chain_color(chain_id)
        return dict(self.chain_colors)

    def get_sequence(self):
        """
        Get list of protein and amino acid sequences

        The result is cached as the sequence_annotations attribute.

        Returns:
            pandas dataframe of chain_id and sequence
        """
        data = self.to_dataframe()
        c = data["frame"] == 0
        data = data[c]
        c0 = biotite.structure.filter_amino_acids(self.atom_array)
        c1 = data["res_name"] == "HOH"
        
        # extract amino acid sequences
        out = []
        for chain_id in data.loc[c0, "chain_id"].drop_duplicates().values:
            c2 = data["chain_id"] == chain_id
            
            # get the indices
            c = c0 & ~c1 & c2
            atom_id = Location(data.loc[c, "atom_id"].values)
            res_id = Location(data.loc[c, "res_id"].values)
            
            # get the sequence
            x = data.loc[~c1 & c2, ["res_name","res_id","hetero"]].drop_duplicates()
            seqlist = x["res_name"].values
            for i in range(len(seqlist)):
                try:
                    seqlist[i] = biotite.sequence.ProteinSequence.convert_letter_3to1(seqlist[i])
                except:
                    seqlist[i] = "-"
            seq = "".join(seqlist)
            
            out.append([atom_id, chain_id, res_id, seq])
        cols = ["atom_id", "chain_id","res_id","sequence"]
        df = pd.DataFrame(out, columns=cols)
        
        # add info about number of poses or frames
        if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
            df["n_frames"] = len(self.atom_array)
        else:
            df["n_frames"] = 0

        # assign a default color to each chain
        df["color"] = [self.get_chain_color(c) for c in df["chain_id"].values]

        self.sequence_annotations = df
        return df
        
    def get_ligands(self):
        """
        Get list of ligands in crystal structure

        The result is cached as the ligand_annotations attribute.

        Returns:
            pandas dataframe of ligands
        """
        data = self.to_dataframe()
        c = data["frame"] == 0
        data = data[c]
        c1 = biotite.structure.filter_amino_acids(self.atom_array)
        c1|= biotite.structure.filter_nucleotides(self.atom_array)
        data = data[~c1]
        
        out = []
        for chain_id, res_name in data[["chain_id","res_name"]].drop_duplicates().values:
            c0 = data["chain_id"] == chain_id
            c0&= data["res_name"] == res_name
            atom_id = Location(data.loc[c0, "atom_id"].values)
            res_id = Location(data.loc[c0, "res_id"].values)
            out.append([atom_id, chain_id, res_name, res_id, len(res_id)])
        
        cols = ["atom_id", "chain_id", "res_name","res_id","count"]
        df = pd.DataFrame(out, columns=cols)
        
        # add info about number of poses or frames
        if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
            df["n_frames"] = len(self.atom_array)
        else:
            df["n_frames"] = 0

        # assign a default color to each chain
        df["color"] = [self.get_chain_color(c) for c in df["chain_id"].values]

        self.ligand_annotations = df
        return df

    def parse_structure_annotations(self, filename: str) -> pd.DataFrame:
        """
        Parse annotations from the header of a pdb or mmCIF structure file.

        The annotations dataframe stores at least the [name, description]
        columns. Structure-level metadata (title, keywords, entity
        descriptions) and secondary structure records are read when they are
        available in the file.

        Args:
            filename: name of pdb or mmCIF file to read.

        Returns:
            pandas dataframe of annotations
        """
        cols = ["name", "feature", "location", "color", "description"]
        rows = []

        # strip compression extensions before detecting the file type
        name = os.path.basename(str(filename)).lower()
        if name.endswith('.gz'):
            name = name[:-3]

        try:
            if name.endswith(('.cif', '.mmcif', '.bcif')):
                if name.endswith('.bcif'):
                    cif = pdbx.BinaryCIFFile.read(filename)
                else:
                    cif = pdbx.CIFFile.read(filename)
                block = cif.block

                # structure title and keywords
                title = None
                if 'struct' in block:
                    values = _column(block['struct'], 'title')
                    title = _text(values[0]) if values else None
                keywords = None
                if 'struct_keywords' in block:
                    values = _column(block['struct_keywords'], 'pdbx_keywords')
                    keywords = _text(values[0]) if values else None
                if title:
                    rows.append(['title', 'structure', None, None, title])
                if keywords:
                    rows.append(['keywords', 'structure', None, None, keywords])

                # entity descriptions, one per chain
                if 'entity' in block:
                    entity = block['entity']
                    ids = _column(entity, 'id')
                    descriptions = _column(entity, 'pdbx_description')
                    for i, entity_id in enumerate(ids):
                        desc = _text(descriptions[i]) if i < len(descriptions) else None
                        if desc:
                            rows.append([f'entity {_text(entity_id)}', 'chain',
                                         None, None, desc])

                # helix and turn annotations
                if 'struct_conf' in block:
                    conf = block['struct_conf']
                    ids = _column(conf, 'id')
                    types = _column(conf, 'conf_type_id')
                    asym = _column(conf, 'beg_label_asym_id')
                    beg = _column(conf, 'beg_label_seq_id')
                    end = _column(conf, 'end_label_seq_id')
                    for i, conf_id in enumerate(ids):
                        label = _text(types[i]) if i < len(types) else 'secondary structure'
                        chain = _text(asym[i]) if i < len(asym) else None
                        start = _text(beg[i]) if i < len(beg) else None
                        stop = _text(end[i]) if i < len(end) else None
                        desc = f'{label} {chain}:{start}-{stop}'
                        rows.append([_text(conf_id) or f'conf {i+1}',
                                     'secondary_structure', None, None, desc])

                # beta sheet annotations
                if 'struct_sheet_range' in block:
                    sheet = block['struct_sheet_range']
                    ids = _column(sheet, 'id')
                    asym = _column(sheet, 'beg_label_asym_id')
                    beg = _column(sheet, 'beg_label_seq_id')
                    end = _column(sheet, 'end_label_seq_id')
                    for i, sheet_id in enumerate(ids):
                        chain = _text(asym[i]) if i < len(asym) else None
                        start = _text(beg[i]) if i < len(beg) else None
                        stop = _text(end[i]) if i < len(end) else None
                        desc = f'sheet {chain}:{start}-{stop}'
                        rows.append([_text(sheet_id) or f'sheet {i+1}',
                                     'secondary_structure', None, None, desc])
            else:
                pdb_file = pdb.PDBFile.read(filename)
                header = pdb_file.get_header()

                title = _text(header.get('title'))
                keywords = _text(header.get('keywords'))
                if title:
                    rows.append(['title', 'structure', None, None, title])
                if keywords:
                    rows.append(['keywords', 'structure', None, None, keywords])
                for key in ['compound', 'source']:
                    value = header.get(key)
                    if value:
                        rows.append([key, 'structure', None, None, str(value)])
        except Exception:
            # annotations are optional, so ignore files without usable metadata
            pass

        return pd.DataFrame(rows, columns=cols)

    def get_protein_sse(self):
        """
        Compute secondary structure of protein sequence using P-SEA algorithm

        Args:
            frame: which slice of the atom_stack to extract from
        
        Returns:
            pandas dataframe of secondary structures
        """
        df = self.get_sequence()
        out = []
        for atom_id in df["atom_id"].values:
            s1 = min(atom_id.positions)
            s2 = max(atom_id.positions)
            x = self.atom_array[s1:s2]
            sse = biotite.structure.annotate_sse(x)
            sse = "".join(sse)
            out.append(sse)
        df["dssp"] = out
        return df
        
    def get_dssp_annotations(self):
        """
        Get secondary structure annotations for protein sequence
        
        Args:
            frame: which slice of the atom_stack to extract from
        
        Returns:
            pandas dataframe of annotations
        """
        data = self.to_dataframe()
        c = data["frame"] == 0
        data = data[c]
        df = self.get_protein_sse()
        
        out = []
        cols = ["chain_id","atom_id","dssp"]
        for chain_id, atom_id, sse in df[cols].values:
            c = data["atom_id"] >= min(atom_id.positions)
            c&= data["atom_id"] < max(atom_id.positions)
            c&= data["chain_id"] == chain_id
            idx1 = data.loc[c, "res_id"].drop_duplicates().values
            sse = np.array([i for i in sse])
            idx2 = np.arange(len(sse))
            
            c0 = sse[:-1] == sse[1:]
            c1 = np.concatenate(([False], c0))
            c2 = np.concatenate((c0, [False]))
            s1 = idx1[~c1]
            s2 = idx1[~c2]
            sse = sse[~c1]
            
            for i in range(len(s1)):
                if sse[i] != "c":
                    # atom index for the secondary structure
                    c1 = data["chain_id"] == chain_id
                    c1&= data["res_id"] >= s1[i]
                    c1&= data["res_id"] <= s2[i]
                    atom_id = Location(data.loc[c1, "atom_id"].values)
                    res_id = Location(data.loc[c1, "res_id"].values)
                    L = s2[i] - s1[i]
                    out.append([atom_id, chain_id, res_id, sse[i], L])
        cols = ["atom_id", "chain_id", "res_id", "dssp", "length"]
        df = pd.DataFrame(out, columns=cols)
        return df
        
    def get_annotations(self):
        """
        Get list of annotations for pdb structure
        
        Returns:
            pandas dataframe of annotations
        """
        return self.annotations
        
    def annotate_by_residue_id(self, name: str, start: int, end: int, color: str):
        """
        Add annotation to pdb structure
        
        Args:
            name: name of feature
            start: starting residue index
            end: stopping residue index
            color: color to use
        
        Returns:
            copy of pdb record object
        """
        out = self.copy()
        loc = Location([start, end])
        row = pd.DataFrame(
            [[name, 'residue', loc, color, f'{start}-{end}']],
            columns=['name', 'feature', 'location', 'color', 'description'],
        )
        out.annotations = pd.concat([out.annotations, row], ignore_index=True)
        return out
        
    def annotate_by_atom_id(self, name: str, start: int, end: int, color: str):
        """
        Add annotation to pdb structure
        
        Args:
            name: name of feature
            start: starting atom index
            end: stopping atom index
            color: color to use
        
        Returns:
            copy of pdb record object
        """
        out = self.copy()
        loc = Location([start, end])
        row = pd.DataFrame(
            [[name, 'atom', loc, color, f'{start}-{end}']],
            columns=['name', 'feature', 'location', 'color', 'description'],
        )
        out.annotations = pd.concat([out.annotations, row], ignore_index=True)
        return out
        
        
    def add_bonds(self):
        """
        Add bondlist using distance and charge information to the AtomArray
        """
        df = self.to_dataframe()
        c = df["frame"] == 0
        df = df[c]
        c0 = biotite.structure.filter_amino_acids(self.atom_array)
        c0|= biotite.structure.filter_nucleotides(self.atom_array)
        
        atom_list = []
        for chain_id in df["chain_id"].drop_duplicates():
            c1 = df["chain_id"] == chain_id
            c_poly = c0 & c1
            c_lig = ~c0 & c1
            if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
                poly_atom_array = self.atom_array[:, c_poly]
                lig_atom_array = self.atom_array[:, c_lig]
            else:
                poly_atom_array = self.atom_array[c_poly]
                lig_atom_array = self.atom_array[c_lig]
            
            # bonds for residues
            bondlist = biotite.structure.connect_via_residue_names(poly_atom_array)
            poly_atom_array.bonds = bondlist
            atom_list.append(poly_atom_array)
            
            # bonds for ligands
            if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
                bondlist = biotite.structure.connect_via_distances(lig_atom_array[0])
            else:
                bondlist = biotite.structure.connect_via_distances(lig_atom_array)
            lig_atom_array.bonds = bondlist
            atom_list.append(lig_atom_array)
        
        # concatenate atom_array
        for i in range(1, len(atom_list)):
            atom_list[0] = atom_list[0] + atom_list[i]
        
        # return copy
        out = self.copy()
        out.atom_array = atom_list[0]
        return out

    def merge(self, rec):
        """
        Merges atoms from one Record into another
        
        Args:
            rec: Record object
        
        Returns:
            returns merged copy of self
        """
        out = self.copy()
        
        # make both atom_arrays the same type
        c1 = isinstance(out.atom_array, biotite.structure.AtomArray)
        c2 = isinstance(out.atom_array, biotite.structure.AtomArrayStack)
        c3 = isinstance(rec.atom_array, biotite.structure.AtomArray)
        c4 = isinstance(rec.atom_array, biotite.structure.AtomArrayStack)
        
        if c2 & c4:
            assert len(out.atom_array) == len(rec.atom_array)
        elif c1 & c4:
            L = len(rec.atom_array)
            out.atom_array = biotite.structure.stack([out.atom_array]*L)
        elif c2 & c3:
            L = len(out.atom_array)
            rec.atom_array = biotite.structure.stack([rec.atom_array]*L)
        
        df1 = self.to_dataframe()
        df2 = rec.to_dataframe()
        df2 = df2.rename(columns={"atom_id":"atom_id_src"})
        
        # figure out which atoms to replace or add
        cols = ["chain_id","res_id","res_name","atom_name","element"]
        df2 = df2[cols+["atom_id_src"]].drop_duplicates()
        df2 = df2.merge(df1[cols+["atom_id"]].drop_duplicates(), on=cols, how="left")
        c5 = df2["atom_id"].notna()
        atoms_to_replace = df2.loc[c5, ["atom_id_src","atom_id"]].astype(int).values
        atoms_to_add = df2.loc[~df2["atom_id_src"].isin(atoms_to_replace[:,0]), "atom_id_src"].astype(int).values
        
        if c1 & c3:
            for atom_id_src, atom_id in atoms_to_replace:
                out.atom_array[atom_id] = rec.atom_array[atom_id_src].copy()
            out.atom_array+= rec.atom_array[atoms_to_add]
        else:
            for frame in range(len(out.atom_array)):
                for atom_id_src, atom_id in atoms_to_replace:
                    out.atom_array[frame][atom_id] = rec.atom_array[frame][atom_id_src].copy()
            rec_atom_array = [rec.atom_array[i][atoms_to_add] for i in range(len(rec.atom_array))]
            out.atom_array+= biotite.structure.stack(rec_atom_array)
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
        result = biotite.structure.superimpose(
                    fixed=rec.atom_array, 
                    mobile=out.atom_array,
                    atom_mask=atom_mask)
        (fitted_arr, Aff) = result
        out.atom_array = fitted_arr
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
        result = biotite.structure.superimpose_homologs(
                    fixed=rec.atom_array, 
                    mobile=out.atom_array, 
                    min_anchors=min_anchors)
        (fitted_arr, Aff, idx1, idx2) = result
        out.atom_array = fitted_arr
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
        result = biotite.structure.superimpose_structural_homologs(
                    fixed=rec.atom_array, 
                    mobile=out.atom_array)
        (fitted_arr, Aff, idx1, idx2) = result
        out.atom_array = fitted_arr
        return out
        
    def __len__(self):
        """
        Number of atoms
                
        Returns:
            length of atom_array
        """
        return len(self.atom_array)
        
    def __repr__(self):
        """
        Gets representation of underlying data
        
        Returns:
            string representation of the underlying data
        """
        sequences = self.get_sequence()
        sequences = sequences.drop(columns=["sequence"])
        ligands = self.get_ligands()
        
        out = f"List of Sequences:\n{sequences}\n\n"
        out+= f"List of Ligands:\n{ligands}\n\n"
        return out

    def __delitem__(self, key):
        """
        Deletes a list of atoms in the AtomArray
        
        Args:
            key: slice or list of bools or strings
        
        Returns:
            modified copy of self
        """
        out = self.copy()
        # slice through dataframe like a string
        if isinstance(key, slice):
            if isinstance(out.atom_array, biotite.structure.AtomArrayStack):
                x1 = out.atom_array[:, :key.start]
                x2 = out.atom_array[:, key.stop:]
            else:
                x1 = out.atom_array[:key.start]
                x2 = out.atom_array[key.stop:]
            out.atom_array = x1 + x2
        
        # return a slice of atom_stack
        elif isinstance(key, int):
            del(out[key:key+1])
        
        # slice via array indices or bool indices
        elif isinstance(key, np.ndarray):
            if isinstance(out.atom_array, biotite.structure.AtomArrayStack):
                idx2 = np.arange(len(out.atom_array[0]))
                idx1 = idx2[key]
                mask = ~np.isin(idx2, idx1)
                out.atom_array = out.atom_array[:, mask]
            else:
                idx2 = np.arange(len(out.atom_array))
                idx1 = idx2[key]
                mask = ~np.isin(idx2, idx1)
                out.atom_array = out.atom_array[mask]
        
        # return annotation values
        elif isinstance(key, str):
            out.atom_array.del_annotation(key)
        
        # catch errors
        else:
            raise RuntimeError(f"Unable to __delitem__ with {key}")
        
        # rewrite values
        self.atom_array = out.atom_array
        
    def __getitem__(self, key):
        """
        Get slices of an AtomArray
        
        Args:
            key: slice or list of bools or strings
        
        Returns:
            copy of the values of the underlying data
        """
        out = self.copy()
        # slice through dataframe like a string
        if isinstance(key, slice):
            if isinstance(out.atom_array, biotite.structure.AtomArrayStack):
                out.atom_array = out.atom_array[:, key.start:key.stop]
            else:
                out.atom_array = out.atom_array[key.start:key.stop]
            return out
        
        # slice via array indices or bool indices
        elif isinstance(key, np.ndarray):
            if isinstance(out.atom_array, biotite.structure.AtomArrayStack):
                out.atom_array = out.atom_array[:, key]
            else:
                out.atom_array = out.atom_array[key]
            return out
        
        # return a slice of atom_stack
        elif isinstance(key, int):
            if isinstance(out.atom_array, biotite.structure.AtomArrayStack):
                out.atom_array = out.atom_array[:, key]
            else:
                out.atom_array = out.atom_array[key]
            return out
        
        # return annotation values
        elif isinstance(key, str):
            return out.atom_array._annot[key]
        
        # catch errors
        else:
            raise RuntimeError(f"Unable to __getitem__ with {key}")
        
    def __setitem__(self, key:tuple, value:str):
        """
        Sets a feature or sequence in AtomArray to some value similar
        to df.loc[:, "column_name"]
        
        Args:
            key: tuple consisting of (array, str)
                str denotes the annotation key to modify
                array denotes the atom_id indices to modify at
            value: string or array
        """
        if isinstance(self.atom_array, biotite.structure.AtomArrayStack):
            num_atoms = len(self.atom_array[0])
        else:
            num_atoms = len(self.atom_array)
        
        # set all indices to a value
        if isinstance(key, str):
            if isinstance(value, int) or isinstance(value, float) or isinstance(value, str):
                value = [value] * num_atoms
            self.atom_array.set_annotation(key, value)
        
        # set target indices to a value
        elif isinstance(key, tuple) and isinstance(key[1], str):
            cur_val = self.atom_array._annot[key[1]]
            cur_val[key[0]] = value
            self.atom_array.set_annotation(key[1], cur_val)
        
        # catch errors
        else:
            raise RuntimeError(f"Unable to __setitem__ with {key}")            
        
    def __radd__(self, other):
        """
        Performs other + self
        
        Args:
            other: Record object
        
        Returns:
            modified copy of object
        """
        out = self.copy()
        # merge two atom_arrays
        if type(out.atom_array) == type(other.atom_array):
            out.atom_array = other.atom_array + out.atom_array
        
        # merge atom_array into atom_array stacks
        elif isinstance(other.atom_array, biotite.structure.AtomArrayStack):
            atom_stack = []
            for i in range(len(other.atom_array)):
                atom_stack.append(other.atom_array[i] + out.atom_array)
            out.atom_array = biotite.structure.stack(atom_stack)
        elif isinstance(out.atom_array, biotite.structure.AtomArrayStack):
            atom_stack = []
            for i in range(len(out.atom_array)):
                atom_stack.append(other.atom_array + out.atom_array[i])
            out.atom_array = biotite.structure.stack(atom_stack)
        return out
        
    def __add__(self, other):
        """
        Performs self + other

        Args:
            other: Record object
        
        Returns:
            modified copy of object
        """
        out = self.copy()
        # merge two atom_arrays
        if type(out.atom_array) == type(other.atom_array):
            out.atom_array = out.atom_array + other.atom_array
        
        # merge atom_array into atom_array stacks
        elif isinstance(other.atom_array, biotite.structure.AtomArrayStack):
            atom_stack = []
            for i in range(len(other.atom_array)):
                atom_stack.append(out.atom_array + other.atom_array[i])
            out.atom_array = biotite.structure.stack(atom_stack)
        elif isinstance(out.atom_array, biotite.structure.AtomArrayStack):
            atom_stack = []
            for i in range(len(out.atom_array)):
                atom_stack.append(out.atom_array[i] + other.atom_array)
            out.atom_array = biotite.structure.stack(atom_stack)
        return out
        
