# science libs
import numpy as np
import pandas as pd
import Bio
import Bio.Seq
import Bio.SeqIO
import Bio.SeqFeature
import Bio.SeqRecord

# system libraries
import copy
import json
import time
import datetime
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional, Union

# location helpers (all wrap-aware)
from .locations import Location, SeqFeatureLocation_to_Location

# sequence search helpers
from .search import (
    match_DNA,
    match_AA,
    search_ORF,
)

# Module-level constants
FWD_COLOR = "ApEinfo_fwdcolor"
REV_COLOR = "ApEinfo_revcolor"
LABEL = "ApEinfo_label"
FEATURE_COLS = ["name", "feature", "location", "color", "description"]
QUALIFIER_CLEAR_KEYS = (FWD_COLOR, REV_COLOR, LABEL, "gene", "locus_tag", "label")

DEFAULT_CMAP = {
    "CDS": "green",
    "intron": "turquoise",
    "ncRNA": "turquoise",
    "promoter": "orange",
    "enhancer": "salmon",
    "RBS": "cyan",
    "protein_bind": "salmon",
    "repeat_region": "turquoise",
    "polyA": "red",
    "rep_origin": "blue",
}
DEFAULT_COLOR = "purple"

# ---------------------------------------------------------------------------
# Helpers (sequence/annotation; location logic lives in locations.py)
# ---------------------------------------------------------------------------

def update_dict(
    new: Optional[Mapping[str, Any]] = None,
    old: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Merge two mappings, skipping ``None`` values in ``new``.

    This behaves like ``dict.update`` except that entries in ``new`` whose
    value is ``None`` do not overwrite existing values in ``old``. That
    makes it convenient for applying optional parameters, where ``None``
    means "keep the existing value".

    Neither input is modified; a new dictionary is returned.

    Args:
        new: Values to apply. ``None`` values are skipped. May be ``None``.
        old: Base values. May be ``None``, in which case an empty mapping
            is used.

    Returns:
        A new dictionary containing ``old`` updated with the non-``None``
        entries from ``new``.

    Examples:
        >>> update_dict({"a": 1, "b": None}, {"a": 0, "c": 3})
        {'a': 1, 'c': 3}
        >>> update_dict(None, {"a": 0})
        {'a': 0}
    """
    out: dict[str, Any] = dict(old) if old is not None else {}
    if new:
        out.update({key: value for key, value in new.items() if value is not None})
    return out

def isArrayType(x, dtype) -> bool:
    """Check whether every element of ``x`` is an instance of ``dtype``.

    Strings, bytes, and bytearrays are *not* considered arrays even though
    they are iterable — they're treated as scalar values.

    Args:
        x: Value to inspect. Non-iterables (and strings/bytes) return ``False``.
        dtype: A type or tuple of types accepted by :func:`isinstance`.

    Returns:
        ``True`` if ``x`` is a non-string iterable and every element matches
        ``dtype`` (vacuously ``True`` for an empty iterable), otherwise ``False``.
    """
    # Strings/bytes are iterable but should not be treated as arrays.
    if isinstance(x, (str, bytes, bytearray)):
        return False

    try:
        items = iter(x)
    except TypeError:
        return False

    return all(
        isinstance(item, dtype)
        or (isinstance(item, np.generic) and isinstance(item.item(), dtype))
        for item in items
    )

def search_key(data: dict, keyword: str, rev: bool = False) -> Optional[Any]:
    """Return the value of the first key that matches ``keyword``.

    Matching is case-insensitive substring matching, evaluated in the
    dictionary's insertion order. Non-string keys are ignored.

    Args:
        data: Mapping to search.
        keyword: Text to match. Converted with :class:`str` before
            comparison.
        rev: Reverse the containment test. When ``False`` (default), a key
            matches if it contains ``keyword``; when ``True``, a key
            matches if it is contained in ``keyword``.

    Returns:
        The value associated with the first matching key, or ``None`` if
        no key matches.

    Examples:
        >>> search_key({"CDS": "green", "promoter": "orange"}, "cds")
        'green'
        >>> search_key({"gfp": "green"}, "GFP_1", rev=True)
        'green'
    """
    needle = str(keyword).lower()
    for key, value in data.items():
        if not isinstance(key, str):
            continue
        haystack = key.lower()
        if rev:
            matched = haystack in needle
        else:
            matched = needle in haystack
        if matched:
            return value
    return None


def read_genbank(fname: str) -> "Record":
    """Read a GenBank file into a :class:`Record`.

    Args:
        fname: Path to the GenBank (``.gb`` / ``.gbk``) file to read. The
            same value is used as the returned record's ``name``.

    Returns:
        A :class:`Record` containing the parsed sequence, features, and
        annotations.

    Raises:
        ValueError: If the file does not contain exactly one GenBank
            record.
    """
    rec = Bio.SeqIO.read(fname, "genbank")
    return Record(rec, name=fname)


def read_snapgene(fname: str) -> "Record":
    """Read a SnapGene file into a :class:`Record`.

    Args:
        fname: Path to the SnapGene (``.dna``) file to read. The same
            value is used as the returned record's ``name``.

    Returns:
        A :class:`Record` containing the parsed sequence, features, and
        annotations.

    Raises:
        ValueError: If the file does not contain exactly one SnapGene
            record.
    """
    rec = Bio.SeqIO.read(fname, "snapgene")
    return Record(rec, name=fname)


def get_default_colors(name: str, feature: str, color: Optional[list] = None) -> list:
    """Get a set of default colors based on feature type.

    Args:
        name: Name/label of genomic feature (unused; kept for API compat).
        feature: Type of genomic feature (e.g., ``'CDS'``, ``'promoter'``).
        color: Current color value or None to use defaults.

    Returns:
        List ``[fwd_color, rev_color]``.
    """
    if color is not None:
        if isinstance(color, str):
            return [color, color]
        if len(color) == 1:
            return [color[0], color[0]]
        return [color[0], color[-1]]

    found = search_key(DEFAULT_CMAP, feature)
    if found is None:
        return [DEFAULT_COLOR, DEFAULT_COLOR]
    if isinstance(found, str):
        return [found, found]
    if len(found) == 1:
        return [found[0], found[0]]
    return [found[0], found[-1]]


def get_feature_name(gene: Bio.SeqFeature.SeqFeature) -> str:
    """Get name of gene from a Bio.SeqFeature (checks known qualifier tags)."""
    for tag in ("ApEinfo_label", "gene", "locus_tag", "label"):
        if tag in gene.qualifiers:
            val = gene.qualifiers[tag]
            if isinstance(val, list) and len(val) > 0:
                return val[0]
            if isinstance(val, str):
                return val
    return "unknown"


def SeqRecord_to_dataframe(rec: Bio.SeqRecord.SeqRecord) -> pd.DataFrame:
    """Convert GenBank record features to a pandas DataFrame.

    Does not mutate the input SeqRecord.

    Args:
        rec: A :class:`Bio.SeqRecord.SeqRecord`.

    Returns:
        DataFrame with columns matching ``FEATURE_COLS``.
    """
    data = []
    if isinstance(rec, Bio.SeqRecord.SeqRecord):
        for gene in rec.features:
            name = get_feature_name(gene)
            feature = gene.type
            loc = SeqFeatureLocation_to_Location(gene.location)
            
            q = gene.qualifiers
            if FWD_COLOR in q and REV_COLOR in q:
                color = [q[FWD_COLOR][0], q[REV_COLOR][0]]
            else:
                color = get_default_colors(name, feature)
            
            description = {k: v for k, v in q.items() if k not in QUALIFIER_CLEAR_KEYS}
            data.append([name, feature, loc, color, description])

    return pd.DataFrame(data, columns=FEATURE_COLS)


def dataframe_to_SeqRecord(df: pd.DataFrame, seq: str) -> Bio.SeqRecord.SeqRecord:
    """Convert features DataFrame + nucleotide sequence to a Bio.SeqRecord.

    Args:
        df: DataFrame with columns matching ``FEATURE_COLS``.
        seq: Nucleotide sequence.

    Returns:
        :class:`Bio.SeqRecord.SeqRecord`.
    """
    SeqRecord = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(seq))

    features = []
    for name, feature, loc, color, qualifiers in df[FEATURE_COLS].values:
        quals = dict(qualifiers) if qualifiers else {}
        quals["label"] = name
        quals[FWD_COLOR] = [color[0]]
        quals[REV_COLOR] = [color[-1]]
        gene = Bio.SeqFeature.SeqFeature(type=feature,
                                         location=loc.to_BioLocation(),
                                         qualifiers=quals)
        features.append(gene)
    SeqRecord.features = features
    return SeqRecord

def format_to_SeqRecord(
    data: Union[str, "Record", Bio.SeqRecord.SeqRecord],
) -> Bio.SeqRecord.SeqRecord:
    """Format a string, Record, or SeqRecord into a SeqRecord.

    Args:
        data: A string, ``Record``, or ``Bio.SeqRecord.SeqRecord``.

    Returns:
        A :class:`Bio.SeqRecord.SeqRecord`.

    Raises:
        TypeError: If ``data`` cannot be converted.
    """
    if isinstance(data, Bio.SeqRecord.SeqRecord):
        return data
    if isinstance(data, Record):
        return dataframe_to_SeqRecord(data.features, data.sequence)
    if isinstance(data, str):
        return Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(data))
    raise TypeError(f"Cannot convert {type(data)!r} to SeqRecord")


# ---------------------------------------------------------------------------
# Record
# ---------------------------------------------------------------------------

class Record:
    """Holds functions pertaining to viewing and editing GenBank files.

    Location manipulation (wrapping, shifting, merging) is delegated to
    ``locations.py``.

    Attributes:
        annotations: GenBank metadata.
        sequence: Nucleotide sequence.
        features: Feature table with columns
            ``[name, feature, location, color, description]``.
        name: Identifier for the record (usually a filename).
    """

    def __init__(self, rec: Union[str, Bio.SeqRecord.SeqRecord], name: str = "."):
        """Initialize the Record object.

        Args:
            rec: ``Bio.SeqRecord`` object or raw nucleotide string.
            name: Name of the plasmid.

        Raises:
            TypeError: If ``rec`` is neither ``str`` nor ``SeqRecord``.
        """
        self.annotations: dict = {
            "molecule_type": "DNA",
            "data_file_division": "circular",
        }
        
        if isinstance(rec, Bio.SeqRecord.SeqRecord):
            self.features = SeqRecord_to_dataframe(rec)
            self.sequence = str(rec.seq)
            ref = rec.annotations
            for k in ("data_file_division", "date", "source", "organism",
                      "taxonomy", "molecule_type"):
                if k in ref:
                    self.annotations[k] = ref[k]
        elif isinstance(rec, str):
            self.sequence = rec
            self.features = pd.DataFrame(columns=FEATURE_COLS)
        else:
            raise TypeError(f"Record expects str or SeqRecord, got {type(rec)!r}")

        self.name = name
        self.reset_index()

    # -- internal helpers ---------------------------------------------------

    def reset_index(self) -> None:
        """Sort features in place by ``(start, length, feature, name)``.

        Returns:
            None. The record's feature table is modified in place.
        """        
        df = self.features
        if len(df) == 0:
            return self
        
        starts = np.array([v.start() for v in df["location"].values])
        ends = np.array([v.end() for v in df["location"].values])
        
        order = pd.DataFrame({
            "start": starts,
            "length": ends - starts,
            "feature": df["feature"].values,
            "name": df["name"].values,
            "_idx": np.arange(len(df)),
        })
        order = order.sort_values(by=["start", "length", "feature", "name"])
        self.features = df.iloc[order["_idx"].values].reset_index(drop=True)
        
        return self
        
    # -- I/O ----------------------------------------------------------------

    def to_genbank(self, filename: str) -> None:
        """Write the record to a GenBank file.

        The first whitespace-delimited token of the output file's stem is
        used as the GenBank record name, and the current date is written
        into the annotations.

        Args:
            filename: Destination path for the GenBank file.

        Returns:
            None. The file is written to disk.
        """
        SeqRecord = dataframe_to_SeqRecord(self.features, self.sequence)
        SeqRecord.annotations = update_dict(self.annotations, SeqRecord.annotations)

        fname = Path(filename).stem.split(" ")[0]
        SeqRecord.name = fname
        SeqRecord.annotations["date"] = (
            datetime.datetime.fromtimestamp(time.time()).strftime("%d-%b-%Y").upper()
        )
        with open(filename, "w") as f:
            Bio.SeqIO.write(SeqRecord, f, "genbank")

    def to_dataframe(self, str_only: bool = True) -> pd.DataFrame:
        """Generate a DataFrame of features with resolved start/end/strand.

        Args:
            str_only: When ``True`` (default), ``location`` and
                ``description`` are serialized as strings. When ``False``,
                ``location`` holds :class:`Location` objects and
                ``description`` holds the qualifier dictionaries.

        Returns:
            Copy of ``self.features`` with added columns ``start``, ``end``,
            ``strand``, ``sequence``, ``length``.
        """
        df = self.features.copy()
        if len(df) == 0:
            for col in ("start", "end", "strand", "sequence", "length"):
                df[col] = []
            return df
        
        pLlen = len(self.sequence)
        df["start"] = [v.start() for v in df["location"].values]
        df["end"] = [v.end() for v in df["location"].values]
        df["strand"] = [v.strand() for v in df["location"].values]
        df["sequence"] = [loc.extract(self.sequence) for loc in df["location"].values]
        df["length"] = [len(s) for s in df["sequence"]]
        if str_only:
            df["location"] = [str(loc) for loc in df["location"].values]
            df["description"] = [json.dumps(d, indent=4) for d in df["description"].values]
        return df

    # -- copies -------------------------------------------------------------

    def copy(self) -> "Record":
        """Return a deep copy of self."""
        return copy.deepcopy(self)

    # -- annotation ---------------------------------------------------------

    def annotate(
        self,
        name: str,
        seq: str,
        feature: str = "unknown",
        color: Optional[list] = None,
        description: Optional[dict] = None,
        circular: bool = True,
    ) -> "Record":
        """Add gene annotations via DNA substring search.

        Args:
            name: Gene name.
            seq: DNA sequence to search for.
            feature: Feature type (e.g. ``'CDS'``, ``'primer_bind'``).
            color: ``[fwd, rev]`` colors or None for defaults.
            description: Extra qualifiers.
            circular: Treat the sequence as circular.

        Returns:
            A copy of ``self`` with the new feature(s) added.
        """
        seq = seq.replace(" ", "")
        if seq == "":
            return self

        if description is None:
            description = {}
        if color is None:
            color = get_default_colors(name, feature)

        out = self.copy()
        seq_u = seq.upper().replace("U", "T")
        pos = match_DNA(seq_u, out.sequence, circular=circular)
        return out.annotate_from_positions(pos, name, feature, color, description)

    def annotate_AA(
        self,
        name: str,
        seq: str,
        feature: str = "unknown",
        color: Optional[list] = None,
        description: Optional[dict] = None,
        circular: bool = True,
    ) -> "Record":
        """Add gene annotations via amino-acid substring search.

        Args:
            name: Gene name.
            seq: Amino-acid sequence to search for.
            feature: Feature type.
            color: ``[fwd, rev]`` colors or None for defaults.
            description: Extra qualifiers.
            circular: Treat the sequence as circular.

        Returns:
            A copy of ``self`` with the new feature(s) added.
        """
        seq = seq.replace(" ", "")
        if seq == "":
            return self

        if description is None:
            description = {}
        if color is None:
            color = get_default_colors(name, feature)

        out = self.copy()
        pos = match_AA(seq, out.sequence, circular=circular)
        return out.annotate_from_positions(pos, name, feature, color, description)
        
    def annotate_from_positions(
        self,
        pos: pd.DataFrame,
        name: str,
        feature: str,
        color: list,
        description: dict,
    ) -> "Record":
        """Append a feature for each row of a match table; mutate in place.

        Args:
            pos: DataFrame of matches with columns ``t_start``, ``t_end``,
                and ``orientation`` (as produced by :func:`match_DNA` or
                :func:`match_AA`).
            name: Feature name/label.
            feature: Feature type (e.g. ``'CDS'``, ``'promoter'``).
            color: ``[forward_color, reverse_color]`` list.
            description: Extra qualifier key/value pairs.

        Returns:
            ``self``, mutated to include the new feature(s).
        """
        if len(pos) == 0:
            return self
        pLlen = len(self.sequence)
        glist = []
        for start, end, strand in pos[["t_start", "t_end", "orientation"]].values:
            loc = Location(start, end, strand)
            glist.append([name, feature, loc, color, dict(description)])
        new_df = pd.DataFrame(glist, columns=FEATURE_COLS)
        self.features = pd.concat([self.features, new_df], ignore_index=True)
        self.reset_index()
        return self
        
    def annotate_by_location(
        self,
        name: str,
        start: int,
        stop: int,
        strand: int,
        feature: str = "unknown",
        color: Optional[list] = None,
        description: Optional[dict] = None,
    ) -> "Record":
        """Add annotations using location information.
    
        Args:
            name: Name of the gene.
            start: Start location of the feature (0-based).
            stop: Stop location of the feature (exclusive).
            strand: Orientation of the feature. ``+1`` is 5'→3', ``-1`` is 3'→5'.
            feature: Type of genetic feature (e.g. ``'CDS'``, ``'mRNA'``,
                ``'primer_bind'``).
            color: ``[fwd_color, rev_color]`` list, or ``None`` to use defaults.
            description: Dictionary of extra qualifiers. ``None`` is treated as
                an empty dict.
    
        Returns:
            A modified copy of ``self``.
        """
        out = self.copy()
        if description is None:
            description = {}
        if color is None:
            color = get_default_colors(name, feature)
        
        loc = Location(start, stop, strand)
        new_row = pd.DataFrame(
            [[name, feature, loc, color, dict(description)]],
            columns=FEATURE_COLS,
        )
    
        out.features = pd.concat([out.features, new_row], ignore_index=True)
        out.reset_index()
        return out

    def get_ORFs(self, thresh: int = 50) -> "Record":
        """Return a copy containing only ORF features above ``thresh`` amino acids."""
        out = self.copy()
        df = search_ORF(out.sequence)
        df = df[df["length"] > thresh]

        out.features = pd.DataFrame(columns=FEATURE_COLS)
        for i in range(len(df)):
            name = f"ORF_{i}"
            seq = df["sequence"].values[i]
            out = out.annotate_AA(name, seq=seq, feature="CDS")
        return out

    # -- dunders ------------------------------------------------------------

    def __len__(self) -> int:
        """Return the number of features in the record.

        Returns:
            The feature count (not the sequence length).
        """
        return len(self.features)

    def __str__(self) -> str:
        return self.sequence

    def get_colored(self, width=0):
        '''
        Print out colored representation of Plasmid
        width = number of characters per line
        '''
        from .terminal import Viewer
        rec = Viewer(self)
        out = rec.colorize(width=width)
        return out

    def __repr__(self) -> str:
        """Return a human-readable multi-line summary of the record.

        The summary includes the object's class and memory address, the
        annotation dictionary, the filename, the sequence length, and a
        table of features (name, type, location, color, and length).

        Returns:
            A multi-line string suitable for printing or logging.
        """
        lines = [
            f"{type(self).__name__} at {hex(id(self))}",
            str(self.annotations),
            f"filename = {self.name}",
            f"length = {len(self.sequence)}",
        ]
        cols = ["name", "feature", "location", "color"]
        df = self.features[cols].copy() if len(self.features) else pd.DataFrame(columns=cols)
        if len(df):
            x = self.to_dataframe()
            df["length"] = x["sequence"].str.len().values
            df["location"] = [str(v) for v in df["location"].values]
        lines.append(str(df))
        return "\n".join(lines)

    # -- slicing ------------------------------------------------------------

    def slice(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> "Record":
        """Return a sub-Record containing only features fully inside ``[start, end)``.

        Args:
            start: Start base (0-based). ``None`` means 0.
            end: Stop base (exclusive). ``None`` means
                ``len(self.sequence)``.

        Returns:
            A copy of ``self`` sliced to ``[start, end)``.
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self.sequence)

        out = self.copy()
        
        L = len(out.sequence)
        if len(out.features) == 0:
            out.sequence = out.sequence[start:end]
            return out
        
        keep: list[bool] = []
        loc_list = []
        for loc in out.features["location"].values:
            loc_start = loc.start()
            loc_end = loc.end()
            if loc_start >= start and loc_end <= end:
                loc_list.append(loc.shift(L, circular=True))
                keep.append(True)
            else:
                loc_list.append(loc)
                keep.append(False)
        out.features["location"] = loc_list
        out.features = out.features[keep]
        out.sequence = out.sequence[start:end]
        out.reset_index()
        return out

    def __getitem__(self, key):
        """Index the record by feature index, name, column, or slice.

        Supported key forms:

        * ``slice`` — returns a sub-:class:`Record` (see :meth:`slice`).
        * ``int`` — returns a sub-:class:`Record` containing the single
          feature at that positional index, trimmed to that feature's span.
        * ``str`` (a feature-table column name) — returns that column as a
          ``pandas.Series``.
        * ``str`` (any other name) — returns that column from
          :meth:`to_dataframe`.
        * ``(index, column_name)`` tuple — returns the scalar cell.
        * Boolean ``numpy.ndarray`` — returns a copy with the selected
          features.
        * Integer array — returns a copy with the selected features.

        Args:
            key: Index, slice, column name, ``(index, column)`` tuple,
                boolean mask, or integer array.

        Returns:
            A :class:`Record`, a ``pandas.Series``, or a scalar, depending
            on the key type.

        Raises:
            TypeError: If ``key`` is of an unsupported type.
        """
        if isinstance(key, slice):
            return self.slice(key.start, key.stop)
        
        if isinstance(key, (int, np.integer)) and not isinstance(key, bool):
            loc = self.features.iloc[key]["location"]
            return self.slice(loc.start(), loc.end())
        
        if isinstance(key, str) or isArrayType(key, str):
            if key in self.features.columns:
                return self.features[key]
            return self.to_dataframe()[key]
        
        if isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], str):
            row, column = key
            return self.features.iloc[row][column]

        if isinstance(key, np.ndarray) and key.dtype == bool:
            out = self.copy()
            out.features = out.features.loc[key]
            out.reset_index()
            return out

        if isArrayType(key, int):
            out = self.copy()
            out.features = out.features.iloc[key]
            out.reset_index()
            return out

        raise TypeError(f"Unsupported key type: {type(key)!r}")

    def __setitem__(self, key, value) -> None:
        """Assign a new sequence (or cell) into the record.

        Supported key forms:

        * ``int`` — replace the span of the feature at that index with
          ``value`` (a sequence string or :class:`Record`).
        * ``slice`` — replace the sliced span with ``value``.
        * ``(index, column_name)`` tuple — set that cell in the feature
          table.

        Mutates ``self.sequence`` / ``self.features`` and re-sorts features.

        Args:
            key: Feature index, slice, or ``(index, column)`` tuple.
            value: Replacement sequence, :class:`Record`, or cell value.

        Returns:
            None. The record is modified in place.

        Raises:
            TypeError: If ``key`` is of an unsupported type.
        """
        if isinstance(key, (int, np.integer)) and not isinstance(key, bool):
            loc = self.features.iloc[key]["location"]
            out = self.replace(loc.start(), loc.end(), value)
            self.sequence = out.sequence
            self.features = out.features
        elif isinstance(key, slice):
            out = self.replace(key.start, key.stop, value)
            self.sequence = out.sequence
            self.features = out.features
        elif isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], str):
            row, column = key
            self.features.iloc[row, self.features.columns.get_loc(column)] = value
        else:
            raise TypeError(f"Unsupported key type: {type(key)!r}")
        self.reset_index()

    def __delitem__(self, key) -> None:
        """Delete feature(s) or a sequence span from the record.

        Supported key forms:

        * ``int`` — drop the feature at that index.
        * Boolean ``numpy.ndarray`` — drop the masked features.
        * Integer array — drop the selected features.
        * ``slice`` — delete the sliced sequence span (features fully
          inside the span are removed via :meth:`replace`).

        Mutates ``self.features`` / ``self.sequence`` and re-sorts features.

        Args:
            key: Feature index, boolean mask, integer array, or slice.

        Returns:
            None. The record is modified in place.

        Raises:
            TypeError: If ``key`` is of an unsupported type.
        """
        if isinstance(key, (int, np.integer)) and not isinstance(key, bool):
            self.features = self.features.drop(self.features.index[int(key)])
        elif isinstance(key, np.ndarray) and key.dtype == bool:
            self.features = self.features.loc[~key]
        elif isArrayType(key, int):
            L = len(self.features)
            keep = np.ones(L, dtype=bool)
            keep[key] = False
            self.features = self.features.loc[keep]
        elif isinstance(key, slice):
            out = self.replace(key.start, key.stop, "")
            self.sequence = out.sequence
            self.features = out.features
        else:
            raise TypeError(f"Unsupported key type: {type(key)!r}")
        self.reset_index()

    # -- replace ------------------------------------------------------------

    def replace(self, start, end, value, circular: bool = True) -> "Record":
        """Replace the sequence span ``[start, end)`` with ``value``.

        Internal helper used by :meth:`__setitem__`, :meth:`__delitem__`,
        and the public ``replace`` API. When ``circular`` is true, the
        record is temporarily re-origined at ``start`` so the span (which
        may wrap the origin) can be replaced cleanly, then re-origined back.

        Args:
            start: Start base (0-based); ``None`` means 0.
            end: Stop base (exclusive); ``None`` means
                ``len(self.sequence)``. May be less than ``start`` to
                indicate a wrap-around span when ``circular`` is true.
            value: Replacement sequence string or :class:`Record`.
            circular: Treat the sequence as circular. Defaults to ``True``.

        Returns:
            A new :class:`Record` with the span replaced. The original
            ``name`` is preserved.
        """
        out = self.copy()
        if start is None:
            start = 0
        if end is None:
            end = len(self.sequence)

        if circular:
            out = out.set_origin(start)
            L = end - start
            if end < start:
                L = end - start + len(self.sequence)
            out = value + out[L:]
            out = out.set_origin(-start)
        else:
            out = out[:start] + value + out[end:]
        out.name = self.name
        return out

    # -- concatenation ------------------------------------------------------

    def __radd__(self, other) -> "Record":
        """Prepend ``other`` to this record: ``other + self``.

        Coerces ``other`` via :func:`format_to_SeqRecord`, concatenates it
        before this record's sequence and features, and merges annotations
        (and, when ``other`` is a :class:`Record`, its name).

        Args:
            other: A string, :class:`Record`, or
                :class:`Bio.SeqRecord.SeqRecord` to prepend.

        Returns:
            A new :class:`Record` whose sequence begins with ``other``.
        """
        rec1 = dataframe_to_SeqRecord(self.features, self.sequence)
        rec2 = format_to_SeqRecord(other)
        rec = Record(rec2 + rec1)
        if isinstance(other, Record):
            rec.annotations = update_dict(other.annotations, dict(self.annotations))
            rec.name = other.name
        else:
            rec.annotations = update_dict(self.annotations)
            rec.name = self.name
        return rec

    def __add__(self, other) -> "Record":
        """Append ``other`` to this record: ``self + other``.

        Coerces ``other`` via :func:`format_to_SeqRecord`, concatenates it
        after this record's sequence and features, and carries over this
        record's annotations and name.

        Args:
            other: A string, :class:`Record`, or
                :class:`Bio.SeqRecord.SeqRecord` to append.

        Returns:
            A new :class:`Record` whose sequence ends with ``other``.
        """
        rec1 = dataframe_to_SeqRecord(self.features, self.sequence)
        rec2 = format_to_SeqRecord(other)
        rec = Record(rec1 + rec2)
        rec.annotations = update_dict(self.annotations)
        rec.name = self.name
        return rec

    # -- transformations ----------------------------------------------------

    def set_origin(self, index:int):
        """
        Sets a new position as the origin
        
        Args:
            index: nucleotide position of the new origin
        
        Returns:
            a copy of the object
        """
        out = self.copy()
        
        slen = len(self.sequence)
        index = index%slen
        
        out.sequence = out.sequence[index:] + out.sequence[0:index]
        out.features["location"] = [loc.shift(-index, slen) for loc in out.features["location"]]
        out.reset_index()
        return out

    # -- sequence utilities -------------------------------------------------

    def translate(self, table: str = "Standard", frame: int = 0) -> str:
        """Translate the nucleotide sequence to amino acids.

        Args:
            table: Codon table to use (default ``"Standard"``).
            frame: Reading frame offset (0, 1, or 2). Default 0.

        Returns:
            Amino-acid string (stops as ``'*'``).
        """
        seq = Bio.Seq.Seq(self.sequence)
        return str(seq[frame:].translate(table=table))

    def revcomp(self) -> "Record":
        """Return a reverse-complemented copy of the record.

        The sequence is reverse-complemented and every feature location is
        flipped to the opposite strand.  Feature names, colors, and
        descriptions are preserved.

        Returns:
            A new :class:`Record` with the reverse-complemented sequence and
            features.
        """
        rec = dataframe_to_SeqRecord(self.features, self.sequence)
        rec = rec.reverse_complement()
        out = Record(rec, name=self.name)
        out.annotations = update_dict(self.annotations)
        out.reset_index()
        return out

    # -- duplicate handling -------------------------------------------------

    def duplicated(self) -> pd.Series:
        """Return a boolean Series marking features that are duplicates.

        Two features are considered duplicates when they have the same
        ``feature`` type, ``location`` (resolved to ``[start, end, strand]``),
        and extracted ``sequence``.  The first occurrence is ``False``;
        subsequent occurrences are ``True``.

        Returns:
            Boolean ``pandas.Series`` aligned with the current feature index.
        """
        df = self.to_dataframe()
        cols = ["feature", "location", "sequence"]
        return df[cols].duplicated()

    def drop_duplicates(self) -> "Record":
        """Return a copy with duplicate features removed.

        See :meth:`duplicated` for the definition of a duplicate.  The first
        occurrence of each group is kept.

        Returns:
            A new :class:`Record` without duplicate features.
        """
        out = self.copy()
        mask = out.duplicated()
        out.features = out.features.loc[~mask.values].reset_index(drop=True)
        out.reset_index()
        return out

    # -- merge annotations from another record ------------------------------

    def merge(self, df):
        """
        Add annotations from one genbank file into another
        """
        out = self.copy()
        ref = df.to_dataframe()
        cols = ["name","feature","sequence","color"]
        for name, feature, seq, color in ref[cols].values:
            out = out.annotate(name, seq, feature, color)
        out.reset_index()
        return out