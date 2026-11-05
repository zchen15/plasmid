
from typing import List

from Bio.Seq import Seq
from Bio.SeqFeature import (
    CompoundLocation,
    Location,
    SimpleLocation,
)

def flatten_location(location: Location) -> List[SimpleLocation]:
    """Recursively flatten a location into its ordered list of parts.

    A :class:`SimpleLocation` is returned as a single-element list. A
    :class:`CompoundLocation` (including nested compounds) is expanded into
    its constituent :class:`SimpleLocation` objects in order.

    Args:
        location: A :class:`SimpleLocation` or :class:`CompoundLocation`
            (possibly nested).

    Returns:
        A list of :class:`SimpleLocation` objects in the original order.

    Example:
        >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
        >>> loc = CompoundLocation([
        ...     SimpleLocation(0, 100, strand=1),
        ...     SimpleLocation(200, 300, strand=1),
        ... ])
        >>> len(flatten_location(loc))
        2
    """
    parts = getattr(location, "parts", None)
    if parts is None:
        return [location]
    out: List[SimpleLocation] = []
    for p in parts:
        out.extend(flatten_location(p))
    return out

def unwrap_location(location: Location, seq_length: int) -> List[int]:
    """Flatten a location into a ``[start, end, strand]`` list.

    Handles locations that wrap the origin of a circular sequence. For a
    forward-strand wrapped feature such as
    ``join(900..1000, 0..100)`` on a 1000 bp sequence, the returned
    ``end`` will be smaller than ``start`` (i.e. ``[900, 100, 1]``) so
    that the wrapped span equals ``end - start + seq_length``.

    Reverse-strand wrapping is detected by inspecting the order of parts
    and their strand values.

    Args:
        location: A :class:`SimpleLocation` or :class:`CompoundLocation`.
        seq_length: Length of the sequence. Currently used for validation
            and future wrap-aware logic; the returned coordinates are
            already modular.

    Returns:
        A list ``[start, end, strand]`` where ``start`` and ``end`` are
        ``int`` and ``strand`` is ``1`` or ``-1``. If the feature wraps
        the origin, ``end`` will be less than or equal to ``start``.

    Raises:
        ValueError: If ``location`` has no parts or has an unsupported
            type.

    Example:
        >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
        >>> wrapped = CompoundLocation([
        ...     SimpleLocation(900, 1000, strand=1),
        ...     SimpleLocation(0, 100, strand=1),
        ... ])
        >>> unwrap_location(wrapped, 1000)
        [900, 100, 1]
    """
    parts = flatten_location(location)
    if not parts:
        raise ValueError("location must contain at least one part")
    
    strand = parts[0].strand
    if strand is None:
        raise ValueError("location strand must be +1 or -1, got None")
    
    # Single, non-wrapping part: return as-is.
    if len(parts) == 1:
        return [int(parts[0].start), int(parts[0].end), int(strand)]
    
    first, last = parts[0], parts[-1]
    
    if strand == 1:
        # Forward: wrap if the last part begins at 0 or before the first part.
        if last.start == 0 or last.start < first.start:
            return [int(first.start), int(last.end), 1]
        return [int(first.start), int(last.end), 1]
    
    # Reverse strand: the last part is the 3' end in linear coordinates.
    return [int(last.start), int(first.end), -1]

def SeqFeatureLocation_to_Location(loc):
    """
    Convert biopython SeqFeature.Location to Location object
    """
    loc_list = flatten_location(loc)
    new_loc = Location(loc_list[0].start, loc_list[0].end, loc_list[0].strand)
    for loc in loc_list[1:]:
        loc = Location(loc.start, loc.end, loc.strand)
        new_loc+= loc
    return new_loc

class Location:
    """
    This class holds function related to storing feature locations
    """
    def __init__(self, start:int, stop:int, strand:int):
        self.positions = [[start, stop, strand]]
        self.columns = ["start", "end", "strand"]
        
    def to_BioLocation(self):
        """
        Converts location to biopython SeqFeature Location objects
        """
        if len(self.positions) == 1:
            [start, stop, strand] = self.positions[0]
            return SimpleLocation(start, stop, strand)
        
        parts = [SimpleLocation(start, stop, strand)
                 for start, stop, strand in self.positions]
        if len(parts) == 1:
            return parts[0]
        return CompoundLocation(parts)
        
    def copy(self):
        """
        Returns a copy of self
        """
        return copy.deepcopy(self)
        
    def __str__(self):
        """
        Format list of indices into a string
        """
        out = []
        for start, stop, strand in self.positions:
            if strand==1:
                strand = "+"
            elif strand==-1:
                strand = "-"
            else:
                strand = "0"
            out.append([f"[{start}:{stop}]({strand})"])
        out = ", ".join(out)
        return out
        
    def __repr__(self):
        """
        Returns the representation of the string
        """
        return self.__str__()
        
    def merge_pos(self, plist:list):
        if len(plist) == 1:
            return plist
        
        out = []
        [s1, e1, sd1] = plist[0]
        for i in range(1, len(plist)):
            [s2, e2, sd2] = plist[i]
            
            # merge forward stands
            if sd1 == 1 and sd2 == 1 and e1 == s2:
                e1 = e2
            # merge reverse strands
            elif sd1 == -1 and sd2 == -1 and s1 == e2:
                s1 = s2

            # todo fuck introns/exons
            
            # check strands
            if len(plist)-1 <= i and sd1 != sd2:
                out.append([s1, e1, sd1])
                s1 = s2
                e1 = e2
                sd1 = sd2
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
        new_pos = self.merge_pos(out.positions + other.positions)
        out.positions = new_pos
        return out

    def __radd__(self, other):
        """
        Performs other + self
        
        Args:
            other: Record object
        
        Returns:
            modified copy of object
        """
        out = self.__add__(other)
        return out

    def __get_item__(self, key):
        """
        Slice through list of positions
        """
        out = self.copy()
        out.positions = out.positions[key]
        return out
        
    def slice(self, start, stop):
        """
        Slice through list of positions
        """
        if start > stop:
            raise ValueError("start > stop")
        
        out = self.copy()
        new_pos = []
        for pos in out.positions:
            [s1, e1, sd1] = pos
            if s1 < start:
                s1 = start
            if e1 > stop:
                e1 = stop
            new_pos.append([s1, e1, sd1])
        out.positions = new_pos
        return out
        
    def shift(self, shift_len:int, seq_len:int=None):
        """
        Shift positions
        """
        out = self.copy()
        for i in range(len(out.positions)):
            out.positions[i][0] = out.positions[i][0] + shift_len
            out.positions[i][1] = out.positions[i][1] + shift_len
        
        # circular shifts
        if seq_len is not None:
            # modulo
            for i in range(len(out.positions)):
                out.positions[i][0]%= out.positions[i][0]
                out.positions[i][1]%= out.positions[i][1]
            
            [start, stop, strand] = out.positions[0]
            
            new_loc = None
            for start, stop, strand in out.positions:
                if start > stop:
                    loc = Location(start, seq_len, strand) + Location(0, stop, strand)
                else:
                    loc = Location(start, stop, strand)

                if new_loc is None:
                    new_loc = loc
                else:
                    new_loc+= loc
            return new_loc
        else:
            return out
        
    def start(self):
        """
        Get the lowest start positions
        """
        return min([p[0] for p in self.positions])

    def end(self):
        """
        Get the highest stop positions
        """
        return max([p[1] for p in self.positions])

    def strand(self):
        """
        Get the strand of the location.

        Returns:
            int: ``1`` if all entries in ``positions`` have strand ``1``,
            ``-1`` if all entries have strand ``-1``, and ``0`` if the
            strands are mixed or unknown.

        Examples:
            >>> Location(100, 200, 1).strand()
            1
            >>> Location(100, 200, -1).strand()
            -1
        """
        if len(self.positions) == 0:
            return 0
        strands = [entry[2] for entry in self.positions]
        if all(strand == 1 for strand in strands):
            return 1
        if all(strand == -1 for strand in strands):
            return -1
        return 0

    def extract(self, sequence: str) -> str:
        """
        Extracts the sequence based on locations.

        Works like biopython's location extraction: returns the
        subsequence of ``sequence`` spanned by this location, reverse
        complementing it when the strand is negative.

        Args:
            sequence (str): The nucleotide sequence to extract from.

        Returns:
            str: The extracted feature sequence.
        """
        seq = str(sequence)
        parts = []
        for start, stop, strand in self.positions:
            segment = seq[int(start):int(stop)]
            if strand == -1:
                segment = str(Seq(segment).reverse_complement())
            parts.append(segment)
        if len(parts) == 1:
            return parts[0]
        return "".join(parts)