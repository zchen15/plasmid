import copy
from typing import List, Union

from Bio.Seq import Seq
from Bio.SeqFeature import (
    CompoundLocation,
    SimpleLocation,
)

def flatten_location(location: Union[SimpleLocation, CompoundLocation]) -> List[SimpleLocation]:
    """Recursively flatten a location into its ordered list of parts.

    A :class:`SimpleLocation` is returned as a single-element list. A
    :class:`CompoundLocation` (including nested compounds) is expanded into
    its constituent :class:`SimpleLocation` objects in order.

    Args:
        location: A :class:`SimpleLocation` or :class:`CompoundLocation`
            (possibly nested).

    Returns:
        A list of :class:`SimpleLocation` objects in the original order.
    """
    parts = getattr(location, "parts", None)
    if parts is None:
        return [location]
    out: List[SimpleLocation] = []
    for p in parts:
        out.extend(flatten_location(p))
    return out


def merge_pos(plist: list) -> list:
    """Merges adjacent positions on the same strand.

    Takes a list of position triplets and merges consecutive positions that are
    exactly adjacent on the same strand. For forward strand positions, merging
    occurs when the end of one position equals the start of the next (e1 == s2).
    For reverse strand positions, merging occurs when the start of one position
    equals the end of the next (s1 == e2).

    Importantly, this function preserves the original order of the input list
    and does NOT merge overlapping intervals.

    Args:
        plist: A list of position triplets, where each position is a list or
            tuple of [start, stop, strand]. Start and stop are integer
            coordinates, and strand is 1 (forward), -1 (reverse), or 0 (unknown).

    Returns:
        A new list of merged position triplets. If no merging is possible,
        returns a copy of the input list. Returns an empty list if the input
        is empty.
    """
    if not plist:
        return []
    
    out = []
    s1, e1, sd1 = plist[0]
    
    for i in range(1, len(plist)):
        s2, e2, sd2 = plist[i]
        can_merge = False
        
        # merge forward strands (adjacent only)
        if sd1 == 1 and sd2 == 1 and e1 == s2:
            e1 = e2
            can_merge = True
        # merge reverse strands (adjacent only)
        elif sd1 == -1 and sd2 == -1 and s1 == e2:
            s1 = s2
            can_merge = True

        # If we couldn't merge, push the current segment and reset
        if not can_merge:
            out.append([s1, e1, sd1])
            s1, e1, sd1 = s2, e2, sd2
            
    # Append the final accumulated segment
    out.append([s1, e1, sd1])
    return out


def SeqFeatureLocation_to_Location(loc: Union[SimpleLocation, CompoundLocation]) -> "Location":
    """Convert a Biopython SeqFeature location to a Location object.

    This function takes a Biopython :class:`SimpleLocation` or
    :class:`CompoundLocation` and converts it to a custom :class:`Location`
    object. The conversion is optimized to extract all positions first and
    merge them in a single pass, avoiding the overhead of creating and
    combining multiple Location objects.

    Args:
        loc: A Biopython :class:`SimpleLocation` or :class:`CompoundLocation`
            object (possibly nested).

    Returns:
        A :class:`Location` object with merged positions.

    Raises:
        ValueError: If the input location is empty (has no parts).
    """
    loc_list = flatten_location(loc)
    if not loc_list:
        raise ValueError("Cannot convert an empty location.")
        
    # Extract all positions into a single list
    positions = [[part.start, part.end, part.strand] for part in loc_list]
    
    # Merge them all at once in a single O(N) pass
    merged_positions = merge_pos(positions)
    
    # Create a single Location object and assign the merged positions
    first = merged_positions[0]
    new_loc = Location(first[0], first[1], first[2])
    new_loc.positions = merged_positions
    
    return new_loc


class Location:
    """A class for storing and manipulating genomic feature locations.

    This class provides a unified interface for working with feature locations,
    supporting both simple (single interval) and compound (multiple interval)
    locations. It can convert to and from Biopython's SeqFeature location
    objects, perform operations like merging, truncating, and shifting
    coordinates, and extract sequences from genomic DNA.

    Locations are stored internally as a list of position triplets, where each
    triplet contains [start, stop, strand]. The start and stop are integer
    coordinates (0-based, half-open), and strand is 1 (forward), -1 (reverse),
    or 0 (unknown).

    Attributes:
        positions: A list of position triplets [start, stop, strand].
    """
    
    def __init__(self, start: int, stop: int, strand: int):
        """Initialize a Location with a single position.

        Args:
            start: The start coordinate (0-based, inclusive).
            stop: The stop coordinate (0-based, exclusive).
            strand: The strand (1 for forward, -1 for reverse, 0 for unknown).
        """
        self.positions = [[start, stop, strand]]
        
    def to_BioLocation(self) -> Union[SimpleLocation, CompoundLocation, None]:
        """Convert this Location to a Biopython SeqFeature location object.

        Returns a :class:`SimpleLocation` if there is only one position, or a
        :class:`CompoundLocation` if there are multiple positions.

        Returns:
            A Biopython :class:`SimpleLocation` or :class:`CompoundLocation`
            object, or None if the location has no positions.
        """
        if not self.positions:
            return None
            
        parts = [SimpleLocation(start, stop, strand) for start, stop, strand in self.positions]
        if len(parts) == 1:
            return parts[0]
        return CompoundLocation(parts)
        
    def copy(self) -> "Location":
        """Return a deep copy of this Location.

        Returns:
            A new :class:`Location` object with the same positions.
        """
        return copy.deepcopy(self)
        
    def __str__(self) -> str:
        """Format the location as a human-readable string.

        Returns:
            A string representation showing each position as [start:stop](strand),
            separated by commas. Strand is shown as '+' for forward, '-' for
            reverse, or '0' for unknown.
        """
        out = []
        for start, stop, strand in self.positions:
            if strand == 1:
                s = "+"
            elif strand == -1:
                s = "-"
            else:
                s = "0"
            out.append(f"[{start}:{stop}]({s})")
        return ", ".join(out)
        
    def __repr__(self) -> str:
        """Return the string representation of the Location.

        Returns:
            The same string as :meth:`__str__`.
        """
        return self.__str__()
        
    def __add__(self, other: "Location") -> "Location":
        """Combine two locations by concatenating and merging their positions.

        Creates a new Location containing all positions from both locations,
        with adjacent positions on the same strand merged together.

        Args:
            other: Another :class:`Location` object to add.

        Returns:
            A new :class:`Location` object with merged positions.
        """
        out = self.copy()
        new_pos = merge_pos(out.positions + other.positions)
        out.positions = new_pos
        return out

    def __radd__(self, other: "Location") -> "Location":
        """Support for adding from the right side.

        This method is called when the left operand doesn't support addition
        with the right operand. It also handles the special case of adding 0,
        which allows using Python's built-in :func:`sum()` function on a list
        of Location objects.

        Args:
            other: Another :class:`Location` object, or 0.

        Returns:
            A new :class:`Location` object with merged positions.
        """
        if other == 0:
            return self.copy()
        return self.__add__(other)

    def __getitem__(self, key):
        """Slice or index into the list of positions.

        Returns a new Location containing the selected positions. If an integer
        index is provided, returns a Location with a single position. If a slice
        is provided, returns a Location with the sliced positions.

        Args:
            key: An integer index or slice object.

        Returns:
            A new :class:`Location` object with the selected positions.
        """
        out = self.copy()
        # Ensure we always return a 2D list, even if key is an integer
        if isinstance(key, int):
            out.positions = [out.positions[key]]
        else:
            out.positions = out.positions[key]
        return out
        
    def truncate(self, start: int, stop: int) -> "Location":
        """Truncate the location to fit within a specified range.

        Clips all positions to fit within the given start and stop coordinates.
        Positions that fall completely outside the range are removed. Positions
        that partially overlap are clipped to the boundaries.

        Args:
            start: The start coordinate of the truncation range (inclusive).
            stop: The stop coordinate of the truncation range (exclusive).

        Returns:
            A new :class:`Location` object with truncated positions.

        Raises:
            ValueError: If start is greater than stop.
        """
        if start > stop:
            raise ValueError("start > stop")
        
        out = self.copy()
        new_pos = []
        for s1, e1, sd1 in out.positions:
            if e1 <= start or s1 >= stop:
                continue
            s1 = max(s1, start)
            e1 = min(e1, stop)
            
            # Prevent 0-length intervals
            if s1 < e1:
                new_pos.append([s1, e1, sd1])
            
        out.positions = new_pos
        return out
        
    def shift(self, shift_len: int, seq_len: int = None) -> "Location":
        """Shift all positions by a specified offset.

        Adds the shift_len to all start and stop coordinates. If seq_len is
        provided, treats the coordinates as circular and wraps them around
        using modulo arithmetic.

        Args:
            shift_len: The number of bases to shift (can be negative).
            seq_len: Optional sequence length for circular genomes. If provided,
                coordinates are wrapped around using modulo arithmetic.

        Returns:
            A new :class:`Location` object with shifted positions.
        """
        out = self.copy()
        for i in range(len(out.positions)):
            out.positions[i][0] += shift_len
            out.positions[i][1] += shift_len
        
        if seq_len is not None:
            for i in range(len(out.positions)):
                out.positions[i][0] %= seq_len
                out.positions[i][1] %= seq_len
            
            new_loc = None
            for start, stop, strand in out.positions:
                if start > stop:
                    # Prevent creating a 0-length interval if stop wraps exactly to 0
                    if stop == 0:
                        loc = Location(start, seq_len, strand)
                    else:
                        loc = Location(start, seq_len, strand) + Location(0, stop, strand)
                else:
                    loc = Location(start, stop, strand)

                if new_loc is None:
                    new_loc = loc
                else:
                    new_loc += loc
            return new_loc
            
        return out
        
    def start(self) -> int:
        """Get the minimum start position across all positions.

        Returns:
            The smallest start coordinate among all positions.

        Raises:
            ValueError: If the location has no positions.
        """
        if not self.positions:
            raise ValueError("Location has no positions")
        return min(p[0] for p in self.positions)

    def end(self) -> int:
        """Get the maximum stop position across all positions.

        Returns:
            The largest stop coordinate among all positions.

        Raises:
            ValueError: If the location has no positions.
        """
        if not self.positions:
            raise ValueError("Location has no positions")
        return max(p[1] for p in self.positions)

    def strand(self) -> int:
        """
        Returns the strand if all positions have the same strand. Returns 0 if
        the strands are mixed or unknown.

        Returns:
            1 if all positions are on the forward strand, -1 if all positions
            are on the reverse strand, or 0 if strands are mixed or unknown.
        """
        if not self.positions:
            return 0
        strands = [entry[2] for entry in self.positions]
        if all(s == 1 for s in strands):
            return 1
        if all(s == -1 for s in strands):
            return -1
        return 0

    def extract(self, sequence: str) -> str:
        """
        Extracts the subsequence from the given sequence based on the location's
        coordinates.
        
        Args:
            sequence: The nucleotide sequence to extract from (as a string or
                any object that can be converted to a string).

        Returns:
            The extracted feature sequence as a string.
        """
        seq = str(sequence)
        parts = []
        
        for start, stop, strand in self.positions:
            segment = seq[int(start):int(stop)]
            if strand == -1:
                segment = str(Seq(segment).reverse_complement())
            parts.append(segment)
                
        return "".join(parts)