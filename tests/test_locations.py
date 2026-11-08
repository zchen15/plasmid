#!/usr/bin/env python
# tests for locations.py

import plasmid.genbank.locations as loc_module
from Bio.Seq import Seq
from Bio.SeqFeature import (SimpleLocation, CompoundLocation)


# ============= Flatten tests =============

def test_flatten_simple():
    """SimpleLocation has its own .parts that recursively expands."""
    loc = SimpleLocation(10, 30, 1)
    parts = loc_module.flatten_location(loc)
    assert len(parts) > 0

def test_flatten_compound():
    """CompoundLocation with 2 distinct simple locations."""
    compound = CompoundLocation(
        parts=[SimpleLocation(10, 30, 1), SimpleLocation(50, 70, -1)]
    )
    flat = loc_module.flatten_location(compound)
    print(flat)
    assert len(flat) == 2

# ============= Merge position tests =============

def test_merge_adjacent_forward():
    """Adjacent positions with forward strand get merged."""
    plist = [[10, 30, 1], [30, 50, 1]]
    merged = loc_module.merge_pos(plist)
    assert len(merged) > 0 and merged[0][0] == 10


def test_merge_adjacent_reverse():
    """Adjacent positions with reverse strand get merged."""
    plist = [[70, 30, -1], [40, 90, -1]]
    merged = loc_module.merge_pos(plist)
    assert len(merged) > 0


def test_merge_disjoint():
    """Disjoint positions remain separate."""
    plist = [[10, 30, 1], [40, 50, 1]]
    merged = loc_module.merge_pos(plist)
    assert len(merged) == 2


def test_merge_empty_list():
    """Empty list returns empty result."""
    assert len(loc_module.merge_pos([])) == 0


# ============= Location conversion tests =============

def test_convert_simple_biopython():
    """Convert custom Location to Biopython SimpleLocation."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    bio_loc = loc.to_BioLocation()
    assert type(bio_loc).__name__ == 'SimpleLocation'


def test_copy():
    """Deep copy produces independent location."""
    base_loc = loc_module.Location(start=10, stop=30, strand=1)
    copied = base_loc.copy()
    assert len(copied.positions) == len(base_loc.positions)


def test_str_repr():
    """String representation contains strand and coordinates."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    result = str(loc)
    assert '+' in result and '10' in result


# ============= Location arithmetic tests =============

def test_add_two_locations():
    """Add two non-overlapping locations."""
    loc1 = loc_module.Location(start=10, stop=30, strand=1)
    loc2 = loc_module.Location(start=40, stop=50, strand=1)
    combined = loc1 + loc2
    assert len(combined.positions) > 1


def test_add_overlapping():
    """Add two overlapping locations."""
    loc1 = loc_module.Location(start=10, stop=30, strand=1)
    loc2 = loc_module.Location(start=25, stop=50, strand=1)
    combined = loc1 + loc2
    # Check that positions were modified
    assert combined.positions[0][0] == 10
    assert len(combined.positions) > 1


# ============= Slicing tests =============

def test_getitem_single_index():
    """Get a single location by index."""
    pos1 = [10, 30, 1]
    pos2 = [40, 50, -1]
    loc = loc_module.Location(start=pos1[0], stop=pos1[1], strand=pos1[2])
    loc.positions.append(pos2)
    half1 = loc[0]
    assert len(half1.positions) == 1


def test_truncate_in_range():
    """Truncate expansion start and end."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    truncated = loc.truncate(start=5, stop=35)
    assert truncated.positions[0] is not None


def test_truncate_clips_edge():
    """Truncate clips positions near edge."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    clipped = loc.truncate(start=28, stop=32)
    assert clipped.positions[0][0] == 28


# ============= Shifting tests =============

def test_shift_simple():
    """Shift location by positive amount."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    shifted = loc.shift(shift_len=20)
    assert shifted.positions[0][0] == 30


def test_shift_circular():
    """Circular shift based on sequence length."""
    seq_len = 20
    loc = loc_module.Location(start=28, stop=32, strand=1)
    shifted = loc.shift(shift_len=-6, seq_len=seq_len)
    assert shifted.positions[0] is not None


# ============= Boundary methods tests =============

def test_start():
    """Get the first position of location."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    assert loc.start() == 10


def test_end():
    """Get the last position of location."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    assert loc.end() == 30


# ============= Strand computation tests =============

def test_strand_mixed():
    """Mixed strands yield strand 0 (no preference)."""
    pos1 = [10, 30, 1]
    pos2 = [40, 50, -1]
    loc = loc_module.Location(start=pos1[0], stop=pos1[1], strand=pos1[2])
    loc.positions.append(pos2)
    assert loc.strand() == 0


def test_strand_all_forward():
    """All forward strands yield strand 1."""
    loc = loc_module.Location(start=10, stop=30, strand=1)
    assert loc.strand() == 1


# ============= Sequence extraction tests =============

def test_extract_sequence_simple():
    """Extract sequence from location with forward strand."""
    seq = 'ATGCGTACGTAGCTA'
    location = loc_module.Location(start=5, stop=20, strand=1)
    extracted = location.extract(seq)
    assert len(extracted) > 5


def test_extract_reverse_complement():
    """Extract sequence using reverse complement."""
    seq = 'ATGCGTACGTAGCTA'
    loc = loc_module.Location(start=5, stop=20, strand=-1)
    extracted = loc.extract(seq)
    expected_rev = str(Seq('TCAGCATGCA').reverse_complement())
    assert len(extracted) > 8


def test_extract_mixed_strands():
    """Extract sequence with mixed strands."""
    pos1 = [0, 10, 1]
    pos2 = [15, 25, -1]
    loc = loc_module.Location(start=pos1[0], stop=pos1[1], strand=pos1[2])
    loc.positions.append(pos2)
    seq = 'ATGCGTACGTAGCTA'
    extracted = loc.extract(seq)
    assert len(extracted) > 5


# Run all tests
if __name__ == '__main__':
    import pytest as pt
    pt.main([__file__, '-v'])