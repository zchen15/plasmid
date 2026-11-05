#!/usr/bin/env python
# tests for locations.py

import plasmid.genbank.locations as loc_module
from Bio.Seq import Seq
from Bio.SeqFeature import (SimpleLocation, CompoundLocation)


# --- Helper function tests ---
def test_flatten_simple():
    loc = SimpleLocation(10, 30, 1)
    parts = loc_module.flatten_location(loc)
    assert len(parts) == 1
    assert type(parts[0]).__name__ == 'SimpleLocation'


def test_flatten_compound():
    parts_bpo = [SimpleLocation(10, 30, 1), SimpleLocation(50, 70, -1)]
    loc = CompoundLocation(parts_bpo)
    flat = loc_module.flatten_location(loc)
    assert len(flat) == 2


def test_flatten_nested():
    inner = [SimpleLocation(10, 30, 1), SimpleLocation(50, 70, -1)]
    outer = CompoundLocation([CompoundLocation(inner)])
    flat = loc_module.flatten_location(outer)
    assert len(flat) == 2


def test_merge_adjacent_forward():
    # Merging adjacent forward strand positions (e1 == s2)
    plist = [[10, 30, 1], [30, 50, 1]]
    merged = loc_module.merge_pos(plist)
    assert len(merged) == 1
    assert merged[0] == [10, 50, 1]


def test_merge_adjacent_reverse():
    plist = [[50, 30, -1], [70, 40, -1]]
    merged = loc_module.merge_pos(plist)
    assert len(merged) == 1
    assert merged[0] == [50, 70, -1]


def test_merge_disjoint():
    plist = [[10, 30, 1], [40, 50, 1]]
    merged = loc_module.merge_pos(plist)
    assert len(merged) == 2


def test_merge_empty():
    assert len(loc_module.merge_pos([])) == 0


# --- Location class tests ---
class TestLocation:

    def test_convert_simple_biopython(self):
        loc = loc_module.Location(10, 30, 1)
        bio_loc = loc.to_BioLocation()
        assert type(bio_loc).__name__ == 'SimpleLocation'

    def test_copy(self):
        loc = loc_module.Location(10, 30, 1)
        copied = loc.copy()
        assert len(copied.positions) == len(loc.positions)
        copied.positions.append([40, 50, 1])
        assert len(loc.positions) == 1

    def test_str_repr(self):
        loc = loc_module.Location(10, 30, 1)
        assert str(loc) == '[10:30](+)' or str(loc) == '[10:30](0)'

    def test_add_two_locations(self):
        loc1 = loc_module.Location(10, 30, 1)
        loc2 = loc_module.Location(40, 50, 1)
        combined = loc1 + loc2
        assert len(combined.positions) == 1

    def test_add_overlapping(self):
        loc1 = loc_module.Location(10, 30, 1)
        loc2 = loc_module.Location(25, 50, 1)
        combined = loc1 + loc2
        assert combined.positions[0][0] == 10
        assert combined.positions[0][1] == 50

    def test_getitem_single_index(self):
        loc = loc_module.Location([[10, 30, 1], [40, 50, -1]])
        half1 = loc[0]
        assert len(half1.positions) == 1

    def test_getitem_slice(self):
        loc = loc_module.Location([[10, 30, 1], [40, 50, -1]])
        first_half = loc[1]
        assert len(first_half.positions) == 2

    def test_truncate_in_range(self):
        loc = loc_module.Location(10, 30, 1)
        truncated = loc.truncate(5, 35)
        assert truncated.positions[0][0] == 5

    def test_shift_simple(self):
        loc = loc_module.Location(10, 30, 1)
        shifted = loc.shift(20)
        assert shifted.positions[0][0] == 30
        assert shifted.positions[0][1] == 50

    def test_shift_circular(self):
        loc = loc_module.Location(28, 32, 1)
        seq_len = 32
        shifted = loc.shift(-6, seq_len)

    def test_start(self):
        loc = loc_module.Location(10, 30, 1)
        assert loc.start() == 10

    def test_end(self):
        loc = loc_module.Location(10, 30, 1)
        assert loc.end() == 30

    def test_strand_mixed(self):
        loc = loc_module.Location([[10, 30, 1], [40, 50, -1]])
        assert loc.strand() == 0

    def test_strand_all_forward(self):
        loc = loc_module.Location(10, 30, 1)
        assert loc.strand() == 1

    def test_extract_sequence_simple(self):
        seq = 'ATGCGTACGTAGCTA'
        location = loc_module.Location(5, 20, 1)
        extracted = location.extract(seq)
        assert len(extracted) == 15 or len(extracted) == 16

    def test_extract_reverse_complement(self):
        seq = 'ATGCGTACGTAGCTA'
        loc = loc_module.Location(5, 20, -1)
        extracted = loc.extract(seq)
        expected_rev = str(Seq('TCAGCATGCA').reverse_complement())
        assert extracted == expected_rev

    def test_extract_mixed_strands(self):
        seq = 'ATGCGTACGTAGCTA'
        loc = loc_module.Location([[0, 10, 1], [15, 25, -1]])
        extracted = loc.extract(seq)
        assert len(extracted) == 35 or len(extracted) == 36


if __name__ == '__main__':
    import pytest as pt
    pt.main([__file__, '-v'])