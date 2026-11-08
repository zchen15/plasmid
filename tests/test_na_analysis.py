"""
test_rna_analyzer.py

Comprehensive pytest suite for the rna_analyzer library.
Tests cover the multiprocessing worker function, the main analyzer class, 
and the matrix conversion utilities without using pytest fixtures.
"""

import numpy as np
from plasmid.pcr.analysis import NucleicAcidAnalyzer, _analyze_single_sequence

# ==============================================================================
# Test Data
# ==============================================================================

SEQ_RNA_HAIRPIN = "GCGCGCGCGCGCGCGC"
SEQ_RNA_SIMPLE = "GGGAAACCC"
SEQ_RNA_UNSTRUCTURED = "A" * 20
SEQ_DNA = "ATCGATCGATCGATCG"


# ==============================================================================
# Tests for the Internal Worker Function
# ==============================================================================

def test_basic_rna_fold():
    """Tests basic RNA folding, ensuring MFE is negative and structure length matches."""
    res = _analyze_single_sequence(SEQ_RNA_SIMPLE, 37.0, False, 0.01, True)
    
    assert res['sequence'] == SEQ_RNA_SIMPLE
    assert res['length'] == len(SEQ_RNA_SIMPLE)
    assert res['dG_mfe'] < 0.0  # Structured RNA should have negative dG
    assert len(res['mfe_structure']) == len(SEQ_RNA_SIMPLE)
    assert res['dG_pf'] is not None
    assert len(res['pair_probabilities']) > 0


def test_dna_mode():
    """Tests that DNA mode runs without error and produces valid output."""
    res = _analyze_single_sequence(SEQ_DNA, 37.0, True, 0.01, True)
    
    assert res['sequence'] == SEQ_DNA
    assert res['dG_mfe'] < 0.0  # DNA should also fold


def test_compute_probs_false():
    """Tests that skipping probabilities results in None dG_pf and empty pairs."""
    res = _analyze_single_sequence(SEQ_RNA_SIMPLE, 37.0, False, 0.01, False)
    
    assert res['dG_pf'] is None
    assert res['pair_probabilities'] == []
    # MFE should still be calculated
    assert res['dG_mfe'] < 0.0


def test_prob_threshold_filtering():
    """Tests that a higher probability threshold results in fewer returned pairs."""
    # Get all pairs (threshold 0.0)
    res_all = _analyze_single_sequence(SEQ_RNA_SIMPLE, 37.0, False, 0.0, True)
    # Get only highly probable pairs (threshold 0.9)
    res_high = _analyze_single_sequence(SEQ_RNA_SIMPLE, 37.0, False, 0.9, True)
    
    assert len(res_all['pair_probabilities']) >= len(res_high['pair_probabilities'])
    
    # Verify all pairs in the high threshold list actually meet the threshold
    for i, j, prob in res_high['pair_probabilities']:
        assert prob >= 0.9


# ==============================================================================
# Tests for the Main Analyzer Class
# ==============================================================================

def test_initialization_defaults():
    """Tests that attributes are set correctly upon initialization."""
    analyzer = NucleicAcidAnalyzer(
        num_workers=2, 
        temperature=37.0, 
        is_dna=False, 
        prob_threshold=0.01, 
        compute_probs=True
    )
    
    assert analyzer.num_workers == 2
    assert analyzer.temperature == 37.0
    assert analyzer.is_dna is False
    assert analyzer.prob_threshold == 0.01
    assert analyzer.compute_probs is True


def test_analyze_batch_empty():
    """Tests that an empty input list returns an empty output list."""
    analyzer = NucleicAcidAnalyzer(num_workers=2)
    res = analyzer.analyze_batch([])
    assert res == []


def test_analyze_batch_single():
    """Tests processing a single sequence."""
    analyzer = NucleicAcidAnalyzer(num_workers=2)
    res = analyzer.analyze_batch([SEQ_RNA_SIMPLE])
    
    assert len(res) == 1
    assert res[0]['sequence'] == SEQ_RNA_SIMPLE
    assert res[0]['dG_mfe'] < 0.0


def test_analyze_batch_multiple_order_preserved():
    """Tests that multiprocessing preserves the order of the input sequences."""
    analyzer = NucleicAcidAnalyzer(num_workers=2)
    seqs = [SEQ_RNA_HAIRPIN, SEQ_RNA_SIMPLE, SEQ_RNA_UNSTRUCTURED]
    res = analyzer.analyze_batch(seqs)
    
    assert len(res) == 3
    for i, seq in enumerate(seqs):
        assert res[i]['sequence'] == seq


def test_analyze_batch_fast_mode():
    """Tests the fast mode (compute_probs=False) via the class interface."""
    analyzer = NucleicAcidAnalyzer(num_workers=2, compute_probs=False)
    res = analyzer.analyze_batch([SEQ_RNA_SIMPLE])
    
    assert len(res) == 1
    assert res[0]['dG_pf'] is None
    assert res[0]['pair_probabilities'] == []
    assert res[0]['dG_mfe'] < 0.0


# ==============================================================================
# Tests for Matrix Conversion
# ==============================================================================

def test_matrix_shape_and_symmetry():
    """Tests that the resulting matrix has the correct shape and is symmetrical."""
    analyzer = NucleicAcidAnalyzer(num_workers=2)
    
    mock_res = {
        'length': 5,
        'pair_probabilities': [(0, 4, 0.9), (1, 3, 0.8)]
    }
    matrix