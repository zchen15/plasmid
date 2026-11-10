import RNA
import numpy as np
import multiprocessing
from typing import List, Dict, Any, Optional, Tuple
import os


# ==============================================================================
# Multiprocessing Worker Function
# ==============================================================================

def analyze_sequences(
    sequences: List[str],
    temperature: float,
    is_dna: bool,
    prob_threshold: float,
    compute_probs: bool
) -> List[Dict[str, Any]]:
    """
    Internal top-level worker function to calculate thermodynamic properties
    for a list of nucleic acid sequences.

    Sequences (for example, a batch of up to 100) are processed in order and
    the results are returned as a list of result dictionaries, which reduces
    the number of multiprocessing tasks and the IPC overhead. This function
    MUST remain at the top level of the module to be picklable by the
    multiprocessing module.

    Args:
        sequences (List[str]): The nucleic acid sequences (A, U, C, G for
                        RNA; A, T, C, G for DNA).
        temperature (float): Temperature in degrees Celsius for the 
                             thermodynamic calculations.
        is_dna (bool): If True, uses DNA thermodynamic parameters instead 
                       of RNA.
        prob_threshold (float): The minimum base pair probability required 
                                to include a pair in the output list. 
                                Filters out noise and saves IPC bandwidth.
        compute_probs (bool): If True, computes the Partition Function and 
                              extracts base pair probabilities. If False, 
                              skips this step for significantly faster execution.

    Returns:
        List[Dict[str, Any]]: A list of result dictionaries, one per input
        sequence, preserving the input order. Each dictionary contains:
            - 'sequence' (str): The input sequence.
            - 'length' (int): Length of the sequence.
            - 'dG_mfe' (float): Minimum Free Energy in kcal/mol.
            - 'mfe_structure' (str): Dot-bracket notation of the MFE structure.
            - 'dG_pf' (float or None): Partition Function free energy in kcal/mol 
                                       (None if compute_probs is False).
            - 'pair_probabilities' (List[Tuple[int, int, float]]): A sparse list 
              of base pairs (i, j, prob) using 0-based indexing.
    """
    # Initialize ViennaRNA model details
    md = RNA.md()
    md.temperature = temperature
    if is_dna:
        md.DNA = 1  # Switch to DNA parameters if requested

    results = []
    for sequence in sequences:
        # Create a fold compound (isolated C-struct, safe for multiprocessing)
        fc = RNA.fold_compound(sequence, md)
        
        # 1. Calculate Minimum Free Energy (dG) - Always computed
        mfe_structure, dG_mfe = fc.mfe()
        
        # 2. Conditionally Calculate Partition Function and Pair Probabilities
        dG_pf = None
        pair_probs = []
        
        if compute_probs:
            dG_pf = fc.pf() 
            
            # ViennaRNA uses 1-based indexing; we convert to 0-based for Python/NumPy
            for i in range(1, len(sequence) + 1):
                for j in range(i + 1, len(sequence) + 1):
                    prob = fc.pr(i, j)
                    if prob >= prob_threshold:
                        pair_probs.append((i - 1, j - 1, prob))
                    
        results.append({
            "sequence": sequence,
            "dG_mfe": dG_mfe,
            "mfe_structure": mfe_structure,
            "dG_pf": dG_pf,
            "pair_probabilities": pair_probs
        })

    return results

def pp_to_ppmatrix(
    analysis_result: Dict[str, Any],
) -> np.ndarray:
    """
    Converts the sparse list of base pair probabilities from a single 
    analysis result into a dense 2D NumPy probability matrix.
    
    This is highly useful for downstream machine learning tasks, structural 
    visualization, or generating heatmaps using libraries like matplotlib 
    or seaborn.
    
    Args:
        analysis_result (Dict[str, Any]): A single result dictionary 
            returned by `NucleicAcidAnalyzer.analyze_batch`.
    
    Returns:
        np.ndarray: A 2D NumPy array of shape (N, N) where N is the 
        sequence length. The matrix is symmetrical, with values representing 
        the base pair probabilities. Unpaired positions remain 0.0.
    """
    seq_len = len(analysis_result["sequence"])
    matrix = np.zeros((seq_len, seq_len), dtype=np.float32)
    
    for i, j, prob in analysis_result["pair_probabilities"]:
        matrix[i, j] = prob
        matrix[j, i] = prob  # Ensure the matrix is symmetrical
    
    return matrix


# ==============================================================================
# Main Library Class
# ==============================================================================

class NucleicAcidAnalyzer:
    """
    A multiprocessing-enabled library for high-throughput nucleic acid 
    thermodynamic analysis using the ViennaRNA package.
    
    This class manages a pool of worker processes to evaluate batches of 
    sequences in parallel. It is designed to handle the heavy computational 
    load of RNA/DNA folding while minimizing Inter-Process Communication (IPC) 
    overhead by returning sparse data structures.
    
    Attributes:
        num_workers (int): Number of CPU cores utilized by the process pool.
        temperature (float): Default temperature (°C) for calculations.
        is_dna (bool): Default flag indicating whether to use DNA parameters.
        prob_threshold (float): Default minimum probability threshold for pairs.
        compute_probs (bool): Default flag indicating whether to compute the 
                              Partition Function.
    """
    
    def __init__(
        self, 
        temperature: float = 37.0, 
        is_dna: bool = False,
        prob_threshold: float = 0.01,
        compute_probs: bool = True,
        batch_size:int = 100,
    ):
        """
        Initializes the NucleicAcidAnalyzer with the specified computational 
        and thermodynamic parameters.

        Args:
            num_workers (int, optional): Number of CPU cores to use for the 
                multiprocessing pool. Defaults to the total number of logical 
                cores available on the machine (os.cpu_count()).
            temperature (float, optional): Temperature in Celsius for 
                thermodynamic calculations. Defaults to 37.0.
            is_dna (bool, optional): If True, uses DNA thermodynamic parameters 
                instead of RNA. Defaults to False.
            prob_threshold (float, optional): Minimum base pair probability 
                (0.0 to 1.0) to include in the output. Lowering this increases 
                memory/IPC overhead. Defaults to 0.01.
            compute_probs (bool, optional): If True, calculates the Partition 
                Function and base pair probabilities. If False, only calculates 
                the MFE (dG), resulting in a massive speedup. Defaults to True.
        """
        self.num_workers = os.cpu_count()
        self.temperature = temperature
        self.is_dna = is_dna
        self.prob_threshold = prob_threshold
        self.compute_probs = compute_probs
        self.batch_size = batch_size
        
    def analyze_batch(self, sequences: List[str]) -> List[Dict[str, Any]]:
        """
        Analyzes a batch of nucleic acid sequences in parallel using a 
        multiprocessing Pool.

        This method splits the sequences into batches of `batch_size` and
        dispatches one worker task per batch to the process pool. The per-batch
        result lists are then concatenated into a single list. Because it uses
        the 'spawn' start method, it is safe to run on Windows, macOS, and Linux
        without risking C-extension deadlocks.

        Args:
            sequences (List[str]): A list of nucleic acid sequences to analyze. 
                Sequences should only contain valid characters (A, U, C, G for 
                RNA; A, T, C, G for DNA).

        Returns:
            List[Dict[str, Any]]: A list of dictionaries containing the analysis 
            results for each sequence, preserving the order of the input list. 
            See `analyze_sequences` for the dictionary structure.
        """

        if not sequences:
            return []

        # Prepare one task per batch of up to `batch_size` sequences
        args = [
            (sequences[start:start+self.batch_size], self.temperature, self.is_dna,
             self.prob_threshold, self.compute_probs)
            for start in range(0, len(sequences), self.batch_size)
        ]

        # Create the Pool and execute the batches
        ctx = multiprocessing.get_context('spawn')
        with ctx.Pool(processes=self.num_workers) as pool:
            # starmap unpacks the tuple arguments into the worker function
            batches = pool.starmap(analyze_sequences, args)

        # Concatenate the per-batch result lists, preserving the input order
        return [result for batch in batches for result in batch]
