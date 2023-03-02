"""
Run simulations to determine resulting fragment length distributions
"""
from .params import *
import numpy as np

def get_breaks_to_try(cleavage_array: np.ndarray, breaks_per_nt: float = 1./150.):
    '''Input: the breakage rate as breaks per base pair 
           Output: Number of breaks to attempt on a structure to observe the given breakage rate on average  
    '''
    
    """
    Computes the number of breaks to attempt on a structure to observe the given breakage rate, on average

    Parameters
    ----------
    cleavage_array : np.ndarray
        Probability of cleavage corresponding to each nucleotide position.
    
    breaks_per_nt: float
        Intended number of breaks per nucleotide.

    Returns
    -------
    breaks_to_try : int
         Number of breaks to attempt on the simulated nucleotide strand to observe the given breakage rate on average.
    expected_breaks: float
        Expected value when `breaks_to_try` is attempted on the `cleavage_array`
    """
    nts = len(cleavage_array) 
    aggregate_prob = np.sum(cleavage_array)/nts
    expected_breaks = breaks_per_nt * nts
    breaks_to_try = int(expected_breaks / aggregate_prob)
    return breaks_to_try, expected_breaks

def get_fld(cleavage_prob, trials= 100, break_rate = 150.): 
    '''Input: Array of cleavage probabilites and other parameters
              'break_rate' = 1 break per this many nucleotides
       Output: Array of fragment lengths from all trials (the fragment length distribution)
    '''
    bpbp = 1./break_rate
    bps = len(cleavage_prob)
    breaks_to_try, exp_breaks = breaksToTry(cleavage_prob, breaks_per_bp = bpbp)
    # List to store fragment lengths
    frag_lens_all_trials = []
    # Lsit to store midpoints of the fragment
    midpts_all_trials = []
    for i in range(0,trials):
        # random.seed(10)
        locations_to_attempt_cut = [random.randint(0, bps-1) for i in range(0,breaks_to_try)]
        probs_at_sites_2_cut = cleavage_prob[locations_to_attempt_cut]
        # random.seed(10)
        cuts = [random.random() < i for i in probs_at_sites_2_cut]
        frags, midpts = getFragLens(pot_cut_locs = locations_to_attempt_cut, cut_bool = cuts)
        frag_lens_all_trials.append(frags)
        midpts_all_trials.append(midpts)
    #Flatten list of arrays
    frag_lens_all_trials = np.concatenate(frag_lens_all_trials).ravel() 
    midpts_all_trials = np.concatenate(midpts_all_trials).ravel() 
    xmin = 50
    idxs = np.where(frag_lens_all_trials>xmin)
    frag_lens_all_trials = frag_lens_all_trials[idxs]
#     subset to same indices as frag_lens_all_trials
    midpts_all_trials = midpts_all_trials[idxs]
    return frag_lens_all_trials, midpts_all_trials