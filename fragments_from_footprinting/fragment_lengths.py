"""
Run simulations to determine resulting fragment length distributions
"""
from .params import *
import numpy as np
import random
import pandas as pd
import seaborn as sns

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

def get_frag_lens(pot_cut_locs, cut_bool):

    """
    Get fragment lengths and locations for a single trial

    Parameters
    ----------
    pot_cut_locs : np.ndarray
        Array of nucleotide positions on which to attempt breaks.
    
    cut_bool : np.ndarray
       Array of booleans indicating whether attempted break was successful

    Returns
    -------
    fragments : np.ndarray
        Array of fragment lengths from one trial.

    midpoints : np.ndarray
        The center location of these fragments.

    TO DO
    -----
    Can create an option that uses the location and the lag location to make a contact map as opposed to a v-plot

    """
    # Convert input into dataframe
    cut_opportunity_df_w_dups = pd.DataFrame({'loc': pot_cut_locs, 'cuts': cut_bool}).sort_values(by='loc')
    # Account for duplicates caused by two attempted cuts at the same location
    cut_opportunity_df = cut_opportunity_df_w_dups.groupby(['loc']).max()
    cut_opportunity_df = cut_opportunity_df.reset_index()
    # subset to successful breaks
    cuts_df = cut_opportunity_df[cut_opportunity_df['cuts']==1]
    
    # Make new column showing next sucessful break location 
    lag_series = np.append(np.array(cuts_df['loc'])[1:], np.nan)
    cuts_df_w_lag = cuts_df.copy()
    cuts_df_w_lag['lag_loc'] = lag_series
    
    # Make columns with fragment length - this is the diff between the location and the lag location - and location (midpoint)
    cuts_df_w_lag['frag_lens'] = cuts_df_w_lag['lag_loc'] - cuts_df_w_lag['loc']
    cuts_df_w_lag['mid_idx'] = np.round((cuts_df_w_lag['lag_loc'] + cuts_df_w_lag['loc'])/2.0,0)
    fragments = np.array(cuts_df_w_lag.frag_lens)[:-1] #index to remove nan that is at the end
    midpoints = np.array(cuts_df_w_lag.mid_idx)[:-1] #index to remove nan that is at the end
    return fragments, midpoints


def get_fld(cleavage_prob: np.ndarray, trials: int = 10000, break_rate: int = 150, xmin: int = 50): 
    """
    Generates fragment length distribution

    Parameters
    ----------
    cleavage_prob : np.ndarray
        Probability of cleavage corresponding to each nucleotide position.
    
    trials : int
        Number of trials (i.e. number of times the given number of breaks are attempted on the simulated nucleotide array).

    break_rate : int
        1 break per this many nucleotides.
    
    xmin : int
        Minimum fragment length to consider. Default is 50nt to compare to RICC-seq simulated data.

    Returns
    -------
    frag_lens_all_trials : np.ndarray
        Array of fragment lengths from all trials (the fragment length distribution).

    midpts_all_trials : np.ndarray
        The center location of these fragments relative to the simulated nucleotide array. 

    TO DO
    -----
    Maybe get rid of minimum fragment length requirement.

    """
    # Breaks per nucleotide
    bpnt = 1./break_rate
    nts = len(cleavage_prob)
    breaks_to_try, exp_breaks = get_breaks_to_try(cleavage_prob, breaks_per_nt = bpnt)
    # List to store fragment lengths
    frag_lens_all_trials = []
    # Lsit to store midpoints of the fragment
    midpts_all_trials = []
    for i in range(0,trials):
        # random.seed(10)
        locations_to_attempt_cut = [random.randint(0, nts-1) for i in range(0,breaks_to_try)]
        probs_at_sites_2_cut = cleavage_prob[locations_to_attempt_cut]
        # random.seed(10)
        # Draw from a uniform random distribution [0,1). If that number is < the probability of a cut, the cut is successful.  
        cuts = [random.random() < i for i in probs_at_sites_2_cut]
        frags, midpts = get_frag_lens(pot_cut_locs = locations_to_attempt_cut, cut_bool = cuts)
        frag_lens_all_trials.append(frags)
        midpts_all_trials.append(midpts)
    # Flatten list of arrays
    frag_lens_all_trials = np.concatenate(frag_lens_all_trials).ravel() 
    midpts_all_trials = np.concatenate(midpts_all_trials).ravel() 
    idxs = np.where(frag_lens_all_trials>xmin)
    frag_lens_all_trials = frag_lens_all_trials[idxs]
    # subset to same indices as frag_lens_all_trials
    midpts_all_trials = midpts_all_trials[idxs]
    np.save('frag_lens.npy', frag_lens_all_trials)
    np.save('frag_midpts.npy', midpts_all_trials)
    return frag_lens_all_trials, midpts_all_trials

def bin_frags(frag_lens: np.ndarray, midpts: np.ndarray, bin_width: float = 10.):
    """
    Make fragment lengths and locations into a pandas dataframe. 
    Limit to only fragment lengths of interest and add optional binning for sparse data.

    Parameters
    ----------
    frag_lens : np.ndarray
        Array of fragment lengths from all trials (the fragment length distribution).

    midpts : np.ndarray
        The center location of these fragments relative to the simulated nucleotide array. 

    bin_width : float
        Number of nucleotides to aggregate the midpoints data over for easy visualization purposes and faster runtime.        

    Returns
    -------
    frags_and_mids : pd.DataFrame
        Dataframe with fragment lengths, locations, and binned coordinates
    """

    # Round to nearest 10 to make the bin mins whole numbers when bin width is set to 10
    min_range = round(np.min(midpts), -1)
    max_range = round(np.max(midpts), -1)
    bin_boundaries = list(np.linspace(min_range,max_range, 1+int((max_range-min_range)/bin_width)))
    #label bin based on bin min value
    bin_labels = bin_boundaries[:-1] 
    frags_and_mids = pd.DataFrame({'frag_len': frag_lens,
                      'midpoints': midpts})
    frags_and_mids['bin_mins'] = pd.cut(frags_and_mids['midpoints'], bins=bin_boundaries, labels=bin_labels)
    
    #Save data
    frags_and_mids.to_csv('fragment_lens_and_locations.csv')
    return frags_and_mids

"""
move to differnet module
"""
# def plot_fld(frags_and_mids, width = 2):
#     max_y_bin = 1000
#     min_y = 0#50
#     data_for_hist = frags_and_mids[frags_and_mids.frag_len<max_y_bin].frag_len
#     stored_histogram = sns.histplot(data=data_for_hist, binwidth=width, stat= 'probability')
#     plt.show()
#     return




