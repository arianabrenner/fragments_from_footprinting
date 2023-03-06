import sys
import pytest
# from .params import *
import numpy as np
import random
import fragments_from_footprinting as ff

def test_get_breaks_to_try():
    """
    Test that introducing the given number of attempted breaks yields the intended number of actual breaks

    TO DO:
    -----
    Load in the actual used parameters to ff.generate_cleav_prob() to test if the proper number of 
    breaks are attempted for the situation actaully being attempted.
    """
    example_cp = ff.generate_cleav_prob()
    breaks_to_try, expected_breaks = ff.get_breaks_to_try(example_cp)
    trials = 1000
    num_breaks_all_trials = []
    nts = len(example_cp)
    for i in range(0,trials):
        # random.seed(10)
        locations_to_attempt_cut = [random.randint(0, nts-1) for i in range(0,breaks_to_try)]
        probs_at_sites_2_cut = example_cp[locations_to_attempt_cut]
        # random.seed(10)
        cuts = [random.random() < i for i in probs_at_sites_2_cut]
        # Number of successful breaks for a given trial
        tot_breaks = np.sum(cuts)
        num_breaks_all_trials.append(tot_breaks)    
    # Compute average number of oberved breaks across all trials
    avg_num_breaks = np.mean(num_breaks_all_trials)
    diff_amount = np.round(np.absolute(expected_breaks-avg_num_breaks),2)
    testresult = diff_amount < 1
    print("Average Num. Breaks from Simulation and Expected Num. Breaks are within 1bp of each other: "+str(testresult))
    print("Observed Difference is : "+str(diff_amount)+'bp')
    assert diff_amount < 1

def test_get_fld():
    example_cp = ff.generate_cleav_prob()
    frags, mids = ff.get_fld(example_cp, trials = 10, save_data = 0)
    assert len(frags) == len(mids)
    # Midpoints must be greater than 0 and less than total number of nucs (to do: import this number)
    assert np.min(mids) > 0
    # Add assertion that the max frag len is less than total numb nucleotides

def test_vplot_data():
    """
    Test that the vplot_data returns an array of the correct size
    """
    example_cp = ff.generate_cleav_prob()
    frags, mids = ff.get_fld(example_cp, trials = 10, save_data = 0)
    fm_df = ff.frag_mid_df(frags, mids)
    max_frag_ = 300
    dist_from_center_ = 500
    bin_lens_ = 3 
    bin_locs_ = 20
    vplot_arr = ff.vplot_data(fm_df, max_frag = max_frag_, dist_from_center = dist_from_center_, bin_lens = bin_lens_, bin_locs = bin_locs_, save_data = 0)
    expected_rows = (2. * dist_from_center_) / bin_locs_
    expected_cols = max_frag_ / bin_lens_
    observed_shape = np.shape(vplot_arr)
    assert (expected_rows == observed_shape[0]) & (expected_cols == observed_shape[1])
