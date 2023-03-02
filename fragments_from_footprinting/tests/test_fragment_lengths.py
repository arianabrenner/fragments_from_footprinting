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