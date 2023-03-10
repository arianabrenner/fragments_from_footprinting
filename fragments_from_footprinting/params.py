"""
Provides parameters for the simulation. Parameters are saved in params.csv file that is passed as input
"""
import pandas as pd
import numpy as np 

def set_parameters(params_file: str = "params.csv"):
    """
    Reads parameter values in from params.csv. Must run package in directory with params.csv for this to work.

    Parameters
    ----------
    params_file : str
        Filepath to params.csv

    Returns
    -------
    nrl : int
        Nucleosome repeat length in nucleotides.
    
    wrap : int
        Number of nucleotides considered to be wrapped around the nucleosome.
    
    num_nucs : int
        Number of nucleosomes to simulate.

    max_fragment_length : int
        Largest fragment length tp consider (in nucleotides).

    distance_from_frag_center : int 
        Number of nucleotides away from the simulated fragment's center to consider. 
    """
    params_df = pd.read_csv(params_file)
    
    # Fiber parameters
    nrl = int(np.array(params_df.nrl))
    wrap = int(np.array(params_df.wrap))
    num_nucs = int(np.array(params_df.num_nucs))

    # Simulation parameters
    break_rate = int(np.array(params_df.break_rate))
    num_trials = int(np.array(params_df.num_trials))
    dyad_bool = int(np.array(params_df.dyad_bool))
    # dyad_width = int(np.array(params_df.dyad_width))
    
    # Variables for plotting
    max_fragment_length = int(np.array(params_df.max_fragment_length))
    distance_from_frag_center = int(np.array(params_df.distance_from_frag_center))

    return nrl, wrap, num_nucs, max_fragment_length, distance_from_frag_center, break_rate, num_trials, dyad_bool

nrl, wrap, num_nucs, max_fragment_length, distance_from_frag_center, break_rate, num_trials, dyad_bool = set_parameters()


dyad_width = 15
# Secondary Parameters
link_len = nrl - wrap
fiber_length = nrl * num_nucs + link_len
fiber_midpoint = fiber_length / 2
