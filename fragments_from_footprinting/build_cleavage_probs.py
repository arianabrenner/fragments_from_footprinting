"""
Includes functions that create and save the numpy array of cleavage probabilities by nucleotide position.
"""
from .params import *
import numpy as np


def make_dyad_array(dyad_prob: float = 1.0, nuc_prob: float = 0.0, wrap_bp: int = wrap, dyad_width: int = dyad_width):
    """
    Generate cleavage probability array for nucleosome

    Parameters
    ----------
    dyad_prob : float
        Maximim probability of dyad cleavage
    nuc_prob : float or np.ndarray
        default: 0% chance of cleavage at nucleosome (0.0)
    dyad_width : int 
        number of nucleotides in the dyad 
    wrap_bp : int
        default: set in params.py
        Number of base pairs wrapped around the nucleosome.
    """
    nuc_prob_arr = np.repeat(nuc_prob, wrap)
    # print("Dyad Probability: "+str(dyad_prob), flush=True)
    # np.linspace(... dyad_prob-nuc_prob,...) subtract nuc prob because you later add the dyad_diff to nuc_prob_arr
    gradient = np.linspace(0, dyad_prob-nuc_prob, int(1+np.round(dyad_width/2)))[1:]
    num_zeros_to_add = (wrap-2*len(gradient)+1)/2
    zeros_to_add = np.zeros(int(num_zeros_to_add))
    dyad_diff = np.concatenate((zeros_to_add,gradient,np.flip(gradient)[1:],zeros_to_add))
    nuc_prob_arr_dyad = nuc_prob_arr + dyad_diff
    return nuc_prob_arr_dyad

def generate_cleav_prob(link_prob: float = 1.0, nuc_prob: float = 0.0, linker_length: int = link_len, wrap_bp: int = wrap, dyad_bool = dyad_bool, dyad_width: int = dyad_width) -> np.ndarray:

    """
    Generate an array where each item is the cleavage probability of the corresponding base pair.

    Parameters
    ----------
    link_prob : float
        default: 100% chance of cleavage at linker (1.0)
    nuc_prob: float or np.ndarray
        default: 0% chance of cleavage at nucleosome (0.0)
    linker_length: int
        default: set in params.py
        Number of base pairs in the DNA linking neighboring nucleosomes. 
    wrap_bp : int
        default: set in params.py
        Number of base pairs wrapped around the nucleosome.
    Returns
    -------
    cleavage_prob : np.ndarray
        Probability of cleavage corresponding to each nucleotide position. 

    TO DO
    -----
    Make the nuc_prob and link_prob variables in the params.csv that are drawn in.
    """

    linker_cleavage_prob_arr = np.repeat(link_prob, linker_length)

    # Make dyad array if dyad_bool = True
    if dyad_bool == 1:
        # assume max dyad prob is equal to the linker prob unless otherwise specified
        nuc_prob_arr = make_dyad_array(dyad_prob = link_prob, nuc_prob = nuc_prob, wrap_bp = wrap_bp, dyad_width = dyad_width)
    
    elif type(nuc_prob) == float:
        # If nuc_prob is a number, repeat that number n times where n=num nucleotides wrapped around nucleosome
        nuc_prob_arr =np.repeat(nuc_prob, wrap_bp)

    elif isinstance(nuc_prob, (np.ndarray)):
        # If nuc_prob is an array, use it as is
        nuc_prob_arr = nuc_prob

    else:
        raise ValueError("nuc_prob_arr not properly defined")

    
    # Pass an error if the probability array is the wrong size        
    if len(nuc_prob_arr) != wrap:    
        raise ValueError('Error: nuc_prob_arr not correct length')

    
    #Iterate through num_nucs times to build complete cleavage array
    cleavage_prob = linker_cleavage_prob_arr
    for i in range(0,num_nucs):
        cleavage_prob = np.concatenate((cleavage_prob,nuc_prob_arr))
        cleavage_prob = np.concatenate((cleavage_prob,linker_cleavage_prob_arr))
    #convert from list to array
    cleavage_prob = np.array(cleavage_prob)
    np.save('intermed_data/cleavage_prob.npy', cleavage_prob)
    return cleavage_prob

if __name__ == "__main__":
    print("build_cleavage.py invoked")
    # Do something if this file is invoked on its own