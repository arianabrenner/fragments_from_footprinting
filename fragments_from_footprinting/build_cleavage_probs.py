"""
Includes functions that create and save the numpy array of cleavage probabilities by nucleotide position.
"""
from .params import *
import numpy as np

def canvas(with_attribution=True):
    """
    Placeholder function from cookiecutter use to ensure modules are properly acessed.
    Kept in this code base (at least temporarily) for testing

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote

def generate_cleav_prob(link_prob: float = 1.0, nuc_prob = 0.0, linker_length: int = link_len, wrap_bp: int = wrap) -> np.ndarray:

    """
    Generate an array where each item is the cleavage probability of the corresponding base pair.

    Parameters
    ----------
    link_prob : float
        default: 100% chance of cleavage at linker (1.0)
    nuc_prob: float or np.ndarray
        default: 0% change pf cleavage at nucleosome (0.0)
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
    """

    linker_cleavage_prob_arr = np.repeat(link_prob, linker_length)
    
    if type(nuc_prob) == float:
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
    return cleavage_prob

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
