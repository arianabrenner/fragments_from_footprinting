"""
Includes functions that create and save the numpy array of cleavage probabilities by nucleotide position.
"""
from .params import *
from numpy import np

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

def temp(with_attribution=True):
    print(nrl)
    return nrl


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
