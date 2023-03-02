"""
Unit and regression test for the fragments_from_footprinting package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
# from .params import *
import numpy as np
import fragments_from_footprinting as ff


def test_generate_cleav_prob():
    """
    Test that generate_cleav_prob functions generates a cleavage probability
    array for both a float (constant cleavage probability) or an array (variable cleavage probability)

    TO DO:
    -figure out how to use the wrap values from .params instead of manually entering.
    """
    wrap = 147
    prob_float = 0.5
    prob_array = np.repeat(prob_float, wrap)

    # Test that, for constant probability, passing one float value yields
    # the same result as passing an array of length `wrap` where each element is that same float value
    assert np.array_equal(ff.generate_cleav_prob(nuc_prob = prob_float), ff.generate_cleav_prob(nuc_prob = prob_array))

    # Verify that error is thrown if array is the wrong length using an array that is shorter than `wrap`
    with pytest.raises(ValueError):
        ff.generate_cleav_prob(nuc_prob = prob_array[:-1])