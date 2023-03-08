"""A Python package for predicting fragment length distributions and v-plots resulting from pa protein's ability to shield DNA from radiation induced damage. A particular femphasis is placed on protection specifically by nucleosomes."""

# Add imports here
import os
from .build_cleavage_probs import *
from .fragment_lengths import *
from .params import *
from .plot import *

from ._version import __version__

# Make a directory for the outputs if one does not exist
for DIR in ["intermed_data", "plots"]:
    CHECK_FOLDER = os.path.isdir(DIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(DIR)
        print("created folder : ", DIR)

    else:
        print(DIR, "folder already exists.")