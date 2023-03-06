"""
Provides parameters for the simulation
"""

nrl = 212
wrap = 147
link_len = 212 - wrap
num_nucs = 121
# Variables for plotting
max_fragment_length = 1000
distance_from_frag_center = 1200

# Secondary Parameters to add 
fiber_length = nrl * num_nucs + link_len
fiber_midpoint = fiber_length / 2