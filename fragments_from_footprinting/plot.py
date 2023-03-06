
"""
Includes functions for plotting the fragment length distribution, vplot, and cleavage probability
"""

# Import Modules 
from .fragment_lengths import *
from .params import max_fragment_length, distance_from_frag_center, fiber_midpoint
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
from sklearn import preprocessing
import seaborn as sns 

""" Set visual standards """

plt.rcdefaults()
# Could use Style Guide Instead of Custom
thickness = 2
fsize = 18
mpl.rcParams['lines.linewidth'] = thickness
mpl.rcParams['lines.linestyle'] = '-'
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.titlesize'] = fsize
mpl.rcParams['axes.labelsize'] = fsize

mpl.rcParams['xtick.labelsize'] = fsize-4
mpl.rcParams['ytick.labelsize'] = fsize-4

# These settings ensure that textboxes are readable in illustrator once exported
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#Set Border Width
mpl.rcParams['axes.linewidth'] = 2

#Tick Mark Settings
mpl.rcParams['xtick.major.size'] = thickness
mpl.rcParams['xtick.major.width'] = thickness
mpl.rcParams['ytick.major.size'] = thickness
mpl.rcParams['ytick.major.width'] = thickness

# Set Font
mpl.rcParams['font.sans-serif'] = 'Helvetica'

# Set color-blindness friendly color palette
mpl.rcParams['axes.prop_cycle'] = cycler(color=['#0072B2', '#D55E00', '#009E73', '#CC79A7','darkgrey', '#56B4E9','#E69F00','#F0E442']) # can add black: '#000000'

def plot_vplot(vplot_data: np.ndarray):
    """
    Generate a v_plot pdf that contains the associated cleavage probability. Count data is min-max normalized
    
    Parameters
    ----------
    vplot_data : np.ndarray
    2D Numpy array that contains the data to be plotted.

    Returns
    -------
    Generates .pdf
    """

    # Load cleavage probability array 
    # Add functionilty to return error if file does not exist
    clave_prob = np.load('cleavage_prob.npy')
    min_range = -1. * distance_from_frag_center
    max_range = distance_from_frag_center
    midpoint = fiber_midpoint
    # Subset cleavage probability array to values to plot
    y_vals = clave_prob[int(min_range)+int(midpoint):int(max_range)+int(midpoint)]
    x_vals = np.linspace(min_range, max_range, len(y_vals))

    # Process v-plot count data for plotting
    # imshow orients plots differently than the 2D histogram code used to generate array
    # Here we reorient the data for ploting with imshow
    vplot_data_rotated = np.flip(vplot_data.T)
    # Min-max normalize counts
    min_max_scaler = preprocessing.MinMaxScaler()
    array_to_plot = min_max_scaler.fit_transform(vplot_data_rotated)


    # Plot 
    fig, ax = plt.subplots(2,1,figsize = (9,6), gridspec_kw={'height_ratios': [1, 4]}, sharex=True)
    
    # First Subplot: Cleavage probability
    ax[0].plot(x_vals, y_vals)
    ax[0].set_ylabel('Cleavage \n Prob.')
    ax[0].set_xlim((min_range,max_range))
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    # Second Subplot: V-Plot
    im = ax[1].imshow(array_to_plot, cmap='seismic', extent = [min_range, max_range, 0, max_fragment_length])
    ax[1].set_xlabel('Distance from Fiber Midpoint to Fragment Center (nt)')
    ax[1].set_ylabel('Fragment \n Length (nt)')
    # Position and format color bar
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.22, 0.025, 0.33])
    cbar = fig.colorbar(im, ax=ax[1], cax=cbar_ax)#shrink=0.35)
    cbar.set_label('Relative Counts')
    # plt.tight_layout()
    plt.savefig('vplot_w_cleavage_prob.pdf')
    return


def plot_fld(fld: np.ndarray = None):
    """
    Generate a fragment length distribution .pdf
    
    Parameters
    ----------
    fld : np.ndarray, default None
    1D Numpy array that contains the fragment length resulting from all trials.
    Load saved 'frag_lens.npy' file if no fld is passed.

    Returns
    -------
    Generates .pdf
    """

    frag_lens = np.load('frag_lens.npy') # TO DO write if statement if the file does not exist.
    stored_histogram = sns.histplot(data=frag_lens, binwidth=10, stat= 'probability')
    plt.title('Fragment Length Distribution')
    plt.xlabel('Fragment Length (nt)')
    plt.tight_layout()
    plt.savefig('fld.pdf')
    plt.show() 
    return

def plot_composite(vplot_data: np.ndarray, fld: np.ndarray = None):