
"""
Includes functions for plotting the fragment length distribution, vplot, and cleavage probability
"""

# Import Modules 
from .fragment_lengths import *
from .params import max_fragment_length, distance_from_frag_center, fiber_midpoint, nrl, num_nucs, break_rate
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

def process_vplot_data(vplot_data: np.ndarray):
    """
    Normalize vplot data and rotate the matrix for later plotting with plt.imshow()
    
    Parameters
    ----------
    vplot_data : np.ndarray
    2D Numpy array that contains the data to be plotted.

    Returns
    -------
    array_to_plot : np.ndarray
    Normalized and rotated vplot_data
    """

    # Process v-plot count data for plotting
    # imshow orients plots differently than the 2D histogram code used to generate array
    # Here we reorient the data for ploting with imshow
    vplot_data_rotated = np.flip(vplot_data.T)
    # Min-max normalize counts
    min_max_scaler = preprocessing.MinMaxScaler()
    array_to_plot = min_max_scaler.fit_transform(vplot_data_rotated)
    np.save('intermed_data/vplot_norm.npy', array_to_plot) # This feature is new
    return array_to_plot


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
    cleave_prob = np.load('intermed_data/cleavage_prob.npy')
    min_range = -1. * distance_from_frag_center
    max_range = distance_from_frag_center
    midpoint = fiber_midpoint
    # Subset cleavage probability array to values to plot
    y_vals = cleave_prob[int(min_range)+int(midpoint):int(max_range)+int(midpoint)]
    x_vals = np.linspace(min_range, max_range, len(y_vals))

    array_to_plot = process_vplot_data(vplot_data)  

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
    fig.suptitle(str(num_nucs) + " Nucleosome Fiber with "+str(nrl)+" NRL; 1 break per "+str(break_rate)+" nucleotides")
    plt.savefig('plots/vplot_w_cleavage_prob.pdf')
    plt.show() #for running in jupyter notebook 
    plt.close()
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

    TO DO
    -----
    subset to only the fragment lengths less than max considered.
    """

    frag_lens_pre = np.load('intermed_data/frag_lens.npy') # TO DO write if statement if the file does not exist.
    # Subset to relevant fragments
    frag_lens = frag_lens_pre[frag_lens_pre < max_fragment_length]
    fig, ax = plt.subplots()
    stored_histogram = sns.histplot(data=frag_lens, binwidth=10, stat= 'probability')
    plt.title('Fragment Length Distribution')
    plt.xlabel('Fragment Length (nt)')
    # plt.tight_layout()
    fig.suptitle(str(num_nucs) + " Nucleosome Fiber with "+str(nrl)+" NRL; 1 break per "+str(break_rate)+" nucleotides")
    plt.savefig('plots/fld.pdf')
    plt.show() 
    plt.close()
    return

def plot_composite(vplot_data: np.ndarray):#, fld: np.ndarray = None):
    """
    Generate a v_plot pdf that contains the associated cleavage probability. Count data is min-max normalized
    
    Parameters
    ----------
    vplot_data : np.ndarray
    2D Numpy array that contains the data to be plotted.

    Returns
    -------
    Generates .pdf

    TO DO
    -----
    - Improve handling of reading in fld
    - Need to scale FLD and vplot to accomodate max_fragment_sizes other than 1000
    """

    # Load cleavage probability array 
    # Add functionilty to return error if file does not exist
    clave_prob = np.load('intermed_data/cleavage_prob.npy')
    min_range = -1. * distance_from_frag_center
    max_range = distance_from_frag_center
    midpoint = fiber_midpoint
    # Subset cleavage probability array to values to plot
    y_vals = clave_prob[int(min_range)+int(midpoint):int(max_range)+int(midpoint)]
    x_vals = np.linspace(min_range, max_range, len(y_vals))

    array_to_plot = process_vplot_data(vplot_data)  

    # Plot 
    # fig, ax = plt.subplots(2,1,figsize = (9,6), gridspec_kw={'height_ratios': [1, 4]}, sharex=True)
    fig, (ax1, ax2) = plt.subplots(2,2,figsize = (11,6), gridspec_kw={'height_ratios': [1, 4], 'width_ratios': [9, 2]})
    # First Subplot: Cleavage probability
    ax1[0].plot(x_vals, y_vals)
    ax1[0].set_ylabel('Cleavage \n Prob.')
    ax1[0].set_xlim((min_range,max_range))
    ax1[0].set_ylim((0,1))
    # Get rid of unwanted axes
    ax1[0].spines['right'].set_visible(False)
    ax1[0].spines['top'].set_visible(False)
    ax2[1].spines['right'].set_visible(False)
    ax2[1].spines['top'].set_visible(False)
    #Clear top right plot
    ax1[1].clear()
    ax1[1].axis("off")
    ax1[1].set_visible(False)
    ax1[1].remove()

    # Second Subplot: V-Plot
    im = ax2[0].imshow(array_to_plot, cmap='seismic', extent = [min_range, max_range, 0, max_fragment_length], aspect = 1)
    ax2[0].set_xlabel('Distance from Fiber Midpoint to Fragment Center (nt)')
    ax2[0].set_ylabel('Fragment \n Length (nt)')
    # Position and format color bar
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.22, 0.025, 0.33])
    cbar = fig.colorbar(im, ax=ax2[0], cax=cbar_ax)#shrink=0.35)
    cbar.set_label('Relative Counts')
    # Set shared axes
    ax1[0].sharex(ax2[0])
    ax2[1].sharey(ax2[0])
    # Got scaled value by guessing and checking.
    ax2[1].set_box_aspect(1.875)

    ax2[1].set_ylabel(' ')
    # Plot fld
    frag_lens = np.load('intermed_data/frag_lens.npy') # TO DO -- account for when fld does not exist
    midpts = np.load('intermed_data/frag_midpts.npy') 
    fm = frag_mid_df(frag_lens, midpts)
    data_for_hist = pd.DataFrame(fm[fm.frag_len<max_fragment_length].frag_len)
    sns.histplot(data=data_for_hist, y = 'frag_len', bins=100, stat= 'probability', ax=ax2[1])
    fig.suptitle(str(num_nucs) + " Nucleosome Fiber with "+str(nrl)+" NRL; 1 break per "+str(break_rate)+" nucleotides")
    plt.savefig('plots/vplot_w_fld_and_cleavprob.pdf')
    plt.show()
    plt.close()

    return