#! /usr/bin/env python
'''
Plot the DQ regions of a given BPIXTAB and/or GSAGTAb, with the option of only
plotting SDQ values.

Usage:
To only plot BPIXTAB:
    python plot_dqmap.py -bpix yae1249sl_bpix.fits

To only plot BPIXTAB SDQ regions:
    python plot_dqmap.py -bpix yae1249sl_bpix.fits -sdq

To only plot BPIXTAB SDQ regions with custom SDQFLAGS value:
    python plot_dqmap.py -bpix yae1249sl_bpix.fits -sdq -sdqflags 64

To plot BPIXTAB and GSAGTAB:
    python plot_dqmap.py -bpix yae1249sl_bpix.fits -gsag zbn1927gl_gsag.fits
'''

__author__ = "Jo Taylor"
__date__ = "03-10-2017"
__maintainer__ = "Jo Taylor"
__email__ = "jotaylor@stsci.edu"

import argparse
import matplotlib
from matplotlib import pyplot as pl
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
from astropy.io import fits as pf
import numpy as np

DQs = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384])

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def define_cmap(cmap_name="gist_ncar"):
    '''
    Define a custom linear colormap, based off an existing matplotlib colormap.

    Parameters:
    -----------
        cmap_name : str
            Colormap to modify for use. Default is nipy_spectral.

    Returns:
    --------
        cmap : matplotlib.colors.LinearSegmentedColormap object
            Modified cmap for use.
        bounds : array-like object
            Boundaries of the cmap.
        norm : matplotlib.colors.BoundaryNorm object
            Normalization of cmap
    '''
    
    # Get existing cmap colors
    cmap = pl.get_cmap(cmap_name)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    
    # Create a new segmented cmap based off the existing one.
    cmap = LinearSegmentedColormap.from_list("dq", cmaplist, cmap.N)

    # Define the boundaries for each color, in this case halfway between DQs
    b = [0]
    for i in np.arange(len(DQs)):
        if i == len(DQs)-1:
            b.append(DQs[i]+1000)
            break
        b.append((DQs[i+1]+DQs[i])/2)
    bounds = np.array(b)

    # Normalize the cmap 
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    return cmap, bounds, norm

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def get_reffile_info(bpixtab, gsagtab, sdq, sdqflags, date, hvs):
    '''
    Get the relevant information from BPIXTAB and/or GSAGTAB.

    Parameters:
    -----------
        bpixtab : str
            The BPIXTAB to plot.
        gsagtab : str
            The GSAGTAB to plot.
        sdq : Bool
            Switch to plot only SDQ regions.
        sdqflags : int
            Serious Data Quality Flags
        date : float
            MJD date to use when selecting DQ regions in the GSAGTAB.
        hvs: dict
            Dictionary with information on HVLEVELA & HVLEVELB for use
            when selecting the correct GSAGTAB extension.
    
    Returns:
    --------
        data_dict : dictionary
            Dictionary with information describing the LX, LY, DX, DY, and DQ
            for each segment. 
        date : float
            Updated MJD date to use when selecting DQ regions in the GSAGTAB.
        hvs: dict
            Updated dictionary with information on HVLEVELA & HVLEVELB for use
            when selecting the correct GSAGTAB extension.
    '''

    # If both ref files, parse them both and combine results into one dict.
    if bpixtab and gsagtab:
        bpix_dict = parse_bpixtab(bpixtab, sdq, sdqflags)
        gsag_dict, hvs, date = parse_gsagtab(gsagtab, sdq, sdqflags, date, hvs)
        data_dict = {}
        for segment in bpix_dict.keys():
            data_dict[segment] = {} 
            for key in bpix_dict[segment]:
                data_dict[segment][key] = list(bpix_dict[segment][key])+list(gsag_dict[segment][key])
#        data_dict = {key:bpix_dict[key]+gsag_dict[key] for key in bpix_dict.keys()}
    elif bpixtab:
        data_dict = parse_bpixtab(bpixtab, sdq, sdqflags)
    elif gsagtab:
        data_dict, hvs, date = parse_gsagtab(gsagtab, sdq, sdqflags, date, hvs)
    else:
        print("You must specify either a BPIXTAB or GSAGTAB!")

    return data_dict, hvs, date

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def parse_bpixtab(bpixtab, sdq, sdqflags):
    '''
    Open a BPIXTAB and get the DQ information for each segment.

    Parameters:
    -----------
        bpixtab : str
            The BPIXTAB to plot.
        sdq : Bool
            Switch to plot only SDQ regions.
        sdqflags : int
            Serious Data Quality Flags
    
    Returns:
    --------
        data_dict : dictionary
            Dictionary with information describing the LX, LY, DX, DY, and DQ
    '''

    print("Using {0}".format(bpixtab))
    # Open file and get the 1st data extension.
    with pf.open(bpixtab) as hdulist:
        data = hdulist[1].data
    segments = list(set(data["segment"]))
    
    # Curate dictionary with DQ information for each segment.
    data_dict = {}
    for segment in segments:
        data_dict[segment] = {"lx": [], 
                              "ly": [],
                              "dx": [],
                              "dy": [],
                              "dq": []}
        inds = np.where(data["segment"] == segment)[0]

        # Get only DQ value sin SDQFLAGS if applicable.        
        if sdq:
            inds = np.where(data[inds]["dq"]&sdqflags != 0)[0]
        
        for col in ["lx", "ly", "dx", "dy", "dq"]:
            data_dict[segment][col] = data[inds][col]

    return data_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def parse_gsagtab(gsagtab, sdq, sdqflags, date, hvs):
    '''
    Open a GSAGTAB and get the DQ information for each segment, given a
    specific date, HVLEVELA, and HVLEVELB. 
    
    Parameters:
    -----------
        gsagtab : str
            The GSAGTAB to plot.
        sdq : Bool
            Switch to plot only SDQ regions.
        sdqflags : int
            Serious Data Quality Flags
        date : float
            MJD date to use when selecting DQ regions in the GSAGTAB.
        hvs: dict
            Dictionary with information on HVLEVELA & HVLEVELB for use
            when selecting the correct GSAGTAB extension.
    
    Returns:
    --------
        data_dict : dictionary
            Dictionary with information describing the LX, LY, DX, DY, and DQ
        date : float
            Updated MJD date to use when selecting DQ regions in the GSAGTAB.
        hvs: dict
            Updated dictionary with information on HVLEVELA & HVLEVELB for use
            when selecting the correct GSAGTAB extension.
    '''

    print("Using {0}".format(gsagtab))
    hdulist = pf.open(gsagtab)

    # If user did not specify HVLEVLA/B, set to default values.
    if not hvs["FUVA"]:
        hvs["FUVA"] = 167
        print("WARNING: no HVLEVELA specified, using default value of {0}".format(hvs["FUVA"]))
    else:
        print("Using HVLEVELA = {0}".format(hvs["FUVA"]))
    if not hvs["FUVB"]:
        hvs["FUVB"] = 175
        print("WARNING: no HVLEVELA specified, using default value of {0}".format(hvs["FUVB"]))
    else:
        print("Using HVLEVELB = {0}".format(hvs["FUVB"]))

    # While loop that goes through each data extension of the GSAGTAB until it finds
    # matching values for specified HVLEVELA/B..
    exts = {}
    i = 1
    while len(exts.keys()) != 2:
        for key in ["FUVA", "FUVB"]:
            if key not in exts.keys():
                try:
                    h_key = "hvlevel{0}".format(key[-1].lower())
                    hv = hdulist[i].header[h_key]
                    if hv == hvs[key]:
                        exts[key] = [i,hv]
                # When you're looking at the wrong segment.
                except KeyError:
                    pass
                # When you've reached the end of the extensions.    
                except IndexError:
                    print("Something went wrong, only found {0} for HVLEVELA={1}, HVLEVELB={2}".format(exts, hvs["FUVA"], hvs["FUVB"]))
                    break
        # If the for loop exited normally
        else:
            i += 1
            continue
        # If the for loop did not exit normally
        break

    # If user did not specify date, set to default value.    
    if not date:
        date = 57822.0
        print("WARNING: no date specified, using default value of {0}".format(str(date)))
    else:
        print("Using date (MJD) = {0}".format(date))

    # Curate dictionary with DQ information for each segment.
    data_dict = {}
    for segment in exts.keys():
        data_dict[segment] = {"lx": [], 
                              "ly": [],
                              "dx": [],
                              "dy": [],
                              "dq": []}
        ext = exts[segment][0]
        data = hdulist[ext].data
        # Select data by date.
        inds = np.where(data["date"] <= date)[0]
        
        # Get only DQ value sin SDQFLAGS if applicable.        
        if sdq:
            inds = np.where(data[inds]["dq"]&sdqflags != 0)[0]
        
        for col in ["lx", "ly", "dx", "dy", "dq"]:
            data_dict[segment][col] = data[inds][col]

    hdulist.close()
    
    return data_dict, hvs, date 
                        
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def plot_reffile(data_dict, bpixtab, gsagtab, hvs, date, sdq, sdqflags, cmap, bounds, norm, save):
    '''
    Plot the bad pixel regions (BPIXTAB and/or GSAGTAB) for each segment.

    Parameters:
    -----------
        bpixtab : str
            The BPIXTAB to plot.
        gsagtab : str
            The GSAGTAB to plot.
        sdq : Bool
            Switch to plot only SDQ regions.
        sdqflags : int
            Serious Data Quality Flags
        hvs: dict
            Dictionary with information on HVLEVELA & HVLEVELB for use
            when selecting the correct GSAGTAB extension.
        date : float
            MJD date to use when selecting DQ regions in the GSAGTAB.
        cmap : matplotlib.colors.LinearSegmentedColormap object
            Modified cmap for use.
        bounds : array-like object
            Boundaries of the cmap.
        norm : matplotlib.colors.BoundaryNorm object
            Normalization of cmap
        save : Bool
            Switch to save figure as png
             
    Returns:
    --------
        Nothing
    '''

    # Initialize figure and define the subplots using gridspec
    fig = pl.figure(figsize=(17,12))
    pl.subplots_adjust(hspace=0.1)    
    gs = gridspec.GridSpec(2,1, height_ratios=[1,1])
    
    # Loop over each segment (could be >2 if NUV)
    # Curate list of each patch (rectangular DQ region) and color (DQ flag)).
    for i, segment in enumerate(data_dict):
        patches = []
        colors = []
        ax = pl.subplot(gs[i])
    
        for j in np.arange(len(data_dict[segment]["lx"])):
            r = Rectangle((data_dict[segment]["lx"][j], data_dict[segment]["ly"][j]),
                          data_dict[segment]["dx"][j],
                          data_dict[segment]["dy"][j],)
            patches.append(r)
            colors.append(data_dict[segment]["dq"][j])
    
        # For each segment, create a PatchCollection with all patches
        # (DQ regions) and plot it.
        p = PatchCollection(patches, edgecolors="none", cmap=cmap, norm=norm, alpha=0.6)
        p.set_array(np.array(colors))
        ax.add_collection(p) 
        ax.set_xlim(500, 15500)
        ax.set_ylim(200, 830)
        if bpixtab and gsagtab:
            refname = "{0} and {1} (date={2})".format(bpixtab, gsagtab, date)
        elif bpixtab:
            refname = bpixtab
        elif gsagtab:
            refname = "{0} (date={1})".format(gsagtab, date)
        ax.set_title("{0} {1} DQ Flags".format(refname, segment))
        ax.set_ylabel("YFULL")
        
        # Print SDQFLAGS value if applicable.
        if sdq:
            print("Using SDQFLAGS = {0}".format(sdqflags))
            ax.annotate("SDQFLAGS={0}".format(sdqflags), 
                        xy=(0.01, 0.9),
                        xycoords="axes fraction", 
                        fontsize=14)
        # Print HVLEVELA & HVELEVB if applicable
        if gsagtab:
            ax.annotate("HVLEVEL{0}={1}".format(segment[-1], hvs[segment]),
                        xy=(0.01, 0.85),
                        xycoords="axes fraction",
                        fontsize=14)

        # Set axis labels
    ax_list = fig.axes
    ax_list[-1].set_xlabel("XFULL")
    
    # Create an axis for the colorbar so it encompasses *all* subplots.
    cax = pl.axes([0.92, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(p, cax=cax, boundaries=bounds, ticks=DQs)
    cax.set_ylabel("DQ Flag")

    if save:
        figname = "bad_regions.png"
        fig.savefig(figname, bbox_inches="tight", dpi=200)
        print("Saved {0}".format(figname))
    else:   
        pl.show()

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-bpix", 
                        dest="bpixtab", 
                        default=None,
                        help="BPIXTAB to plot.")
    parser.add_argument("-gsag", 
                        dest="gsagtab", 
                        default=None,
                        help="GSAGTAB to plot.")
    parser.add_argument("-sdq", 
                        dest="sdq", 
                        action="store_true", 
                        default=False, 
                        help="Switch to plot only SDQ")
    parser.add_argument("-sdqflags", 
                        dest="sdqflags", 
                        default=8376,
                        help="Value of SDQFLAGS")
    parser.add_argument("-date",
                        dest="date",
                        default=None,
                        help="Date of interest in GSAGTAB")
    parser.add_argument("-hva",
                        dest="hva",
                        default=None,
                        help="HVLEVELA value of interest in GSAGTAB")
    parser.add_argument("-hvb",
                        dest="hvb",
                        default=None,
                        help="HVLEVELB value of interest in GSAGTAB")
    parser.add_argument("-save",
                        dest="save",
                        action="store_true",
                        default=False,
                        help="Switch to save figure.")
    
    args = parser.parse_args()
    bpixtab = args.bpixtab
    gsagtab = args.gsagtab
    sdq = args.sdq
    sdqflags = int(args.sdqflags)
    if args.date: 
        date = float(args.date)
    else:
        date = args.date
    if args.hva:
        hva = int(args.hva)
    else:
        hva = args.hva
    if args.hvb:
        hvb = int(args.hvb)
    else:
        hvb = args.hvb
    hvs = {"FUVA": hva, "FUVB": hvb}
    save = args.save

    # Define colormap
    cmap, bounds, norm = define_cmap()

    data_dict, hvs, date = get_reffile_info(bpixtab, gsagtab, sdq, sdqflags, date, hvs)
    plot_reffile(data_dict, bpixtab, gsagtab, hvs, date, sdq, sdqflags, cmap, bounds, norm, save)
