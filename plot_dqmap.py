#! /usr/bin/env python
'''
Plot the DQ regions of a given BPIXTAB and/or gsagtab, with the option of only
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

def define_cmap(cmap_name="nipy_spectral"):
    # Get existing cmap colors
    cmap = pl.get_cmap(cmap_name)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    
    #Create a new segmented cmap based off the existing one.
    cmap = LinearSegmentedColormap.from_list("dq", cmaplist, cmap.N)

    # Defint the boundaries for each color, in this case halfway between DQs
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

def get_reffile_info(bpixtab, gsagtab, sdq, sdqflags, date, hva, hvb):
    if bpixtab and gsagtab:
        bpix_dict = parse_bpixtab(bpixtab, sdq, sdqflags)
        gsag_dict = parse_gsagtab(gsagtab, sdq, sdqflags, date, hva, hvb)
        data_dict = {}
        for segment in bpix_dict.keys():
            data_dict[segment] = {} 
            for key in bpix_dict[segment]:
                data_dict[segment][key] = list(bpix_dict[segment][key])+list(gsag_dict[segment][key])
#        data_dict = {key:bpix_dict[key]+gsag_dict[key] for key in bpix_dict.keys()}
    elif bpixtab:
        data_dict = parse_bpixtab(bpixtab, sdq, sdqflags)
    elif gsagtab:
        data_dict = parse_gsagtab(gsagtab, sdq, sdqflags, date, hva, hvb)
    else:
        print("You must specify either a BPIXTAB or GSAGTAB!")

    return data_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def parse_bpixtab(bpixtab, sdq, sdqflags):
    print("Using {0}".format(bpixtab))
    with pf.open(bpixtab) as hdulist:
        data = hdulist[1].data
    segments = list(set(data["segment"]))
    
    data_dict = {}
    for i in np.arange(len(segments)):
        data_dict[segments[i]] = {"lx": [], 
                                  "ly": [],
                                  "dx": [],
                                  "dy": [],
                                  "dq": []}
        inds = np.where(data["segment"] == segments[i])[0]
        
        if sdq:
            inds = np.where(data[inds]["dq"]&sdqflags != 0)[0]
        
        for col in ["lx", "ly", "dx", "dy", "dq"]:
            data_dict[segments[i]][col] = data[inds][col]

    return data_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def parse_gsagtab(gsagtab, sdq, sdqflags, date, hva, hvb):
    print("Using {0}".format(gsagtab))
    hdulist = pf.open(gsagtab)

    if not hva:
        hva = 167
        print("WARNING: no HVLEVELA specified, using default value of {0}".format(hva))
    if not hvb:
        hvb = 175
        print("WARNING: no HVLEVELA specified, using default value of {0}".format(hvb))
    
    exts = {}
    i = 1
    while len(exts.keys()) != 2:
        for key in ["FUVA", "FUVB"]:
            if key not in exts.keys():
                try:
                    h_key = "hvlevel{0}".format(key[-1].lower())
                    hv = hdulist[i].header[h_key]
                    exts[key] = [i,hv]
                except KeyError:
                    pass
                except IndexError:
                    print("Something went wrong, only found {0}".format(exts))
                    break
        i += 1
    
    if not date:
        date = 57822.0
        print("WARNING: no DATE specified, using default value of {0}".format(str(date)))

    data_dict = {}
    for segment in exts.keys():
        data_dict[segment] = {"lx": [], 
                                 "ly": [],
                                 "dx": [],
                                 "dy": [],
                                 "dq": []}
        ext = exts[segment][0]
        data = hdulist[ext].data
        inds = np.where(data["date"] <= date)[0]
        
        if sdq:
            inds = np.where(data[inds]["dq"]&sdqflags != 0)[0]
        
        for col in ["lx", "ly", "dx", "dy", "dq"]:
            data_dict[segment][col] = data[inds][col]

    return data_dict 
                        
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def plot_reffile(data_dict, bpixtab, gsagtab, sdq, sdqflags):
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
            refname = "{0} and {1}".format(bpixtab, gsagtab)
        elif bpixtab:
            refname = bpixtab
        elif gsagtab:
            refname = gsagtab
        ax.set_title("{0} {1} DQ Flags".format(refname, segment))
        ax.set_ylabel("YFULL")
        
        # Print SDQFLAGS value if applicable.
        if sdq:
            ax.annotate("SDQFLAGS={0}".format(sdqflags), 
                        xy=(0.01, 0.9),
                        xycoords="axes fraction", 
                        fontsize=14)
        
        # Set axis labels
    ax_list = fig.axes
    ax_list[-1].set_xlabel("XFULL")
    
    # Create an axis for the colorbar so it encompasses *all* subplots.
    cax = pl.axes([0.92, 0.2, 0.02, 0.6])
    cbar = fig.colorbar(p, cax=cax, boundaries=bounds, ticks=DQs)
    cax.set_ylabel("DQ Flag")
    
    pl.show()
    this = input("press enter to continue")

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
                        default=False, help="Switch to plot only SDQ")
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
    
    args = parser.parse_args()
    bpixtab = args.bpixtab
    gsagtab = args.gsagtab
    sdq = args.sdq
    sdqflags = int(args.sdqflags)
    date = args.date
    hva = args.hva
    hvb = args.hvb

    # Define colormap
    cmap, bounds, norm = define_cmap()

    data_dict = get_reffile_info(bpixtab, gsagtab, sdq, sdqflags, date, hva, hvb)
    plot_reffile(data_dict, bpixtab, gsagtab, sdq, sdqflags)
