#! /usr/bin/env python

from __future__ import print_function

import pdb
import numpy as np
from astropy.io import fits as pf
from six import string_types
import os
import argparse

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

class DataDiff(object):
    def __init__(self, filename, avals, bvals):
        self.filename = filename
        self.avals = avals
        self.bvals = bvals

    def diff_keys(self):
        for ext in self.avals.keys():
            if not "Keywords" in self.avals[ext].keys():
                print("No keywords in {0} Extension {1}".format(self.filename, ext))
                continue 

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def diff_tables(fd_dict, dir1, dir2):
    """
    Calculate quantitative differences in HDU tables.  

    Arguments:
        fd_dict (dict): Dictionary with fitsdiff info.
        dir1 (str): Directory 1 of data to compare.
        dir2 (str): Directory 2 of data to compare.

    Returns:
        fd_dict (dcit): Updated dictionary with fitsdiff info including 
            quantitative table differences.
    """
    
    for filename in fd_dict.keys():
        opened = False
        for str_ext in fd_dict[filename]:
            ext = int(str_ext[-1])
            if "Columns" in fd_dict[filename][str_ext]:
                if not opened:
                    file1 = os.path.join(dir1, filename)
                    file2 = os.path.join(dir2, filename)
                    hdu1 = pf.open(file1)
                    hdu2 = pf.open(file2)
                    hdr1 = hdu1[1].header
                    opened = True

                # If the exposure time is 0, don't continue!
                try:
                    if hdr1["exptime"] == 0.0:
                        continue
                except KeyError:
                    pass
                
                for column in fd_dict[filename][str_ext]["Columns"].keys():
                    if len(hdu1[ext].data[column]) != len(hdu2[ext].data[column]):
                        print("SOMETHING BAD HAPPENED IN DIFF_TABLES")
                        continue
                    diff = hdu1[ext].data[column] - hdu2[ext].data[column]
                    fd_dict[filename][str_ext]["Columns"][column]["MaxDiff"] = max(diff.flatten())
                    fd_dict[filename][str_ext]["Columns"][column]["MinDiff"] = min(diff.flatten())
                    fd_dict[filename][str_ext]["Columns"][column]["MedDiff"] = np.median(diff.flatten())
                    fd_dict[filename][str_ext]["Columns"][column]["MeanDiff"] = np.average(diff.flatten())
                    percdiff = max(abs((hdu2[ext].data[column].flatten() - hdu1[ext].data[column].flatten())/((hdu2[ext].data[column].flatten() + hdu2[ext].data[column].flatten())/2.)))
                    fd_dict[filename][str_ext]["Columns"][column]["PercDiff"] = percdiff 
                    fd_dict[filename][str_ext]["Columns"][column]["NumDiff"] = len(np.where(abs(diff.flatten()) > 0.0))
                    fd_dict[filename][str_ext]["Columns"][column]["Mean1"] = np.average(hdu1[ext].data[column])
                    fd_dict[filename][str_ext]["Columns"][column]["Mean2"] = np.average(hdu2[ext].data[column])

        if opened:        
            hdu1.close()
            hdu2.close()

    return fd_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def diff_images(fd_dict, dir1, dir2):
    """
    Calculate quantitative differences in counts or flt images. 

    Arguments:
        fd_dict (dict): Dictionary with fitsdiff info.
        dir1 (str): Directory 1 of data to compare.
        dir2 (str): Directory 2 of data to compare.

    Returns:
        fd_dict (dcit): Updated dictionary with fitsdiff info including 
            quantitative image differences.
    """

    for filename in fd_dict.keys():
        opened = False
        for str_ext in fd_dict[filename].keys():
            ext = int(str_ext[-1])
            if "ImageDiff" in fd_dict[filename][str_ext].keys():
                if not opened:
                    file1 = os.path.join(dir1, filename)
                    file2 = os.path.join(dir2, filename)
                    hdu1 = pf.open(file1)
                    hdu2 = pf.open(file2)
                    hdr1 = hdu1[1].header
                    opened = True
            
                if hdr1["exptime"] == 0.0:
                    continue
            
                if np.shape(hdu1[ext].data) != np.shape(hdu2[ext].data):
                    print("SOMETHING BAD HAPPENED IN DIFF_TABLES")
                    continue
           
                diff = hdu1[ext].data - hdu2[ext].data
                fd_dict[filename][str_ext]["ImageDiff"]["MaxDiff"] = max(diff.flatten())
                fd_dict[filename][str_ext]["ImageDiff"]["MinDiff"] = min(diff.flatten())
                fd_dict[filename][str_ext]["ImageDiff"]["MedDiff"] = np.median(diff.flatten())
                fd_dict[filename][str_ext]["ImageDiff"]["MeanDiff"] = np.average(diff.flatten())
                fd_dict[filename][str_ext]["ImageDiff"]["PercDiff"] = max(abs((hdu2[ext].data.flatten() - hdu1[ext].data.flatten())/((hdu2[ext].data.flatten() + hdu2[ext].data.flatten())/2.)))
                fd_dict[filename][str_ext]["ImageDiff"]["NumDiff"] = len(np.where(abs(diff.flatten()) > 0.0))
                fd_dict[filename][str_ext]["ImageDiff"]["Mean1"] = np.average(hdu1[ext].data)
                fd_dict[filename][str_ext]["ImageDiff"]["Mean2"] = np.average(hdu2[ext].data)

        if opened:
            hdu1.close()
            hdu2.close()

    return fd_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def parse_fd(logfile):
    """
    Parse fitsdiff output, returning a dictionary sorted by filename.

    Arguments:
        logfile (str): Path of text output from fitsdiff.

    Returns:
        fd_dict (dict): Dictionary with fitsdiff info.
        dir1 (str): Directory 1 of data being compared.
        dir2 (str): Directory 2 of data being compared.
    """
    
    filename = "dummy"
    fd_dict = {}
    lines = open(logfile).read().splitlines()
    dir1 = os.path.dirname(lines[2].split()[1])
    dir2 = os.path.dirname(lines[3].split()[1])
    
    for i in range(len(lines)-1):
        # First of all, determine the filename.
        if dir1 in lines[i]:
        #if "a: " in lines[i]:
            full_filename = lines[i].split()[1]
            # If file doesn't exist, it migiht be because it's been zipped.
            if not os.path.exists(full_filename):
                full_filename += ".gz"
                if not os.path.exists(full_filename):
                    print("WARNING: File does not exist (zipped or unzipped): {0}".format(full_filename))
            filename = os.path.basename(full_filename)
            continue
         
        # Skip trailer files.    
        if "trl.fits" in filename:
            continue
        
        # Keep track of which FITS extension we're looking at. 
        if "HDU" in lines[i]:
            if "Primary HDU" in lines[i]:
                ext = "Ext0"
            else:
                extnum = lines[i].split()[-1][0]
                ext = "Ext" + str(extnum)
            if not filename in fd_dict.keys():
                fd_dict[filename] = {}
            fd_dict[filename][ext] = {}
            continue

        # Check for keywords that are only present in one file
        if "Extra keyword" in lines[i]:
            if "Keywords" not in fd_dict[filename][ext].keys():
                fd_dict[filename][ext]["Keywords"] = {}
            keyword = lines[i].split("'")[1]
            extra_file = lines[i].split()[4]
            value = lines[i].split(": ")[-1]
            if value.startswith("'") and value.endswith("'"):
                value = value.split("'")[1]
            if extra_file == "a:":
                fd_dict[filename][ext]["Keywords"][keyword] = [value, "#NOT PRESENT#"]
            else:
                fd_dict[filename][ext]["Keywords"][keyword] = ["#NOT PRESENT#", value]

        # Check for keywords that differ between the two files. 
        if "Keyword " in lines[i]:
            if "Keywords" not in fd_dict[filename][ext].keys():
                fd_dict[filename][ext]["Keywords"] = {}
           
 
        # If there is only info for one file, e.g.:
#     Keyword GLOBRT_A has different comments:
#        b> global count rate
#     Keyword GLOBRT_B has different comments:
            if "comments" in lines[i]:
                if "Comments" not in fd_dict[filename][ext].keys():
                    fd_dict[filename][ext]["Comments"] = {}
                if "        b>" not in lines[i+2]:
                    present_file = lines[i].split()[0]
                    present_comm = lines[i+1].split("> ")[1]
                    if present_file == "b>":
                        bval = present_comm
                        aval = "#NOT PRESENT#"
                    else:
                        aval = present_comm
                        bval = "#NOT PRESENT#"
                else:
                    aval = lines[i+1].split("> ")[1]
                    bval = lines[i+2].split("> ")[1]

            # Given a keyword that differs, loop through lines until you
            # find the a and b values. Determine if they are floats or not.
            else:
                keyword = lines[i].split()[1]
                j = 1
                while "a>" not in lines[i+j]:
                    j += 1
                else:
                    try:
                        aval = float(lines[i+j].split("a> ")[1])
                    except ValueError:
                        aval = lines[i+j].split("a> ")[1]
                while "b>" not in lines[i+j]:
                    j += 1
                else:
                    try:
                        bval = float(lines[i+j].split("b> ")[1])
                    except ValueError:
                        bval = lines[i+j].split("b> ")[1]
                try:
                    fd_dict[filename][ext]["Keywords"][keyword] = [aval,bval]
                except:
                    print("ruh roh")
                continue
        
        # If Columns differ, make a key for it and it will be quantified
        # later in diff_table().
        if "Column" in lines[i]:
            col = lines[i].split()[1]
            if not "Columns" in fd_dict[filename][ext].keys():
                fd_dict[filename][ext]["Columns"] = {}
            elif col not in fd_dict[filename][ext]["Columns"].keys():
                fd_dict[filename][ext]["Columns"][col] = {}
            continue
        
        # If images differ, make a key for it and it will be quantified
        # later in diff_image().
        if "pixels" in lines[i]:
            fd_dict[filename][ext]["ImageDiff"] = {}
        
    return fd_dict, dir1, dir2

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def get_diffs(logfile):
    """
    Run workhorse functions. 

    Arguments:
        logfile (str): Path of text output from fitsdiff.

    Returns:
        fd_dict (dict): Dictionary with fitsdiff info.
        dir1 (str): Directory 1 of data being compared.
        dir2 (str): Directory 2 of data being compared.
    """
    
    fd_dict, dir1, dir2 = parse_fd(logfile) 
    fd_dict = diff_tables(fd_dict, dir1, dir2)
    fd_dict = diff_images(fd_dict, dir1, dir2)
    
    return fd_dict, dir1, dir2

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", dest="logfile", default="new_retrieval.log",
                        help="Name of output ASCII files")
    args = parser.parse_args()
    logfile = args.logfile
   
    fd_dict, dir1, dir2 = get_diffs(logfile) 
