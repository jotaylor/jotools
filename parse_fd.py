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
#            for keyword in self.avals[ext]["Keywords"]:


#avals = {"Ext0": {"Keywords": {a: aval, b: bval, c: cval...},
#                  "Image": [....]
#                  "Columns": {a: ..., b: ..., c: ...}},
#                  ...}
# For each loop of parse_fd it should output a list that looks like
#  [filename, {}, {}] where the first item is the filename, the 2nd is a dict
#  of all a vals, the 3rd a dict of all b vals.
#  I would feed this list into my class DataDiff.
#    .a[ext]...
#    .b[ext]...
#  each DataDiff object gets curated into master list of all DataDiff objects
#  The I would have methods that do things to each object like give max diffs 
#  etc.

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#def diff_tables(filenames, columns, dir1, dir2):
#    # Ensure everything is in list form.
#    if isinstance(filenames, string_types):
#        filenames = [filenames]
#    if isinstance(columns[0], string_types):
#        columns = [columns]
#
#    for i in range(len(filenames)):
#        file1 = os.path.join(dir1, filenames[i])
#        file2 = os.path.join(dir2, filenames[i])
#        hdu1 = pf.open(file1)
#        hdu2 = pf.open(file2)
#        hdr1 = hdu1[1].header
#        
#        # If the exposure time is 0, don't continue!
#        try:
#            if hdr1["exptime"] == 0.0:
#                continue
#        except KeyError:
#            pass
#        
#        for ext in range(1,len(hdu1)):
#            for column in columns[i]:
#                if column in hdu1[ext].data.columns.names:
#                    if not any(badtype in hdu1[ext].columns[column].format for badtype in ["L","A"]):
#                        print(hdu1[ext].columns[column].format)
#                        diff = hdu1[ext].data[column].flatten() - hdu2[ext].data[column].flatten()
#                        maxdiff = max(abs(diff))
#                        print("Max diff for {0:26s} column {1:21s} is {2}".format(filenames[i],column,maxdiff))
#                    else:
#                        for j in range(len(hdu1[ext].data[column].flatten())):
#                            if hdu1[ext].data[column].flatten()[j] != hdu2[ext].data[column].flatten()[j]:
#                                print("Max diff for {0:26s} column {1:21s} is {2} vs. {3}".format(filenames[i],column,hdu1[ext].data[column].flatten()[j],hdu2[ext].data[column].flatten()[j]))
#
#        hdu1.close()
#        hdu2.close()

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def diff_tables(fd_dict, dir1, dir2):
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

        if opened:        
            hdu1.close()
            hdu2.close()

    return fd_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#def diff_images(filenames, dir1, dir2):
#    if isinstance(filenames, string_types):
#        filenames = [filenames]
#    
#    files_seen = []
#    for i in range(len(filenames)):
#        if filenames[i] not in files_seen:
#            files_seen.append(filenames)
#        else:
#            continue
#        
#        file1 = os.path.join(dir1, filenames[i])
#        file2 = os.path.join(dir2, filenames[i])
#        hdu1 = pf.open(file1)
#        hdu2 = pf.open(file2)
#        hdr1 = hdu1[1].header
#        if hdr1["exptime"] == 0.0:
#            continue
#        
#        for ext in range(1,len(hdu1)):
#            data1 = hdu1[ext].data.flatten()
#            data2 = hdu2[ext].data.flatten()
#            diff = abs(data2 - data1)
#            maxdiff = max(diff)
#            percdiff = max(abs((data2 - data1)/((data2 + data1)/2.)))
#            inds = np.where(diff > 0.0)
#            num_diff = len(inds)
#            if maxdiff != 0.0:
#                print("Max diff for {:26s} {} is {}".format(filenames[i],ext,maxdiff))
#                print("Max % diff for {:26s} {} is {}".format(filenames[i],ext,percdiff))
#                print("    Loc of diff elems is {}".format(inds)) 
#                print("    Number of diff elems is {}\n".format(num_diff))
#        
#
##        orig1 = hdu1[1].data
##        orig2 = hdu2[1].data
##        data1_sci = hdu1[1].data.flatten()
##        data2_sci = hdu2[1].data.flatten()
##        diff = data2_sci - data1_sci
##        diff2 = abs(orig2 - orig1)
##        maxdiff = max(abs(diff))
##        percdiff = max(abs((data2_sci - data1_sci)/((data2_sci + data1_sci)/2.)))
##        inds = np.where(diff2 > 0.0)
##        num_diff = len(inds)
##        if maxdiff != 0.0:
##            print("Max diff for {:26s} SCI is {}".format(filenames[i],maxdiff))
##            print("Max % diff for {:26s} SCI is {}".format(filenames[i],percdiff))
##            print("    Loc of diff elems is {}".format(inds)) 
##            print("    Number of diff elems is {}\n".format(num_diff))
##        
##        data1_err = hdu1[2].data.flatten()
##        data2_err = hdu2[2].data.flatten()
##        orig1 = hdu1[2].data
##        orig2 = hdu1[2].data
##        diff = data2_err - data1_err
##        diff2 = abs(orig2 - orig1)
##        maxdiff = max(abs(diff))
##        percdiff = max(abs((data2_err - data1_err)/((data2_err + data1_err)/2.)))
##        inds = np.where(abs(diff2) > 0.0)
##        num_diff = len(inds)
##        if maxdiff != 0.0:
##            print("Max diff for {:26s} ERR is {}".format(filenames[i],maxdiff))
##            print("Max percent diff for {:26s} ERR is {}".format(filenames[i],percdiff))
##            print("    Loc of diff elems is {}".format(inds)) 
##            print("    Number of diff elems is {}\n".format(num_diff))
##
##        data1_dq = hdu1[3].data.flatten()
##        data2_dq = hdu2[3].data.flatten()
##        orig1 = hdu1[3].data
##        orig2 = hdu1[3].data
##        diff = data2_dq - data1_dq
##        diff2 = abs(orig2 - orig1)
##        maxdiff = max(abs(diff))
##        percdiff = max(abs((data2_dq - data1_dq)/((data2_dq + data1_dq)/2.)))
##        inds = np.where(abs(diff2) > 0.0)
##        num_diff = len(inds)
##        if maxdiff != 0.0:
##            print("Max diff for {:26s} DQ is {}".format(filenames[i],maxdiff))
##            print("Max percent diff for {:26s} DQ is {}".format(filenames[i],percdiff))
##            print("    Loc of diff elems is {}".format(inds)) 
##            print("    Number of diff elems is {}\n".format(num_diff))
#        hdu1.close()
#        hdu2.close()        
        
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
def diff_images(fd_dict, dir1, dir2):
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

        if opened:
            hdu1.close()
            hdu2.close()

    return fd_dict

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def parse_fd(logfile):
    filename = "dummy"
    fd_dict = {}
    lines = open(logfile).read().splitlines()
    dir1 = os.path.dirname(lines[2].split()[1])
    dir2 = os.path.dirname(lines[3].split()[1])
    for i in range(len(lines)-1):
        # First of all, determine the filename.
        if dir1 in lines[i]:
        #if "a: " in lines[i]:
            filename = os.path.basename(lines[i].split()[1])
            continue
         
        # Skip trailer files    
        if "trl.fits" in filename:
            continue
        
        # Keep track of which extension we're looking at    
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
        
        if "Column" in lines[i]:
            col = lines[i].split()[1]
            if not "Columns" in fd_dict[filename][ext].keys():
                fd_dict[filename][ext]["Columns"] = {}
            elif col not in fd_dict[filename][ext]["Columns"].keys():
                fd_dict[filename][ext]["Columns"][col] = {}
            continue
        if "pixels" in lines[i]:
            fd_dict[filename][ext]["ImageDiff"] = {}
        
    return fd_dict, dir1, dir2

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def get_diffs(logfile):
    fd_dict, dir1, dir2 = parse_fd(logfile) 
#    for filename in fd_dict.keys():
#        for ext in fd_dict[filename].keys():
#            if "Columns" in fd_dict[filename][ext].keys():
#                diff_tables(filename, fd_dict[filename][ext]["Columns"], dir1, dir2)
#            elif "ImageDiff" in fd_dict[filename][ext].keys():
#                diff_images(filename, dir1, dir2) 

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
   
    get_diffs(logfile) 
