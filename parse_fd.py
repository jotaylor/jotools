#! /usr/bin/env python

from __future__ import print_function
import numpy as np

from astropy.io import fits as pf
from six import string_types
import os
import argparse

pathname = "/grp/hst/cos" 

def diff_tables(filenames, columns, dir1, dir2):
    if isinstance(filenames, string_types):
        filenames = [filenames]
    if isinstance(columns[0], string_types):
        columns = [columns]
    for i in range(len(filenames)):
        file1 = os.path.join(dir1, filenames[i])
        file2 = os.path.join(dir2, filenames[i])
        hdu1 = pf.open(file1)
        hdu2 = pf.open(file2)
        data1 = hdu1[1].data
        data2 = hdu2[1].data
        hdr1 = hdu1[1].header
        if hdr1["exptime"] == 0.0:
            continue
        for column in columns[i]:
            try:
                diff = data2[column].flatten() - data1[column].flatten()
            except KeyError:
                try:
                    diff = hdu2[2].data[column].flatten() - hdu1[2].data[column].flatten()
                except KeyError:
                    try:
                        diff = hdu2[3].data[column].flatten() - hdu1[3].data[column].flatten()
                    except KeyError:
                        print("Something went horribly wrong with {} {}".format(filenames[i],column))
                        continue
            maxdiff = max(abs(diff))
            print("Max diff for {:26s} {:21s} is {}".format(filenames[i],column,maxdiff))
        hdu1.close()
        hdu2.close()


def diff_images(filenames, dir1, dir2):
    if isinstance(filenames, string_types):
        filenames = [filenames]
    files_seen = []
    for i in range(len(filenames)):
        if filenames[i] not in files_seen:
            files_seen.append(filenames)
        else:
            continue
        file1 = os.path.join(dir1, filenames[i])
        file2 = os.path.join(dir2, filenames[i])
        hdu1 = pf.open(file1)
        hdu2 = pf.open(file2)
        orig1 = hdu1[1].data
        orig2 = hdu2[1].data
        hdr1 = hdu1[1].header
        if hdr1["exptime"] == 0.0:
            continue
        data1_sci = hdu1[1].data.flatten()
        data2_sci = hdu2[1].data.flatten()
        diff = data2_sci - data1_sci
        diff2 = abs(orig2 - orig1)
        maxdiff = max(abs(diff))
        percdiff = max(abs((data2_sci - data1_sci)/((data2_sci + data1_sci)/2.)))
        inds = np.where(diff2 > 0.0)
        num_diff = len(inds)
        if maxdiff != 0.0:
            print("Max diff for {:26s} SCI is {}".format(filenames[i],maxdiff))
            print("Max % diff for {:26s} SCI is {}".format(filenames[i],percdiff))
            print("    Loc of diff elems is {}".format(inds)) 
            print("    Number of diff elems is {}\n".format(num_diff))
        
        data1_err = hdu1[2].data.flatten()
        data2_err = hdu2[2].data.flatten()
        orig1 = hdu1[2].data
        orig2 = hdu1[2].data
        diff = data2_err - data1_err
        diff2 = abs(orig2 - orig1)
        maxdiff = max(abs(diff))
        percdiff = max(abs((data2_err - data1_err)/((data2_err + data1_err)/2.)))
        inds = np.where(abs(diff2) > 0.0)
        num_diff = len(inds)
        if maxdiff != 0.0:
            print("Max diff for {:26s} ERR is {}".format(filenames[i],maxdiff))
            print("Max percent diff for {:26s} ERR is {}".format(filenames[i],percdiff))
            print("    Loc of diff elems is {}".format(inds)) 
            print("    Number of diff elems is {}\n".format(num_diff))

        data1_dq = hdu1[3].data.flatten()
        data2_dq = hdu2[3].data.flatten()
        orig1 = hdu1[3].data
        orig2 = hdu1[3].data
        diff = data2_dq - data1_dq
        diff2 = abs(orig2 - orig1)
        maxdiff = max(abs(diff))
        percdiff = max(abs((data2_dq - data1_dq)/((data2_dq + data1_dq)/2.)))
        inds = np.where(abs(diff2) > 0.0)
        num_diff = len(inds)
        if maxdiff != 0.0:
            print("Max diff for {:26s} DQ is {}".format(filenames[i],maxdiff))
            print("Max percent diff for {:26s} DQ is {}".format(filenames[i],percdiff))
            print("    Loc of diff elems is {}".format(inds)) 
            print("    Number of diff elems is {}\n".format(num_diff))
        hdu1.close()
        hdu2.close()        
        

def parse_fd(logfile):
    fd_dict = {}
    lines = open(logfile).read().splitlines()
    dir1 = os.path.dirname(lines[2].split()[1])
    dir2 = os.path.dirname(lines[3].split()[1])
    for i in range(len(lines)):
        if dir1 in lines[i]:
        #if "a: " in lines[i]:
            filename = os.path.basename(lines[i].split()[1])
            continue
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
        if "Keyword " in lines[i]:
            keyword = lines[i].split()[1]
            j = 1
            while "a>" not in lines[i+j]:
                j += 1
            else:
                try:
                    aval = float(lines[i+j].split()[1])
                except ValueError:
                    aval = lines[i+j].split()[1]
            while "b>" not in lines[i+j]:
                j += 1
            else:
                try:
                    bval = float(lines[i+j].split()[1])
                except ValueError:
                    bval = lines[i+j].split()[1]
            try:
                fd_dict[filename][ext][keyword] = [aval,bval]
            except:
                pdb.set_trace()
            continue
        if "Column" in lines[i]:
            col = lines[i].split()[1]
            if not "Columns" in fd_dict[filename][ext].keys():
                fd_dict[filename][ext]["Columns"] = [col]
            elif col not in fd_dict[filename][ext]["Columns"]:
                fd_dict[filename][ext]["Columns"].append(col)
            continue
        if "pixels" in lines[i]:
            fd_dict[filename][ext]["ImageDiff"] = True
        
    return fd_dict, dir1, dir2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", dest="logfile", default="new_retrieval.log",
                        help="Name of output ASCII files")
    args = parser.parse_args()
    logfile = args.logfile
    
    fd_dict, dir1, dir2 = parse_fd(logfile) 
    for filename in fd_dict.keys():
        for ext in fd_dict[filename].keys():
            if "Columns" in fd_dict[filename][ext].keys():
                diff_tables(filename, fd_dict[filename][ext]["Columns"], dir1, dir2)
            elif "ImageDiff" in fd_dict[filename][ext].keys():
                diff_images(filename, dir1, dir2) 
                

