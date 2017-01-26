#! /usr/bin/env python
'''
Run fitsdiff on two sets of data.

Use: 
    python fitsdiff.py -dir /path/to/dir1 -dir2 /path/to/dir2 -p 
'''


import glob
import os
import subprocess
import multiprocessing as mp
import argparse
from functools import partial

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

def parallel(fd_files, fd_options):
    '''
    Parallelize the fitsdiff function
    
    Parameters:
    -----------
        fd_files : list
            List of files that will serve as the basis for comparison. These
            are from dir1.
        fd_options : dict
            Dictionary containing command line options for fitsdiff.
            
    Returns:
    --------
        Nothing
    '''
    
    fd_func = partial(run_fd, fd_options)
    pool = mp.Pool(processes=10)
    pool.map(fd_func, fd_files)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

def run_fd(fd_options, item):
    '''
    Run the fitsdiff command-line command using subprocess.

    Parameters:
    -----------
        fd_options : dict
            Dictionary containing command line options for fitsdiff.
        item : str
            Individual dataset to be used as a basis for comparison in 
            fitsdiff, from dir1.

    Returns:
    --------
        Nothing
    '''
    
    basename = os.path.basename(item)
    logname = "{0}.log".format(basename)
    with open(logname, "w") as log:
        try:
            subprocess.call("fitsdiff -c {0} -k {1} -d {2} {3} {4}".format(fd_options["ignore_comments"], fd_options["ignore_keywords"], fd_options["precision"], item, os.path.join(fd_options["dir2"], basename)), shell=True, stdout=log)
        except Exception as e:
            print("Something went wrong with {}:\n{}\n".format(item,e))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

def concat_files(fd_files):
    '''
    Parameters:
    -----------
    fd_files : list
        List of files that will serve as the basis for comparison. These
        are from dir1.
    
    Returns:
    --------
        Nothing
    '''
    lognames = [os.path.basename(x)+".log" for x in fd_files]
    with open("concat_fd.log", "w") as outfile:
        for logname in lognames:
            with open(logname) as infile:
                outfile.write(infile.read())

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="prl", action="store_true", default=False,
                        help="Switch to parallelize fitsdiff process")
    parser.add_argument("-dir1", dest="dir1", help="First directory to compare")
    parser.add_argument("-dir2", dest="dir2", help="Second directory to compare")
    parser.add_argument("-precision", dest="precision", type=str,
                        default="2.e-7", help="Precision in fitsdiff")
    parser.add_argument("-comments", dest="ignore_comments", type=str,
                        default="gyromode,stdcfff,stdcffp,checksum,datasum",
                        help = "Comments to ignore in fitsdiff")
    parser.add_argument("-keywords", dest="ignore_keywords", type=str,
                        default="filename,date,iraf-tlm,cal_ver,proctime,checksum,datasum,opus_ver,ureka_v",
                        help="Keywords to ignore in fitsdiff") 
    args = parser.parse_args()
    
#    all_files = glob.glob(os.path.join(args.dir1, "*fits"))
    all_exts = ["corrtag", "counts", "lampflash", "x1d", "x1dsum", 
                "x1dsum1", "x1dsum2", "x1dsum3", "x1dsum4", "flt"]
    
    fd_files = [x for x in all_files if os.path.basename(x).split("_")[1].split(".")[0] in all_exts]

    fd_files = all_files
        
    fd_options = {"dir2":args.dir2, "precision":args.precision, "ignore_comments":args.ignore_comments, "ignore_keywords":args.ignore_keywords} 
    
    if args.prl:
        parallel(fd_files, fd_options)
    else:
        for item in fd_files:
            run_fd(fd_options, item)

    concat_files(fd_files)

