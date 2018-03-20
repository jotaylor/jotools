#! /usr/bin/env python
'''
Run fitsdiff on two sets of data.

Use: 
    python fitsdiff.py -dir1 /path/to/dir1 -dir2 /path/to/dir2 
'''


import glob
import os
import subprocess
import multiprocessing as mp
import argparse
from functools import partial
from astropy.io import fits

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
        fd_output : list
            List of output from individual fitsdiff processes. 
    '''
    
    fd_func = partial(run_fd, fd_options)
    pool = mp.Pool(processes=10)
    fd_output = pool.map(fd_func, fd_files)
    return fd_output

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

def run_fd(fd_options, item):
    '''
    Run fitsdiff using astropy.

    Parameters:
    -----------
        fd_options : dict
            Dictionary containing command line options for fitsdiff.
        item : str
            Individual dataset to be used as a basis for comparison in 
            fitsdiff, from dir1.

    Returns:
    --------
        fd.report()
    '''
    
    basename = os.path.basename(item)
    file2 = os.path.join(fd_options["dir2"], basename)
    if not os.path.exists(file2):
        print("Expected {0}, but it does not exist!".format(file2))
        return ""

    try:
        print("Running on {0}".format(item))
        fd = fits.diff.FITSDiff(item, 
                                os.path.join(fd_options["dir2"], basename),
                                ignore_keywords = fd_options["ignore_keywords"],
                                ignore_comments = fd_options["ignore_comments"],
                                tolerance = fd_options["precision"])
    except Exception as e:
        print("Something went wrong with {}:\n{}\n".format(item,e))
        return ""

    return fd.report()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

def write_log(fd_output, outfile):
    '''
    Write output of fitsdiff to logfile.

    Parameters:
    -----------
        fd_output : list
            List of output from individual fitsdiff processes.
        outfile : str
            Name of file to write fitsdiff results to. 
    '''

    if os.path.exists(outfile):
        print("Overwriting {0}!".format(outfile))
        os.remove(outfile)
    else:
        print("Writing output to {0}".format(outfile))

    with open(outfile, "a") as f:
        for diff in fd_output:
            f.write(diff)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="prl", action="store_true", default=False,
                        help="Switch to parallelize fitsdiff process")
    parser.add_argument("--dir1", dest="dir1", help="First directory to compare")
    parser.add_argument("--dir2", dest="dir2", help="Second directory to compare")
    parser.add_argument("--precision", dest="precision", type=str,
                        default="2.e-7", help="Precision in fitsdiff")
    parser.add_argument("--comments", dest="ignore_comments", type=str,
                        default="gyromode,stdcfff,stdcffp,checksum,datasum",
                        help = "Comments to ignore in fitsdiff")
    parser.add_argument("--keywords", dest="ignore_keywords", type=str,
                        default="filename,date,iraf-tlm,cal_ver,proctime,checksum,datasum,opus_ver,ureka_v",
                        help="Keywords to ignore in fitsdiff") 
    parser.add_argument("-o", "--out", dest="outfile", default="fd_output.log",
                        help="Name of output fitsdiff log")
    args = parser.parse_args()
    
    all_files = glob.glob(os.path.join(args.dir1, "*fits"))
    all_exts = ["corrtag", "counts", "lampflash", "x1d", "x1dsum", 
                "x1dsum1", "x1dsum2", "x1dsum3", "x1dsum4", "flt"]
    
    fd_files = [x for x in all_files if os.path.basename(x).split("_")[1].split(".")[0] in all_exts]

    ignore_keywords = [x for x in args.ignore_keywords.split(",")]
    ignore_comments = [x for x in args.ignore_comments.split(",")]
    fd_options = {"dir2":args.dir2, "precision":float(args.precision), "ignore_comments":ignore_comments, "ignore_keywords":ignore_keywords} 
    
    if args.prl:
        fd_output = parallel(fd_files, fd_options)
    else:
        fd_output = ""
        for item in fd_files:
            fd_output += run_fd(fd_options, item)

    write_log(fd_output, args.outfile)
