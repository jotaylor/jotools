#! /user/bin/env python

'''
AUTHOR: Jo Taylor, RIA
DATE: 10/20/2015
NAME: make_regr_csv.py
DESCRIPTION:
Create a spreadsheet with information about the files in the COS regression
suite

''' 

import glob
from astropy.io import fits as pf
import os
import csv
import argparse

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def write_csv(data_dir, outfile, header_keys):
    csv_file = open(outfile, "w")
    writer = csv.writer(csv_file, dialect="excel", quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header_keys[0]+header_keys[1])
    uniq_filenames = []
    exts = list(header_keys.keys())
    exts.sort()
     
    for item in glob.glob(os.path.join(data_dir, "*raw*")):
        csv_line = []
        with pf.open(item) as hdulist:
            hdrs = [hdulist[0].header, hdulist[1].header]
        rootname = hdrs[0]["ROOTNAME"]
        if rootname in uniq_filenames:
            continue
        for ext in exts:
            for keyword in header_keys[ext]:
                try:
                    csv_line.append(hdrs[ext][keyword])
                except KeyError:
                    csv_line.append("N/A")
        writer.writerow(csv_line)
        uniq_filenames.append(rootname)
    
    csv_file.close()
    
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="data_dir", help="Path to regression data")
    parser.add_argument("-o", dest="outfile", help="Name of output CSV file")
    args = parser.parse_args()

    header_keys = {0: ["ROOTNAME", "OPT_ELEM", "CENWAVE", "APERTURE", "EXPTYPE", "OBSTYPE", "OBSMODE", "PROGRAM", "LIFE_ADJ"], 1: ["DATE-OBS"]}
    write_csv(args.data_dir, args.outfile, header_keys)
