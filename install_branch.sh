#!/bin/bash

# Get current working directory
cwd=$(pwd)

# Assign first and second input arguments
direc0=$1
branch=$2
# Get rid of trailing slash if there is one
direc=${direc0%/}

# Go to git directory and checkout branch
cd $direc
git fetch
git checkout $branch
git pull

# Install calcos
rm -rf build/
rm -f lib/calcos/version.py
export PYTHONPATH="${direc}/install/lib/python"
export lref="/grp/hst/cdbs/lref/"
python ./setup.py build
python ./setup.py install --home=${direc}/install

# Return to original directory
cd $cwd
