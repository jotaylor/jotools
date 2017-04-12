#!/bin/bash

# get current working directory
cwd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# first and second input arguments
direc0=$1
branch=$2
# get rid of trailing slash if there is on
direc=${direc0%/}

# go to git directory and checkout branch
cd $direc
git checkout $branch
git pull

# install calcos
rm -rf build/
rm lib/calcos/version.py
export PYTHONPATH="${direc}/install/lib/python"
export lref="/grp/hst/cdbs/lref/"
python ./setup.py build
python ./setup.py install --home=${direc}/install

# Return to original directory
cd $cwd
