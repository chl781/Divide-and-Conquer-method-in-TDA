#!/bin/bash
set -e

tar -xzf R402.tar.gz
tar -xzf packages.tar.gz
tar -xzf functions.tar.gz


export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages

# run your script
Rscript merge.R $1 $2
