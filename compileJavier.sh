#!/bin/bash

cd build

cmake -DNUMODIS_DIR=/home/gonzalez/WORK/Numodis_git/NUMODIS/ -DMEDCOUPLING_ROOT_DIR=/home/gonzalez/WORK/Numodis_git/NUMODIS/packages/medcoupling-install -DELMER_DIR=/home/gonzalez/TOOLS/elmer/install/ -DBOOSTHOME=/home/gonzalez/WORK/Numodis_git/NUMODIS/packages/boost_1_63_0 -DUSE_MKL=OFF -DLAPACKHOME=/home/gonzalez/WORK/Numodis_git/NUMODIS/packages/lapack ..
