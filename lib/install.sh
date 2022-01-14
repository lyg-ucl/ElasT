#!/bin/bash 

tar -xvf blas-3.10.0.tgz
cd BLAS-3.10.0
make
mv blas_LINUX.a libblas.a
cd ..

tar -xvf lapack-3.10.0.tar.gz
cd lapack-3.10.0
cp make.inc.example make.inc
make
cd .. 
