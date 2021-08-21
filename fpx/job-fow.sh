#!/bin/bash
#./job-fow.sh input/inputfile.in number_of_process ( ./job-fow.sh
#input/example.in 8 )

echo $$
dirname=fp$$
outfilename=fp$$.out
mkdir output/${dirname}
cd output/${dirname}
cp ../../$1 .

mpirun -n $2 ../../fp <../../$1 >${outfilename}
