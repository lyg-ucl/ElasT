#!/bin/bash --login

# Select nodes 
#PBS -l select=1
#PBS -l walltime=00:19:00
#PBS -A your-budget-code 
# Move to directory that script was submitted from
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
# Set the number of threads to 1
export OMP_NUM_THREADS=1

echo "Running $jobname"

###############################################################################
# your vasp command
vasp="aprun -n 24 -N 24 /work/n03/n03/yunguo/vasp5.3.5/bin/vasp5"

# prepare structures and inputs
~/bin/elastic 

# run relaxation 
for i in $(ls -d */)
do
  cd ${i%%/}
  if [ ! -f OUTCAR ]; then
    $vasp 
  fi
  cd ..
done

# post process
~/bin/elastic >result
