#!/bin/bash --login
#PBS -N 6

# Select 128 nodes (maximum of 3072 cores)
#PBS -l select=4
#PBS -l walltime=00:59:00

# Replace this with your budget code
#PBS -A n03-lv-users 

#module add vasp5/5.4.1

# Move to directory that script was submitted from
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR

# Set the number of threads to 1
#   This prevents any system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

echo "Running $jobname"

for i in $(ls -d */)               # for each directory
do                                 
  cd ${i%%/}                       # enter dir
  if [ ! -f OUTCAR ]; then
    aprun -n 96 -N 24 /work/n03/n03/yunguo/vasp5.3.5/bin/vasp5
  fi
  cd ..
done
