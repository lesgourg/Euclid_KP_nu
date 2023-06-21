#!/usr/local_rwth/bin/bash

# ask for four tasks (which are 4 MPI ranks)
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --account=rwth1304
##comment --account=rwth1304

# ask for 24 threads per task=MPI rank (which is 1 thread per core on one socket on CLAIX18)
#SBATCH --cpus-per-task=4
#SBATCH --time=23:59:00
#SBATCH --output="mpPSP9.out"
#SBATCH --error="mpPSP9.err"
#SBATCH --job-name="mpPSP9"

#
#################
# ATTENTION !!! #
#################
# Divide the needed memory per task through the cpus-per-task, as slurm requests memory per cpu, not per task !
# Example:
# You need 24 GB memory per task, you have 24 cpus per task ordered
# order 24GB/24 -> 1G memory per cpu (i.e., per thread)

#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)

#change these modules acording to how it was compiled! they must be identical
#modules in this order!
module purge
module load gompi
module load GSL/2.7
module load Python

#which mpicc

##so python knows we are in this custom enviroment with our own packages
##export INSTALLPATH=/work/qh043052/test/install
##export PATH=$INSTALLPATH/bin:/home/qh043052/anaconda3/bin:/home/qh043052/anaconda3/condabin:/opt/MPI/bin:/opt/intel/impi/2018.4.274/compilers_and_libraries/linux/mpi/bin64:/opt/intel/Compiler/19.0/1.144/rwthlnk/bin/intel64:/opt/slurm/bin:/usr/local_host/bin:/usr/local_host/sbin:/usr/local_rwth/bin:/usr/local_rwth/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/opt/singularity/bin:/usr/local/bin:/usr/local/sbin
##export PYTHONPATH=$INSTALLPATH/lib/python3.8/site-packages

#env > enviroment.log
#
#
#

source /home/qh043052/pyenv/montepython/bin/activate
PYTHON=python
echo /opt/MPI/bin/mpiexec" " -np 2
which $MPIEXEC
echo $MPIEXEC" "$FLAGS_MPI_BATCH
####   Before this, machine-specific commands

## Choose paths - usually unmodified if running inside Euclid_KP_nu
rootdir="../../Euclid_KP_nu/"
datadir="data/"
inputdir="input/montepython_mcmc/"
resultsdir="results/montepython_mcmc/"
source "${inputdir}mpfunctions.sh"

check_help $1
process_arguments "$@"

echo "****** Running script with $run ******"

cd ../MontePythons/montepython_kpnu/

## Choose case to run and extra options
CASE="PSP_w0waMN"  # for changes that only involve options in the .param file
EXTRA="pessimistic"       # for changes that only involve options in the .data file

## Print INPUT .param file
INPUT="${rootdir}${inputdir}${CASE}.param"
echo "Running param file: $INPUT"

## Print results directory
CHAINS="${rootdir}${resultsdir}${CASE}_${EXTRA}"
echo "Chains will be saved to: $CHAINS "

# Proposal CovMat to be used
COVMAT="covmat/${CASE}_${EXTRA}.covmat"
#COVMAT="covmat/PSP_MN.covmat"

usecovmat=true #    false or 0
rm_oldchains=true # false or 0
rm_fiducial=true #  false or 0
# MontePython chains default parameters
def_jumping="1.6"
def_Nsteps="100000"
def_upd="100"
def_superupd="20"
# MontePython parameter for info and plots
non_markov=true

# Run functions according to options
check_run_print
check_run_info
check_run_fiducial
check_run_run
