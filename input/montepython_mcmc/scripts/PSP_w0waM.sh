#!/usr/local_rwth/bin/bash

# ask for four tasks (which are 4 MPI ranks)
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --account=rwth1304
##comment --account=rwth1304

# ask for 24 threads per task=MPI rank (which is 1 thread per core on one socket on CLAIX18)
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --output="PSP_w0waM.out"
#SBATCH --error="PSP_w0waM.err"
#SBATCH --job-name="PSP_w0waM"

#
#################
# ATTENTION !!! #
#################
# Divide the needed memory per task through the cpus-per-task, as slurm requests memory per cpu, not per task !
# Example:
# You need 24 GB memory per task, you have 24 cpus per task ordered
# order 24GB/24 -> 1G memory per cpu (i.e., per thread)

#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)

module purge
module restore intelpy
source ~/pyenv/intelpy/bin/activate

PYTHON=python
which python
echo /opt/MPI/bin/mpiexec" " -np 2
which $MPIEXEC
echo $MPIEXEC" "$FLAGS_MPI_BATCH
####   Before this, machine-specific commands

## Choose paths - usually unmodified if running inside Euclid_KP_nu
rootdir="../../Euclid_KP_nu/"
datadir="data/"
inputdir="input/montepython_mcmc/"
scriptsdir="input/montepython_mcmc/scripts/"
resultsdir="results/montepython_mcmc/"
source "${scriptsdir}mpfunctions.sh"

check_help $1
process_arguments "$@"

echo "****** Running script with $run ******"

cd ../MontePythons/montepython_kpnu/

## Choose case to run and extra options
CASE="PSP_w0waM"  # for changes that only involve options in the .param file
EXTRA="superpessimistic"       # for changes that only involve options in the .data file

## Print INPUT .param file
if [ "$EXTRA" = "superpessimistic" ]; then         
	### This is needed because the pessimistic cases are not in the filename of the .param
    INPUT="${rootdir}${inputdir}${CASE}.param"
else
    INPUT="${rootdir}${inputdir}${CASE}_${EXTRA}.param"
fi

echo "Running param file: $INPUT"

## Print results directory
CHAINS="${rootdir}${resultsdir}${CASE}_${EXTRA}"
echo "Chains will be saved to: $CHAINS "

# Proposal CovMat to be used
COVMAT="covmat/${CASE}_${EXTRA}.covmat"
#COVMAT="covmat/PSP_w0waMN.covmat"

usecovmat=true #    false or 0
rm_oldchains=true # false or 0
rm_fiducial=true #  false or 0
# MontePython chains default parameters
def_jumping="1.6"
def_Nsteps="100000"
def_upd="100"
def_superupd="20"
# MontePython parameter for info and plots
## uncomment to overwrite default or line argument
# non_markov=true

# Run functions according to options
check_run_print
check_run_info
check_run_fiducial
check_run_run
