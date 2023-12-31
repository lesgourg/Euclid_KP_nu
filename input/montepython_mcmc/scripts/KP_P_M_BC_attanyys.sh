#!/usr/local_rwth/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --account=rwth1304
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --output="slurm_outerr/KPattanyys_P_M_BC.out"
#SBATCH --error="slurm_outerr/KPattanyys_P_M_BC.err"
#SBATCH --job-name="KPattanyys_P_M_BC"
#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)
#change these modules acording to how it was compiled! they must be identical
module purge
module restore intelpy
source /home/qh043052/pyenv/intelpy/bin/activate
PYTHON=python
#echo `which mpicc`
echo "$MPIEXEC $FLAGS_MPI_BATCH"
# Specific cases to run as likelihoods
export CASE="P_M_BC"
# for changes that only involve options in the .data file
export EXTRA="superpessimistic"       

## Choose paths - usually unmodified if running inside Euclid_KP_nu
inputdir="input/montepython_mcmc/"
resultsdir="results/montepython_mcmc/"
export scriptsdir="input/montepython_mcmc/scripts"


source "$scriptsdir/mpfunctions.sh"

check_help $1
process_arguments "$@"

echo "****** Running script with $run ******"

# MontePython path relative to Euclid_KP_nu
montepythondir="../MontePythons/montepython_kpnu/"
# Paths realtive to MontePython
rootdir="../../Euclid_KP_nu/"
datadir="data/"
covdir="covmat/"


## Print INPUT .param file
if [ "$EXTRA" = "superpessimistic" ]; then         
	### This is needed because the pessimistic cases are not in the filename of the .param
    INPUT="${rootdir}${inputdir}${CASE}.param"
else
    INPUT="${rootdir}${inputdir}${CASE}_${EXTRA}.param"
fi
echo "Running param file: $INPUT"

## Print results directory
export CHAINS="${rootdir}${resultsdir}${CASE}_${EXTRA}"
echo "Chains will be saved to: $CHAINS "

# Proposal CovMat to be used
#export COVMAT="${covdir}${CASE}_${EXTRA}.covmat"
export COVMAT="covmat/P_M_superpessimistic.covmat"

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


#Move to the working directory
cd $montepythondir
# Run functions according to options

check_run_print
check_run_info
check_run_fiducial
check_run_run
