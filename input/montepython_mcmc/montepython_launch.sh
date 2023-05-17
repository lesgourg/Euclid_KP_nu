#source /home/qh043052/pyenv/montepython/bin/activate
#PYTHON=python
#echo /opt/MPI/bin/mpiexec" " -np 2
#which $MPIEXEC
#echo $MPIEXEC" "$FLAGS_MPI_BATCH
####   Before this, machine-specific commands

## Choose paths - usually unmodified if running inside Euclid_KP_nu
rootdir="../Euclid_KP_nu/"
datadir="data/"
inputdir="input/montepython_mcmc/"
resultsdir="results/montepython_mcmc/"
source "${inputdir}mpfunctions.sh"

check_help $1
process_arguments "$@"

echo "****** Running script with $run ******"

cd ../montepython/

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
