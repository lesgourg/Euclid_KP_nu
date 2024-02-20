# Specific cases to run as likelihoods
export CASE="$casevar"
# for changes that only involve options in the .data file
#export EXTRA="superpessimistic"       
export EXTRA="optimisticBC"       

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
#if [ "$EXTRA" = "superpessimistic" ]; then         
	### This is needed because the pessimistic cases are not in the filename of the .param
#INPUT="${rootdir}${inputdir}${CASE}.param"
#else
#INPUT="${rootdir}${inputdir}${CASE}_${EXTRA}.param"
#fi

## .param file does not contain info about pess, opt, superpess. It is contained in .data
INPUT="${rootdir}${inputdir}${CASE}.param"
echo "Running param file: $INPUT"

## Print results directory
export CHAINS="${rootdir}${resultsdir}${CASE}_${EXTRA}"
echo "Chains will be saved to: $CHAINS "

# Proposal CovMat to be used
#export COVMAT="${covdir}${CASE}_${EXTRA}.covmat"
export COVMAT="covmat/P_M_BC_superpessimistic.covmat"

usecovmat=true #    false or 0
rm_oldchains=true # false or 0
rm_fiducial=true #  false or 0

# MontePython chains default parameters
def_jumping="1.6"
def_Nsteps="110000"
def_upd="100"
def_superupd="20"
# MontePython parameter for info and plots
## uncomment to overwrite default or line argument
# non_markov=true


