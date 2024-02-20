#!/usr/local_rwth/bin/bash


# Function to display help information
show_help() {
    echo "Usage: $0 [pathopts] [Case] [action] [extra]"
    echo "Usage: $0 --showjobid JOBID"
    echo
    echo "This script performs various actions based on the provided arguments."
    echo
    echo "Arguments:"
    echo "  pathopts : The script providing location of paths and options. Without the .sh extension."
    echo "  Case     : The MP case to process. Corresponding to the .param file to be run with MP, but leaving the '.param' extension out."
    echo "  action   : The action to perform. Options include 'writejob', 'fiducial', 'info', and 'run'."
    echo "  extra    : Additional arguments for the action."
    echo
    echo "Options:"
    echo "  -h, --h : Show this help message and exit."
    echo "  --showjobid : Show the job information corresponding to a JOBID number."
    echo
    echo "Examples:"
    echo "  $0 mypathopts.sh myparamfile run"
    echo
}

# Check for help option
if [[ "$1" == "-h" || "$1" == "--h" ]]; then
    show_help
    exit 0
fi

if [[ "$1" = "--showjobid" ]]; then
	JOBID=$2
	scontrol show job $JOBID
	exit 0
fi	


check_filedir() {
if [ ! -e "$1" ]; then
  echo "$1 does not exist."
  script=$(readlink -f $0)
  echo "exiting script $script"
  exit 1
fi
}

### Specify a file with proper paths and settings, opt, pess, covmats, etc
pathopts=$1

export scriptsdir="input/montepython_mcmc/scripts"
check_filedir $scriptsdir
export mpathopts="$scriptsdir/${pathopts}.sh"
echo "Checking if file in argument 1: $1 exists..."
check_filedir $mpathopts
export slurmopts="$scriptsdir/slurm_opts.sh"
check_filedir $slurmopts
export runopts="$scriptsdir/run_opts.sh"
check_filedir $runopts
export tempscript="$scriptsdir/tempjob.sh"
export tempopts="$scriptsdir/tempopts.sh"

export mpfunctions="$scriptsdir/mpfunctions.sh"
source $mpfunctions
#
prefix="KP"
Case=$2

action=${3:-"nada2"}
extra=${4:-"nada3"}
allargs=$@
actargs="${@:3}"
echo $allargs
echo $Case
echo $action
echo $extra
echo $actargs


hashjob $Case $prefix 


if [[ "$action" = "writejob" ]]; then
        cp -v $tempscript "${scriptsdir}/${jobscript}"
        rm -f $tempscript
	echo "---> Job written to ${scriptsdir}/${jobscript}"
        if [[ "$extra" = "cat" ]]; then
	$extra "${scriptsdir}/${jobscript}"
        fi
	exit 0
fi
if [[ "$action" = "fiducial" ]]; then 
        cp -v $tempscript "${scriptsdir}/${jobscript}"
        rm -f $tempscript
	export OMP_NUM_THREADS=8
	pwd
	echo "---> Running case ${CASE}_${EXTRA}"
        echo "---> Running locally ${scriptsdir}/${jobscript} $action"
	bash "${scriptsdir}/${jobscript}" fiducial
	exit 0
fi
if [[ "$action" = "info" ]]; then 
        cp -v $tempscript "${scriptsdir}/${jobscript}"
        rm -f $tempscript
	echo "---> Running case ${CASE}_${EXTRA}"
        echo "---> Running locally ${scriptsdir}/${jobscript} $actargs"
	bash "${scriptsdir}/${jobscript}" $actargs
	exit 0
fi
if [[ "$action" = "run" ]]; then
        cp -v $tempscript "${scriptsdir}/${jobscript}"
        rm -f $tempscript
	echo "---> Launching job with name $jobname"
	JOBI=$(sbatch "${scriptsdir}/${jobscript}" run)
	echo "$JOBI"
	JOBID=$(echo $JOBI | egrep -o '[0-9.]+')
	echo "JOBID: $JOBID"
	sleep 5
	scontrol show job ${JOBID}
	exit 0
else 
	echo "---> Running case ${CASE}_${EXTRA} on arguments: $allargs"
	echo "*** Action $allargs not supported, exiting script"
	exit 1
fi
