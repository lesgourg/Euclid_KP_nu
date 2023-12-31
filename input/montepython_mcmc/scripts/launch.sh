#!/usr/local_rwth/bin/bash

export scriptsdir="input/montepython_mcmc/scripts"
export mpathopts="$scriptsdir/path_opts.sh"
export slurmopts="$scriptsdir/slurm_opts.sh"
export runopts="$scriptsdir/run_opts.sh"
export tempscript="$scriptsdir/tempjob.sh"
export tempopts="$scriptsdir/tempopts.sh"

source "$scriptsdir/mpfunctions.sh"
#
prefix="KP"
Case=$1

action=${2:-"nada2"}
extra=${3:-"nada3"}
allargs=$@
actargs="${@:2}"
echo $allargs
echo $Case
echo $action
echo $extra
echo $actargs

if [[ "$Case" = "showjobid" ]]; then
	JOBID=$action
	scontrol show job $JOBID
	exit 0
fi	

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
