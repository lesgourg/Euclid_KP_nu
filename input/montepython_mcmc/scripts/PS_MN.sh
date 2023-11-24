#!/usr/local_rwth/bin/bash

scriptsdir="input/montepython_mcmc/scripts/"
mpathopts="$scriptsdir/path_opts.sh"
slurmopts="$scriptsdir/slurm_opts.sh"
jobscript="$scriptsdir/tempjob.sh"

prefix="Kn_"

Case=$1
action=$2
allargs=$@

if [[ "$Case" = "showjobid" ]]; then
	JOBID=$action
	scontrol show job $JOBID
	exit 0
fi	

sed "s/\$jname/${jobname}/g" $slurmopts >> $job_script
sed "s/\$casevar/${Case}/g" $mpathopts >> $job_script

tempjob=$(cat $jobscript)

shash=`python -c "import hashlib; s='$tempjob'; print(int(hashlib.sha256(s.encode('utf-8')).hexdigest(), 16) % 10**6)"`
tash=`echo "$shash" | tr 0123456789 santiegouy`
jobhash=$tash
echo "Hashed jobname of script:  $jobhash"
jobname="${prefix}${jobhash}"
jobscript="${Case}_${jobhash}.sh"

cp -v $tempjob $jobscript

if [[ "$action" = "writejob" ]]; then
	echo "Job written to $jobscript"
	cat $jobscript
fi
if [[ "$action" = "run" ]]; then
	echo "Launching job with name $jobname"
	JOBID=$(sbatch $jobscript run)
	echo "JOBID: $JOBID"
	sleep 5
	scontrol show job ${JOBID}
fi
if [[ "$action" = "fiducial" ]]; then 
        source $jobscript
        echo "Running locally $jobscript"
	echo "Running locally $CASE_$EXTRA $action"
	bash $jobscript fiducial
	exit 0
else 
        source $jobscript
        echo "Running locally $jobscript"
	echo "Running locally $CASE_$EXTRA on arguments: $allargs"
	bash $jobscript $allargs
	exit 0
fi
