#!/usr/local_rwth/bin/zsh

# Function to display the help message
display_help() {
    echo "Usage: $0 [options] [run|fiducial|info|print] [more options]"
    echo
    echo "This script runs a MontePython job with MPI support."
    echo
    echo "Options:"
    echo "  -h, --help       Display this help message and exit."
    echo
    echo "Modes:"
    echo "  dryrun           Test (dry) Run the MontePython script printing commands to be run"
    echo "  run              Run the MontePython script with MPI support (default)."
    echo "  fiducial         Create the fiducial file and exit."
    echo "  info             Obtain information from existing chains and exit."
    echo "  print            Print configuration details and exit."
    echo
    echo "More options for info: "
    echo "  --second_chain   Path to chains to plot together with the original chains"
    echo
    echo "Default values:"
    echo "  run                  = 'run'"
    echo "  def_jumping          = 1.6"
    echo "  def_Nsteps           = 100000"
    echo "  def_upd              = 100"
    echo "  def_superupd         = 20"
    echo "  usecovmat            = true"
    echo "  rm_oldchains         = true"
    echo "  rm_fiducial          = true"
    echo
    echo "Example:"
    echo "  $0 run"
    exit 1
}

# Check if the user provided the -h flag
check_help() {
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    display_help
    exit
fi
}

export PYTHON=python

export date_time=$(date +"%Y-%m-%d_%H-%M")
echo "Current Date and Time: $date_time"

process_arguments() {
export run="${1:-"run"}"
# Initialize the option variable
export secondchain=""
# Loop through the command line arguments
for arg in "$@"
do
  case $arg in
    --second_chain=*)
      # Extract the option value after the equal sign
      export secondchain="${arg#*=}"
      # Remove the current argument from the positional arguments
      #
      echo $secondchain
      #shift
      ;;
    *)
      # Unknown option, you can handle it here or ignore it
      ;;
  esac
done
}

function check_run_print {
echo "Files used: "
echo "INPUT: $INPUT"
echo "CHAINS: $CHAINS"
echo "CASE: $CASE"
echo "EXTRA: $EXTRA"
echo "COVMAT: $COVMAT"
    if [ "$run" = "print" ]; then
        exit 1
    fi
}

function check_run_info {
    if [[ -d "$CHAINS" && "$run" = "info" ]]; then
        echo "CHAINS folder exists, obtaining info from chains..."
        if [[ -n "$secondchain" ]]; then echo "comparing against second chains at: $secondchain"; fi
        if [[ $non_markov == true ]]; then nmkvopt="--keep-non-markovian"; echo "non-markovian option selected, keeping non-markovian points"; else nmkvopt=""; fi
	chainsnum=${def_Nsteps}
	echo "Analyzing chains:  $CHAINS/$chainsnum "
	pwdir=$(pwd)
	echo $pwdir
        echo "$PYTHON montepython/MontePython.py info --want-covmat --plot-mean $nmkvopt "$CHAINS/\*$chainsnum\*" $secondchain"
        $PYTHON montepython/MontePython.py info --want-covmat --plot-mean $nmkvopt "$CHAINS/"*"$chainsnum"* $secondchain
	echo "Copyng covmat into montepython's covmat folder"
        cp -v "${CHAINS}/${CASE}_${EXTRA}.covmat" "covmat/"
        echo "Info obtained, exit script."
        exit 1
    fi
}

function extract_likes {
    if [[ -e $INPUT ]]; then
    exp_line=$(grep "data.experiments" $INPUT)
    else
    echo "$INPUT file not found, exiting script..."
    exit 0
    fi
    # Remove the prefix and suffix to isolate the comma-separated values
    exp_values=${exp_line#*[}
    exp_values=${exp_values%]*}
    # Split the values into an array
    IFS=',' read -ra EXPERIMENTS <<< "$exp_values"
    # Remove the leading and trailing single quotes from each element
    for i in "${!EXPERIMENTS[@]}"; do
        EXPERIMENTS[$i]=$(echo "${EXPERIMENTS[$i]}" | tr -d "'")
    done

    # You can now access any experiment using its index
    # Indexes start from 0 in bash

    # Print all experiments to check
    for index in "${!EXPERIMENTS[@]}"; do
        echo "Experiment likelihood: $(($index+1)): ${EXPERIMENTS[$index]}"
    done
    export EXPERIMENTS="$(printf '%s,' "${EXPERIMENTS[@]}")"
}

get_fiducial_file() {
  local likename=$1
  local DATA=$2
  local FIDUCIAL_SEARCH="$likename.fiducial_file"
  export FIDUCIAL_FILE=$(grep $FIDUCIAL_SEARCH $DATA | awk -F "$FIDUCIAL_SEARCH *= *" '{print $2}' | sed 's/"//g' | sed "s/'//g")
  #echo "Found: "
  #echo $FIDUCIAL_FILE
}

get_fiducial_files_array() {
#declare -a FIDUCIAL_FILES=()
IFS=',' read -r -a EXPERIMENTS <<< "$EXPERIMENTS"
for ii in "${!EXPERIMENTS[@]}"; do
  likename=${EXPERIMENTS[$ii]}
  echo "Experiment number: $(($ii+1)):  $likename"
  DATA="montepython/likelihoods/${likename}/${likename}.data"
  echo $DATA
  if [[ $likename == *"euclid"* ]]; then
    DATA_EXTRA="montepython/likelihoods/${likename}/${likename}.data.${EXTRA}"
    echo "Data file with extra options: $DATA_EXTRA"
    echo "Copying DATA_EXTRA to DATA"
    cp -v $DATA_EXTRA $DATA
    fi
  get_fiducial_file $likename $DATA
    echo "Fiducial file for $likename: $FIDUCIAL_FILE"
    FIDUCIAL_FILES[$ii]="${FIDUCIAL_FILE}"
    if [[ -e "${datadir}${FIDUCIAL_FILE}" && $rm_fiducial == true ]]; then 
       rm -rv "${datadir}${FIDUCIAL_FILE}"
    else 
       echo "${datadir}${FIDUCIAL_FILE} does not exist yet or the option to remove the existing fiducial is not wanted" 
    fi
  #fi
  export FIDUCIAL_FILES="$(printf '%s,' "${FIDUCIAL_FILES[@]}")"
  export EXPERIMENTS="$(printf '%s,' "${EXPERIMENTS[@]}")"
done
}


function remove_oldchains {
if [[ -d "$CHAINS" && $rm_oldchains == true ]]; then 
	echo "Removing $CHAINS folder after copying to a backup folder at ${CHAINS}_bckp_${date_time}"
	cp -vr $CHAINS "${CHAINS}_bckp_${date_time}"
	rm -Rfv $CHAINS
else 
	echo "Chains $CHAINS folder does not exist yet or option to remove old chains is not wanted"
fi
}

function check_run_fiducial {
echo "Checking for existing chains and fiducial files"
extract_likes
get_fiducial_files_array
remove_oldchains
if [[ -e $COVMAT && $usecovmat == true ]]; then Copt="-c $COVMAT"; else echo "Covmat $COVMAT non-existing or use covmat option not wanted"; Copt=""; fi
echo "Creating fiducial"
if [[ "$run" = "fiducial" || "$run" = "run" ]]; then
  $PYTHON montepython/MontePython.py run -p $INPUT -o $CHAINS -f 0 -N 1 $Copt
  echo "Generated fiducial"
  echo "Running chain with 1 point on fiducial"
  $PYTHON montepython/MontePython.py run -p $INPUT -o $CHAINS -f 0 -N 1 $Copt --display-each-chi2
  if [ "$run" = "fiducial" ]; then
   echo "Only fiducial run requested, exiting script."
   exit 1
  fi
fi
if [ "$run" = "dryrun" ]; then
  echo "$PYTHON montepython/MontePython.py run -p $INPUT -o $CHAINS -f 0 -N 1 $Copt"
  echo "Generated fiducial"
  echo "Running chain with 1 point on fiducial"
  echo "$PYTHON montepython/MontePython.py run -p $INPUT -o $CHAINS -f 0 -N 1 $Copt"
fi
}

function check_run_run {
if [ "$run" = "dryrun" ]; then
	echo "Launching chains in dry run:"
	echo "$MPIEXEC $FLAGS_MPI_BATCH $PYTHON montepython/MontePython.py run -o $CHAINS --conf myconf.conf -f "$def_jumping" -N "$def_Nsteps" --update "$def_upd" --superupdate "$def_superupd" $Copt"
  exit 1
fi
if [ "$run" = "run" ]; then
	echo "Launching chains"
	$MPIEXEC $FLAGS_MPI_BATCH $PYTHON montepython/MontePython.py run -o $CHAINS --conf myconf.conf -f "$def_jumping" -N "$def_Nsteps" --update "$def_upd" --superupdate "$def_superupd" $Copt
else
	echo "run input argument not specified, MP run with chains has not been started"
	exit 1
fi
}

