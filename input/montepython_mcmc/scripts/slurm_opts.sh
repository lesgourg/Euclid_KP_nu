#!/usr/local_rwth/bin/bash
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --account=rwth1304
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --output="slurm_outerr/$jname_$casevar.out"
#SBATCH --error="slurm_outerr/$jname_$casevar.err"
#SBATCH --job-name="$jname_$casevar"
#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)
#change these modules acording to how it was compiled! they must be identical
module purge
module restore intelpy
source /home/qh043052/pyenv/intelpy/bin/activate
PYTHON=python
#echo `which mpicc`
echo "$MPIEXEC $FLAGS_MPI_BATCH"
