#!/bin/bash
# - This script submits jobs to the SLURM batch system (used in TU Dresden computing cluster).
#   The script activates python environment, sets SLURM parameters.
# - The script will ask for N cores in the cluster, where N is the argument of the script
#   `./slurm_submit.sh <nthreads>`


d=$(date +%Y-%m-%d)  # Output file and job name are set to the current date.
job_name="$d" 
OUTFILE="$d.out"
source activate py37 # activate conda environment

# This is a workaround: usually no commands are allowed before #SBATCH.
# NB: Access `$SLURM_CPUS_PER_TASK` by escaping $-sign`\$SLURM_CPUS_PER_TASK`

sbatch << EOT
#!/bin/bash
#SBATCH --time=25:00:00   				# walltime
#SBATCH --nodes=1   					# number of nodes
#SBATCH --ntasks=1    					# limit to one node
#SBATCH --cpus-per-task=$1				# number of processor cores (i.e. threads)
#SBATCH --mem-per-cpu=300M   			# memory per CPU core
#SBATCH -J "$job_name"   				# job name
#SBATCH --mail-user=anton.solovev@tu-dresden.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A p_cilia						# project name

python master_script_escape.py $1  > "$OUTFILE"

exit 0
EOT