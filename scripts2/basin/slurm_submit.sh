#!/bin/bash
# Run: `./batch.sh <nthreads>`
# Normally, no commands are allowed before #SBATCH. In this file - workaround
# Access `$SLURM_CPUS_PER_TASK` by escaping $-sign`\$SLURM_CPUS_PER_TASK`
# output file and job name - by date
# args: number_of_cpus

d=$(date +%Y-%m-%d)
job_name="$d" 
OUTFILE="$d.out"
source activate py37 # activate conda environment

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


python master_basin_loc_backward.py $1  > "$OUTFILE"


exit 0
EOT