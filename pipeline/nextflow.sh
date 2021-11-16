#! /bin/bash
#SBATCH -o /path/to/output/messages/nextflow%j.txt
#SBATCH -e /path/to/error/messages/nextflow-er%j.txt

set -e

// must have a conda virtual environment for dependencies
// create one and change this path to locate it
. /data/SGS/conda/etc/profile.d/conda.sh
conda activate

// load and run nextflow
// use -resume flag to pick up where left off if pipeline fails
module load nextflow
#nextflow run main.nf
nextflow main.nf -resume
