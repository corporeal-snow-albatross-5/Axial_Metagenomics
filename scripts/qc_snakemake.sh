Follow Sarah's instructions on how to set up the snakelike pipeline: https://github.com/shu251/qc-trim, though the
pipeline is pasted below, just in case: 

Set up working environment and directory

(1) Get this repository and build conda environment
Requires an installation of Anaconda. To get more familiar with using conda environments, see this post using conda
with R or read about them here. The conda environment to run this pipeline was build to run Snakemake.

1. clone this repo
git clone https://github.com/shu251/qc-trim.git

2. build conda environment to run snakemake
conda env create --name snake-qc --file envs/snake.yaml

3. Enter this conda environment
conda activate snake-qc

(2) Set up your working directory
To run this pipeline, we need to update the config.yaml file which will tell the Snakemake pipeline where our raw 
fastq files are and where we want the output files to be placed.
Use a text editor like 'nano' or 'tmux' to modify these lines of the config file:
	•	Move all fastq files to qc-trim/raw_data/directory or change the raw_data: line in the config file to the full 
path to your fastq files.
	•	Since I am working on an HPC I use a scratch space to write output files. Change this scratch line to where 
your scratch directory is
	•	Within my scratch directory, the output file will be placed in a new directory of this name. In the example 
below, all my output files will be placed in the directory /vortexfs1/scratch/sarahhu/test-qc/. This way my scratch
directory will be organized as I run multiple projects.
	•	Also change the location of where adapters for your study can be found for the trimmomatic step

1. Change to location of all raw .fastq (or .fastq.gz) files

raw_data: /vortexfs1/omics/huber/shu/qc-trim/raw_data

2. Change to output directory (on HPC, use your scratch)
scratch:  /vortexfs1/scratch/sarahhu

3. Add a project name, output files will be placed in a directory of this name in scratch
proj_name: test-qc

4. Change to output directory to eventually move all trimmed reads to
outputDIR: /vortexfs1/omics/huber/shu/qc-trim

5. Location of illumina adapters, full path
adapters: /vortexfs1/omics/huber/db/adapters/illumina-adapters.fa

(3) Execute dryrun of snakemake pipeline

snakemake -np

# The command 'snakemake' automatically looks for 'Snakefile'

(4) Execute snakemake

snakemake --use-conda

# The '--use-conda' flag is necessary because we enable fastqc, trimmomatic, and multiqc conda environments

(5) Run snakemake with SLURM

Read about executing snakemake on a cluster and another explanation on how to execute with a submit script can be 
found here. Review the submit scripts available in submitscripts. Files in this directory include another config 
file called cluster.yaml, and two scripts to submit your snakemake pipeline to the cluster with and without the dry
run option. First, open and modify the cluster.yaml to fit your machine. Then test run using the provided submit
scripts.

# Make sure you are still in the snake-tagseq conda environment

bash submitscripts/submit-slurm-dry.sh

Outputs should all be green and will report how many jobs and rules snakemake plans to run. This will allow you to 
evaluate any error and syntax messages.
Once ready, use the submit-slurm.sh script to submit the jobs to slurm (see below). Run with screen, tmux, or
nohup.

bash submitscripts/submit-slurm.sh

# This will print jobs submitted as it runs. Monitor these jobs using ```squeue```

#example submit-slurm.sh script below:

#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=QC_HS_April2014       # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=2                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=640000                     # Job memory request
#SBATCH --time=24:00:00               	 # Time limit hrs:min:sec
#SBATCH --output=QC_HS_April2014.log     # Standard output/error
#export OMP_NUM_THREADS=8

module load anaconda/5.1

conda activate snake-qc

snakemake --use-conda

Some notes: with downstream processing (i.e. assembling of metagenomes, removing rRNA from 
metatranscriptome/RNA-SIP data), you will want to use your PAIRED TRIMMED files. They will be in fastq.gz format. 
Most programs (IDBA_UD and SortMeRNA) require paired end, trimmed reads to be merged and turned into a .fasta file.
Below (in the ‘Assembling Metagenomes with IDBA_UD’ section) is how to do this.
