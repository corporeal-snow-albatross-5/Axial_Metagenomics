# **Axial RNA-SIP Sample Processing: From Raw Reads to Mapping mRNA to Concatenated Metagenome**

#### All scripts are located in the 'scripts' folder in this repo. Script name is under each step in the overview. 

### Overview:  
1. Trim Samples with Trimmomatic and QC using FastQC and MultiQC using Sarah Hu's Snakemake pipeline  
-qc_snakemake.sh  
2. Combine all the forward and all the reverse assembly segments using cat command then use fq2fa from IDBA_UD to convert .fastq output of trimming to fasta and merge forward and reverse reads  
-combine_fq2fa.sh  
3. Assemble each metagenome using IDBA_UD  
-idba_assembly.sh  
4. Run assembled metagenomes through MetaQUAST to assess assembly quality  
-metaquast.sh  
5. Annotate each metagenome with the Joint Genome Institute's IMG pipeline   
7. Getting annotated metagenomic files off of IMG  
8. Converting annotated metagenomes to .gff files with python before concatenating them  
gff_conversion.sh  
9. Concatenating .gff files to make one gigantic, metagenome ORF file  
10. Mapping mRNA reads back to concatenated metaG file using Kallisto  
-kallisto.sh (in scripts/kallisto_scripts_and_logs folder)  

### Documentation for each program:
1. Snakemake - https://snakemake.readthedocs.io/en/stable/  
- https://github.com/shu251/qc-trim  
2. Trimmomatic - http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf  
3. FastQC - https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf  
4. MultiQC - https://multiqc.info/docs/  
5. fq2fa - https://denbi-metagenomics-workshop.readthedocs.io/en/latest/assembly/idba_ud.html
- https://github.com/loneknightpy/idba
6. IDBA_UD -  https://denbi-metagenomics-workshop.readthedocs.io/en/latest/assembly/idba_ud.html
- https://github.com/loneknightpy/idba  
7. MetaQUAST - http://quast.sourceforge.net/docs/manual.html
8. JGI-IMG - https://img.jgi.doe.gov/submit/doc/IMGSubmissionUserGuide.pdf
9. Kallisto - https://pachterlab.github.io/kallisto/manual

### Raw files - all located in the sample_list file in this repo. Also see the output tree if you'd like to see all of the samples associated with Axial. The metadata folder also has some useful info about the samples.

#### **Note: all of these programs were run on WHOI's cluster, Poseidon, so the code for that is included. Additionally, all programs were installed using Conda Environments** 

## **1. Trim Samples with Trimmomatic and QC using FastQC and MultiQC using Sarah Hu's Snakemake pipeline**
```
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

```

## **2. Combine all the forward and all the reverse assembly segments using cat command then use fq2fa from IDBA_UD to convert .fastq output of trimming to fasta and merge forward and reverse reads**
```
#This was done on srun:
srun -p scavenger --time=08:00:00 --ntasks-per-node 2 --mem=40gb --pty bash 

#FS891

cat ERR694110_1_trim.fastq.gz ERR694111_1_trim.fastq.gz ERR694112_1_trim.fastq.gz ERR694113_1_trim.fastq.gz 
ERR694114_1_trim.fastq.gz ERR694115_1_trim.fastq.gz ERR694116_1_trim.fastq.gz ERR694117_1_trim.fastq.gz 
ERR694118_1_trim.fastq.gz ERR694119_1_trim.fastq.gz ERR694120_1_trim.fastq.gz ERR694121_1_trim.fastq.gz 
ERR694122_1_trim.fastq.gz ERR694123_1_trim.fastq.gz ERR694124_1_trim.fastq.gz ERR694125_1_trim.fastq.gz 
ERR694126_1_trim.fastq.gz > FS891_merged_cat_forward_reads.fastq.gz

cat ERR694110_2_trim.fastq.gz ERR694111_2_trim.fastq.gz ERR694112_2_trim.fastq.gz ERR694113_2_trim.fastq.gz 
ERR694114_2_trim.fastq.gz ERR694115_2_trim.fastq.gz ERR694116_2_trim.fastq.gz ERR694117_2_trim.fastq.gz 
ERR694118_2_trim.fastq.gz ERR694119_2_trim.fastq.gz ERR694120_2_trim.fastq.gz ERR694121_2_trim.fastq.gz 
ERR694122_2_trim.fastq.gz ERR694123_2_trim.fastq.gz ERR694124_2_trim.fastq.gz ERR694125_2_trim.fastq.gz 
ERR694126_2_trim.fastq.gz > FS891_merged_cat_reverse_reads.fastq.gz

fq2fa --merge FS891_merged_cat_forward_reads.fastq FS891_merged_cat_reverse_reads.fastq FS891_total_combined_reads.fasta


#FS903 

cat ERR694207_1_trim.fastq.gz ERR694208_1_trim.fastq.gz ERR694209_1_trim.fastq.gz ERR694210_1_trim.fastq.gz 
ERR694211_1_trim.fastq.gz ERR694212_1_trim.fastq.gz ERR694213_1_trim.fastq.gz > FS903_cat_forward_reads.fastq.gz

cat ERR694207_2_trim.fastq.gz ERR694208_2_trim.fastq.gz ERR694209_2_trim.fastq.gz ERR694210_2_trim.fastq.gz 
ERR694211_2_trim.fastq.gz ERR694212_2_trim.fastq.gz ERR694213_2_trim.fastq.gz > FS903_cat_reverse_reads.fastq.gz

gunzip FS903_cat_forward_reads.fastq.gz

gunzip FS903_cat_reverse_reads.fastq.gz

fq2fa ——merge FS903_cat_forward_reads.fastq FS903_cat_reverse_reads.fastq FS903_total_combined_reads.fasta


#FS904

cat ERR694199_1_trim.fastq.gz ERR694200_1_trim.fastq.gz ERR694201_1_trim.fastq.gz ERR694202_1_trim.fastq.gz 
ERR694203_1_trim.fastq.gz ERR694204_1_trim.fastq.gz ERR694205_1_trim.fastq.gz ERR694206_1_trim.fastq.gz > 
  FS904_cat_forward_reads.fastq.gz

cat ERR694199_2_trim.fastq.gz ERR694200_2_trim.fastq.gz ERR694201_2_trim.fastq.gz ERR694202_2_trim.fastq.gz 
ERR694203_2_trim.fastq.gz ERR694204_2_trim.fastq.gz ERR694205_2_trim.fastq.gz ERR694206_2_trim.fastq.gz > 
  FS904_cat_reverse_reads.fastq.gz
  
gunzip FS904_cat_forward_reads.fastq.gz 

gunzip FS904_cat_reverse_reads.fastq.gz 

fq2fa ——merge FS904_cat_forward_reads.fastq FS904_cat_reverse_reads.fastq FS904_total_combined_reads.fasta


#FS906

cat ERR1163082_1_trim.fastq.gz ERR1163083_1_trim.fastq.gz ERR1163084_1_trim.fastq.gz ERR1163085_1_trim.fastq.gz 
ERR1163086_1_trim.fastq.gz ERR1163087_1_trim.fastq.gz ERR1163089_1_trim.fastq.gz ERR1163090_1_trim.fastq.gz 
ERR1163091_1_trim.fastq.gz ERR1163092_1_trim.fastq.gz ERR1163093_1_trim.fastq.gz ERR1163094_1_trim.fastq.gz > 
FS906_cat_forward_reads.fastq.gz

cat ERR1163082_2_trim.fastq.gz ERR1163083_2_trim.fastq.gz ERR1163084_2_trim.fastq.gz ERR1163085_2_trim.fastq.gz 
ERR1163086_2_trim.fastq.gz ERR1163087_2_trim.fastq.gz ERR1163089_2_trim.fastq.gz ERR1163090_2_trim.fastq.gz 
ERR1163091_2_trim.fastq.gz ERR1163092_2_trim.fastq.gz ERR1163093_2_trim.fastq.gz ERR1163094_2_trim.fastq.gz > 
FS906_cat_reverse_reads.fastq.gz

gunzip FS906_cat_forward_reads.fastq.gz

gunzip FS906_cat_reverse_reads.fastq.gz

fq2fa ——merge FS906_cat_forward_reads.fastq FS906_cat_reverse_reads.fastq FS906_total_combined_reads.fasta


#FS907 

cat ERR1163105_1_trim.fastq.gz ERR1163106_1_trim.fastq.gz ERR1163107_1_trim.fastq.gz ERR1163108_1_trim.fastq.gz 
ERR1163109_1_trim.fastq.gz ERR1163110_1_trim.fastq.gz ERR1163111_1_trim.fastq.gz ERR1163112_1_trim.fastq.gz 
ERR1163113_1_trim.fastq.gz ERR1163114_1_trim.fastq.gz ERR1163115_1_trim.fastq.gz ERR1163116_1_trim.fastq.gz 
ERR1163117_1_trim.fastq.gz > FS907_cat_forward_reads.fastq.gz

cat ERR1163105_2_trim.fastq.gz ERR1163106_2_trim.fastq.gz ERR1163107_2_trim.fastq.gz ERR1163108_2_trim.fastq.gz 
ERR1163109_2_trim.fastq.gz ERR1163110_2_trim.fastq.gz ERR1163111_2_trim.fastq.gz ERR1163112_2_trim.fastq.gz 
ERR1163113_2_trim.fastq.gz ERR1163114_2_trim.fastq.gz ERR1163115_2_trim.fastq.gz ERR1163116_2_trim.fastq.gz 
ERR1163117_2_trim.fastq.gz > FS907_cat_reverse_reads.fastq.gz

gunzip FS907_cat_forward_reads.fastq.gz

gunzip FS907_cat_reverse_reads.fastq.gz

fq2fa ——merge FS907_cat_forward_reads.fastq FS907_cat_reverse_reads.fastq FS907_total_combined_reads.fasta


#FS908

cat ERR1163125_1_trim.fastq.gz ERR1163126_1_trim.fastq.gz ERR1163127_1_trim.fastq.gz ERR1163128_1_trim.fastq.gz 
ERR1163129_1_trim.fastq.gz ERR1163130_1_trim.fastq.gz ERR1163131_1_trim.fastq.gz ERR1163132_1_trim.fastq.gz 
ERR1163133_1_trim.fastq.gz ERR1163134_1_trim.fastq.gz ERR1163135_1_trim.fastq.gz > FS908_cat_forward_reads.fastq.gz

cat ERR1163125_2_trim.fastq.gz ERR1163126_2_trim.fastq.gz ERR1163127_2_trim.fastq.gz ERR1163128_2_trim.fastq.gz 
ERR1163129_2_trim.fastq.gz ERR1163130_2_trim.fastq.gz ERR1163131_2_trim.fastq.gz ERR1163132_2_trim.fastq.gz 
ERR1163133_2_trim.fastq.gz ERR1163134_2_trim.fastq.gz ERR1163135_2_trim.fastq.gz > FS908_cat_reverse_reads.fastq.gz

gunzip FS908_cat_forward_reads.fastq.gz

gunzip FS908_cat_reverse_reads.fastq.gz

fq2fa ——merge FS908_cat_forward_reads.fastq FS908_cat_reverse_reads.fastq FS908_total_combined_reads.fasta


#FS914
fq2fa ——merge ERR2021505_1_trim.fastq.gz ERR2021505_2_trim.fastq.gz ERR2021505_merged.fa

#FS915
fq2fa ——merge ERR2021507_1_trim.fastq.gz ERR2021507_2_trim.fastq.gz ERR2021507_merged.fa

#FS917
fq2fa ——merge ERR2021503_1_trim.fastq.gz ERR2021503_2_trim.fastq.gz ERR2021503_merged.fa

#Background_1500m_2015_DNA
fq2fa ——merge ERR2021513_1_trim.fastq.gz ERR2021513_2_trim.fastq.gz ERR2021513_merged.fa

#AnemomePlume_2015_DNA
fq2fa ——merge ERR2021511_1_trim.fastq.gz ERR2021511_2_trim.fastq.gz ERR2021511_merged.fa
```

#### **Note: AnemomePlume_2015_DNA, FS914_Anemone_DNA, Background_1500m_2015_DNA, FS917_Marker33_DNA, and FS915_Marker113_DNA only had one forward and reverse read pair associated with each (i.e. they were not split into many different read segments by the SRA), so they were taken directly through IDBA_UD after merging forward and reverse reads with fq2fa above**  


## **3. Assemble each metagenome using IDBA_UD**
```
Assembling Metagenomes with IDBA_UD

#Create environment, so that you can have all of the programs necessary for assembly in one place: 

conda create --name idba-ud_assembly_env 

#Enter environment:

conda activate idba-ud_assembly_env

#If you need to check available environments, type:

conda info --envs

#Load IDBA-UD into conda -- google "conda install idba ud" -- and you will receive this line of code:

conda install -c bioconda idba

#Files must be in one read file (so you have to merge your reads) and they must be in .fasta format. IDBA-UD has a program 
#(called fq2fa) built into its package, where you can easily do this. However, you will need to clone the git repo by doing
#this:

git clone https://github.com/loneknightpy/idba.git

#gunzip your files first to run through merging program:

for files in <type your file path here>; do; gunzip ${files}; done

#ex. for gunzip: for files in /vortexfs1/scratch/selkassas/FS903_Marker113_DNA/test-qc/trimmed/merged_files; do; gunzip 
#${files}; done

#also shown above in the cat and fq2fa code. 

#Slurm example below:

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_gunzip		           # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=2                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_gunzip%j.log               # Standard output/error

export OMP_NUM_THREADS=8

cd /vortexfs1/scratch/selkassas/qc-trim/raw_data/paired_trimmed_files_for_sortmerna

for file in *fastq.gz; do gunzip ${file}; done

#The command to merge files is: 

fq2fa --merge read_1_trim.fq read_2_trim.fq read_merged.fa 

#This is if you want to run in serial 
#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=<enter_job_name>      # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=160000                     # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=<enter_job_name>.log    # Standard output/error
#export OMP_NUM_THREADS=8

#This is if you want to run in parallel (pick either parallel or serial; not both)
#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_<enter_job_name>         # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=4                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_<enter_job_name>%j.log     # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length
#Here is the format your script should be in:

idba_ud -r <path to your merged, trimmed .fasta file (or if you already put in your slurm script to enter a
particular directory with all of your merged files, then that'll work too)> -o #(your output) <path to where you
want your output file to go/your assembly file name> #Then specify max kmer length or include any other flags.
example: --maxk 120 <-- maximum kmer length of 120. 

#example
#FS915_Marker113_DNA -- 1 metagenome - 2015 - read length of 151bp, so chose a max kmer length of 120. 
	#Sample ERR2021507
idba_ud -r /vortexfs1/scratch/selkassas/FS915_Marker113_DNA/test-qc/trimmed/merged_files/ERR2021507_merged.fa -o
/vortexfs1/scratch/selkassas/FS915_Marker113_DNA/ERR2021507_assembly --maxk 120

#To submit to slurm: 
sbatch axial_assemblies.sh (this file, but written in nano on Poseidon) 
#You will get an email that tells you that your job has started, along with a job id number. MAKE SURE YOU LOG YOUR
JOB-ID NUMBER! 
#check that your job has started via squeue
#check YOUR jobs only: squeue -u selkassas

-----------------------

#Full script that I used for each: 

##Assemblies Code

#FS891 

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS891                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS891%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS891_DNA/test-qc/trimmed/merged_files/FS891_total_combined_reads.fasta -o 
/vortexfs1/scratch/selkassas/FS891/idbaud_assemblies/FS891_DNA_combined_assembly --maxk 120

-----------------------

# FS903

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS903                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS903%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS903_Marker113_DNA/test-qc/trimmed/merged_files/FS903_Marker113_DNA_combin
ed.fa -o /vortexfs1/scratch/selkassas/FS903_Marker113_DNA/idbaud_assemblies/FS903_Marker113_DNA_combined_assembly 
--maxk 120

-----------------------

#FS904

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS904                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=4                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=1000000                                # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS904%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS904_Marker33_DNA/test-qc/trimmed/merged_files/FS904_Marker33_DNA_combined
.fa -o /vortexfs1/scratch/selkassas/FS904_Marker33_DNA/idbaud_assemblies/FS904_Marker33_DNA_combined_assembly 
--maxk 120

-----------------------

#FS906

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS906                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS906%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS906_Marker113_DNA/test-qc/trimmed/merged_files/FS906_Marker113_DNA_combin
ed.fa -o /vortexfs1/scratch/selkassas/FS906_Marker113_DNA/idbaud_assemblies/FS906_Marker113_DNA_combined --maxk 120

-----------------------

#FS907

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS907                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS907%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS907_Anemone_DNA/test-qc/trimmed/merged_files/FS907_Anemone_DNA_combined.f
a -o /vortexfs1/scratch/selkassas/FS907_Anemone_DNA/idbaud_assemblies/FS907_Anemone_DNA_combined_assembly --maxk 
120

-----------------------

#FS908

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS908                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS908%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS908_Marker33_DNA/test-qc/trimmed/merged_files/FS908_Marker33_DNA_combined
.fa -o /vortexfs1/scratch/selkassas/FS908_Marker33_DNA/idbaud_assemblies/FS908_Marker33_DNA_combined_assembly 
--maxk 120

-----------------------

#FS915

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS915                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS915%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS915_Marker113_DNA/test-qc/trimmed/merged_files/ERR2021507_merged.fa -o 
/vortexfs1/scratch/selkassas/FS915_Marker113_DNA/idbaud_assemblies/ERR2021507_assembly/assembly_maxkmer_120 maxk 
--120

-----------------------

#FS917

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS917                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS917%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS917_Marker33_DNA/test-qc/trimmed/merged_files/ERR2021503_merged.fa -o 
/vortexfs1/scratch/selkassas/FS917_Marker33_DNA/idbaud_assemblies/ERR2021503_assembly/assembly_maxkmer_120 --maxk 
120

-----------------------

#FS914

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_FS914                    # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_FS914%j.log                # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/FS914_Anemone_DNA/test-qc/trimmed/merged_files/ERR2021505_merged.fa -o 
/vortexfs1/scratch/selkassas/FS914_Anemone_DNA/idbaud_assemblies/ERR2021505_assembly/assembly_maxkmer_120 --maxk 
120

-----------------------  

#Background_1500m_2015_DNA

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_background               # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_background%j.log           # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/Background_1500m_2015_DNA/test-qc/trimmed/ERR2021513_merged.fa -o /vortexfs1/scratch/selkassas/Background_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_120/ --maxk 120

----------------------- 

#AnemonePlume_2015_DNA

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=parallel_anemone_plume            # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=36                           # Number of CPU cores per task
#SBATCH --mem=100gb                                  # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec
#SBATCH --output=parallel_anemone_plume%j.log        # Standard output/error

export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate idba-ud_assembly_env

#Has IDBA_UD installed and its repo contents

#Run each file through IDBA_UD:
#Look at the length of your reads. Ideally you want a kmer that is 80% of your read length

idba_ud -r /vortexfs1/scratch/selkassas/AnemonePlume_2015_DNA/test-qc/trimmed ERR2021511_merged.fa -o 
/vortexfs1/scratch/selkassas/Background_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_120/ --maxk 120

```

## **4. Run assembled metagenomes through MetaQUAST to assess assembly quality**  
```
#Add Quast to one of your existing conda environments, or create a new one. I just added it to my assembly 
environment. 

conda install -c bioconda quast

#In your slurm script, make sure to direct it to each of the contig files and the scaffold! Example below:

python /vortexfs1/home/selkassas/.conda/envs/idba-ud_assembly_env/bin/metaquast 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/scaffold.fa 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/contig-20.fa 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/contig-40.fa 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/contig-60.fa 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/contig-80.fa 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/contig-100.fa -o 
/vortexfs1/scratch/selkassas/AnemonePlume_1500m_2015_DNA/idbaud_assemblies/assembly_maxkmer_100/

#Transferring output of MetaQuast to Desktop (only way to view). You need the -r flag in order to copy the entire 
  directory and its contents to your desktop.

scp -r selkassas@poseidon.whoi.edu:/vortexfs1/scratch/selkassas/Background_1500m_2015_DNA/idbaud_assemblies/assembl
y_maxkmer_120/quast_results/results_2020_08_30_17_36_21 .

#If you don't want all the output files that MetaQUAST generates, then use this shortened version (courtesy of Dr. 
#Lauren Seyler). My example below: 

python /vortexfs1/home/selkassas/.conda/envs/idba-ud_assembly_env/bin/metaquast contig-20.fa contig-40.fa 
contig-60.fa contig-80.fa contig-100.fa contig-120.fa contig.fa scaffold.fa --no-icarus --space-efficient 
--no-read-stats --max-ref-number 0

#Then did this for each of the assembled metagenomes. 

```

### **It is very important to note that IDBA_UD automatically selects the best assembly by testing different kmer sizes up to, in this case, a kmer size of 120. The "contig.fa" file is the file that you must upload to JGI in the next step. It represents the best assembly that IDBA_UD could create given the parameters!**  


## **5. Annotate each metagenome with the Joint Genome Institute's IMG pipeline**

By uploading each assembly to the Department of Energy's Joint Genome Institute website, you can get all of the 
ORF's annotated with additional taxonomy information. That is why we use JGI vs. something like Prokka or RAST-tk. 
All of the Gold Analysis project ID's are located in the "JGI GOLD Analysis Project ID's.csv" in this repo. It is a
tricky process to upload things to JGI, so I will try my best to walk you through it using Julie's powerpoint 
slide descriptions. Full PDF of slideshow with images of the directions (where to click, etc.) is also in this 
repo named **"JGI_IMG_Pipeline"** and full list of analysis project ID's is below and also in the repo as **"JGI GOLD Analysis Project ID's"**   

1. Rename analysis project name  
2. Study should come up with Axial  
3. Put in Assembly Method, etc  
4. Press SUBMIT, will get an analysis project number  
5. Then go to IMG/MER, Submit Data Set  
6. New Submission: “AP ID” is the number you got from GOLD  
7. Need to specify gene calling, mixed microbial community  
8. “Submit Dataset to IMG”  
9. Need to resolve lack of access to Axial study in GOLD  

Sample Name               	JGI GOLD Analysis Project ID  
Background_1500m_2015_DNA	  Ga0455662  
FS891_Anemone_DNA	          Ga0455663  
FS903_Marker113_DNA	        Ga0456187  
FS904_Marker33_DNA	        Ga0456188  
FS906_Marker113_DNA	        Ga0456189  
FS907_Anemone_DNA	          Ga0456190  
FS908_Marker33_DNA	        Ga0456191  
FS914_Anemone_DNA	          Ga0456192  
FS915_Marker113_DNA	        Ga0456193  
FS917_Marker33_DNA	        Ga0454710  


## **6. Getting annotated metagenomic files off of IMG**

This is actually way trickier than it seems like it would be! 

Steps to get annotated metagenomic files off of IMG  
1. Go to https://img.jgi.doe.gov/  
2. Then click the tab that says “IMG/MER”  
3. Then click the pull down option “IMG/M ER  
4. Then on the left of that page you will see a box that says “IMG Content” and at the bottom you will see “My
Private Datasets. If you click there, everything you’ve ever uploaded to IMG should show up!  
5. Then you have to select the sample name by clicking the box on the left then add the selected file to your 
Genome Cart  
6. Then scroll to the top where it says “my analysis carts” and click “# genomes”  
7. Then select the sample name again by clicking the box on the left and scroll up to the tab that says “upload & 
export & save”  
8. Then scroll down to where it says “download genomes”  
9. Click that and it will bring up a page that at the bottom says “you can check the status here:” Click that.  
10. And it will say download request is being processed. It will be on that page for a few minutes until it gives 
you the link to your download.  
--Alternatively you can just wait for the email notification that the download is finished  


## **7. Converting annotated metagenomes to .gff files with python before concatenating them**  
Note: I have uploaded the gff2seqfeatures.py script to the repo.   
```
Converting IMG data to .gff files before concatenating them for mapping mRNA from SIP experiments to the concatenated metagenomic file

*You use Connor Scannerton’s (Victoria Orphan’s lab) script: gff2seqfeatures.py to extract orf information from the
.gff and .fna files from the IMG annotation of the metagenomes (MG):
  *.gff files contain contig and locus Id, start and stop site, strand info
  *.fna files contains the nucleotide sequence of each contig

Output: Assembled.ffn file with the nucleotide sequence of each orf as called in the gff file.

You then concatenate (cat) all of those assembled.ffn file for each sample into one megafile.

1. Enter a conda environment that you want to do the commands on, or alternatively make a new conda environment
conda activate bowtie2 (the environment I always put random installs in)

2. Then install bcbio and seqIO as part of the bio python package. 
conda install -c bioconda bcbio-gff
conda install -c conda-forge biopython

3. Code to convert to .gff files before concatenating all metagenomes
nano gff2seqfeatures.py

Paste the following code and save: had to remove spaces from lines 7 and 14 and change “locus tag” to [‘ID’] and
then it ran ok.

#!/usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO
from BCBio import GFF
import argparse

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("gff")
    p.add_argument('infasta')
    #p.add_argument('outffn')
    args = p.parse_args()

    seq_dict = SeqIO.to_dict(SeqIO.parse(args.infasta, "fasta"))
    with open(args.gff) as gffp:
        for rec in GFF.parse(gffp, base_dict=seq_dict):
            for feature in rec.features:
                feat = feature.extract(rec.seq)
                try:
                    print(">{}".format(feature.qualifiers['ID'][0]))
                    print(feat)
                except KeyError as e:
                    print("skipping feature as it doesn't have a locus tag", file=sys.stderr)
                    print(feature, file=sys.stderr)


Command to convert to gff file (btw a gff file is an annotation file. Fasta files don’t store annotation data). 
python /<path to python file>/gff2seqfeatures.py *a.gff *a.fna > assembled.ffn

NOTE!! The two files you need are the functional_annotation.gff and the contigs.fna files for this code if you are
struggling with the *.a.gff and .a.fna!

example:

python /vortexfs1/scratch/selkassas/gff2seqfeatures.py Ga0456190_functional_annotation.gff Ga0456190_contigs.fna >
FS907_assembled.ffn 

```


## **8. Concatenating all metagenomic .ffn files into one giant file**
```
1. mkdir concatenation_station
2. Move all of your .ffn files to the concatenation_station
3. use the cat command to concatenate all of your files into one. See my actual code below:

cat assembled_AnemonePlume_1500.ffn assembled_CTDBack.ffn assembled_FS906.ffn assembled_FS915.ffn
assembled_FS891.ffn assembled_FS907.ffn assembled_FS917.ffn assembled_FS903.ffn assembled_FS908.ffn
assembled_FS904.ffn assembled_FS914.ffn > axial_data_concatenated_annotated_metagenomic_file.fasta

#could potentially use srun if you find that it is taking too long. 
```


## **9. Mapping mRNA reads back to concatenated metaG file using Kallisto**
I did a funky technique with bash that allowed me to iterate through all of my samples at once, instead of writing a line of code for each. 
```
#But first, index the concatenated metagenome:

Indexing concatenated metagenome for running with kallisto

#!/bin/bash
#SBATCH --partition=compute                  # Queue selection
#SBATCH --job-name=kallistometamap           # Job name
#SBATCH --mail-type=ALL                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu       # Where to send mail
#SBATCH --ntasks=1                           # Run a single task
#SBATCH --cpus-per-task=36                   # Number of CPU cores per task
#SBATCH --mem=100gb                           # Job memory request
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec
#SBATCH --output=kallistometamap_result.log  # Standard output/error
#export OMP_NUM_THREADS=8

module load anaconda/5.1

source activate bowtie2

cd /vortexfs1/scratch/selkassas/cat

kallisto index -i axial_data_concatenated_annotated_metagenomic_index.idx
axial_data_concatenated_annotated_metagenomic_file.ffn --make-unique

____________
Kallisto Quant 
#Since I have 111 files, I want to map them all at once.
#first do: 

for file in *fastq; do echo ${file} | cut -d _ -f 1-5 | sort ; done > cut_file_names_sorted.txt

#to get a big file with all of the file names without the suffix attached. 
#contents of cut_file_names_sorted.txt:
FS891_RNA
FS903_30_12H_AAGACG_L001
FS903_30_12L_GGATGT_L001
FS903_30_13H_CCTCGG_L001
FS903_30_13L_TTCGCT_L001
FS903_55_12H_AAGGGA_L002
FS903_55_12L_GGACCC_L002
FS903_55_13H_CCTTCA_L002
FS903_55_13L_TTCAGC_L002
FS903_80_12H_AAGACG_L002
FS903_80_12L_GGATGT_L002
FS903_80_13H_CCTCGG_L002
FS903_80_13L_TTCGCT_L002
FS903_RNA
FS904_12H_30_AAGGGA_L007
FS904_12H_55_AAGACG_L007
FS904_12H_80_AAGACG_L007
FS904_12L_30_GGACCC_L007
FS904_12L_55_GGATGT_L007
FS904_12L_80_GGATGT_L007
FS904_13H_30_CCTTCA_L007
FS904_13H_55_CCTCGG_L007
FS904_13H_80_CCGCGG_L007
FS904_13H_80_CCTCGG_L007
FS904_13L_30_TTCAGC_L007
FS904_13L_55_TTCGCT_L007
FS904_13L_80_TTCGCT_L007
FS904_RNA
FS906_80_TP1_12H_GGATGT_L005
FS906_80_TP1_12L_TTCGCT_L005
FS906_80_TP1_13H_AAGGGA_L006
FS906_80_TP1_13L_CCTTCA_L006
FS906_RNA_CCTTCA_L004
FS906_RNA_CCTTCA_L004_RUN2
FS906_RNA_CCTTCA_L005
FS907_80_12INC_AAGACG_L001
FS907_80_12L10_CCTTCA_L001
FS907_80_12L9_AAGGGA_L001
FS907_80_13H6_GGACCC_L001
FS907_80_13H8_TTCAGC_L001
FS907-cDNA_S15
FS907_RNA_TTCAGC_L004
FS907_RNA_TTCAGC_L004_RUN2
FS907_RNA_TTCAGC_L005
FS908_80_12C_INC_11_TTCGCT_L003
FS908_80_12C_INC_9_GGATGT_L003
FS908_80_12L10_CCTCGG_L001
FS908_80_13H7_GGATGT_L001
FS908_80_13INC_TTCGCT_L001
FS908_80_TP1_12H_AAGACG_L006
FS908_80_TP1_12L_CCTCGG_L006
FS908_80_TP1_13H_GGATGT_L006
FS908_80_TP1_13L_TTCGCT_L006
FS908_RNA_CCTCGG_L004
FS908_RNA_CCTCGG_L004_RUN2
FS908_RNA_CCTCGG_L005
FS914-cDNA_S11
FS915-cDNA_S12
FS917-cDNA_S13
LVWS1_80_TP1_12H_AAGGGA_L001
LVWS1_80_TP1_12L_CCTCGG_L003
LVWS1_80_TP1_12L_CCTCGG_L004
LVWS1_80_TP1_13H_AAGACG_L003
LVWS1_80_TP1_13H_AAGACG_L004
LVWS1_80_TP1_13L_CCTTCA_L001
LVWS1_80_TP2_12H_GGACCC_L001
LVWS1_80_TP2_12L_TTCGCT_L003
LVWS1_80_TP2_12L_TTCGCT_L004
LVWS1_80_TP2_13H_GGATGT_L003
LVWS1_80_TP2_13H_GGATGT_L004
LVWS1_80_TP2_13L_TTCAGC_L001
LVWS2_TP1_12H_AAGGGA_L007
LVWS2_TP1_12L_GGACCC_L007
LVWS2_TP1_13H_CCTTCA_L007
LVWS2_TP1_13L_TTCAGC_L007
LVWS4_55_TP1_12H_AAGGGA_L004
LVWS4_55_TP1_12L_CCTTCA_L004
LVWS4_55_TP1_13H_GGACCC_L004
LVWS4_55_TP1_13L_TTCAGC_L004
LVWS4_55_TP2_12L_AAGACG_L003
LVWS4_55_TP2_13H_CCTCGG_L003
LVWS5_30_TP1_12H_AAGGGA_L005
LVWS5_30_TP1_12L_GGACCC_L006
LVWS5_30_TP1_13H_TTCAGC_L006
LVWS5_30_TP1_13L_CCTTCA_L005
LVWS5_30_TP2_12H_AAGGGA_L003
LVWS5_30_TP2_12L_CCTTCA_L003
LVWS5_30_TP2_13H_GGACCC_L003
LVWS5_30_TP2_13L_TTCAGC_L003
LVWS5_55_TP1_12H_GGATGT_L003
LVWS5_55_TP1_12L_AAGACG_L001
LVWS5_55_TP1_13H_CCTCGG_L001
LVWS5_55_TP1_13L_TTCGCT_L003
LVWS5_55_TP2_12H_AAGACG_L003
LVWS5_55_TP2_12L_CCTCGG_L003
LVWS5_55_TP2_13H_GGATGT_L003
LVWS5_55_TP2_13L_TTCGCT_L003
LVWS6_30_TP1_12H_GGACCC_L005
LVWS6_30_TP1_12L_AAGACG_L004
LVWS6_30_TP1_13H_CCTCGG_L004
LVWS6_30_TP1_13L_TTCAGC_L005
LVWS6_30_TP2_12L_AAGGGA_L003
LVWS6_30_TP2_13H_CCTTCA_L003
LVWS6_55_TP1_12H_AAGACG_L005
LVWS6_55_TP1_12L_GGATGT_L004
LVWS6_55_TP1_13H_TTCGCT_L004
LVWS6_55_TP1_13L_CCTCGG_L005
LVWS6_55_TP2_12H_GGACCC_L003
LVWS6_55_TP2_12L_TTCAGC_L003
LVWS6_55_TP2_13H_AAGACG_L003
LVWS6_55_TP2_13L_CCTCGG_L003

#Then do: 

#!/bin/bash
#SBATCH --partition=compute                  # Queue selection
#SBATCH --job-name=kallistoAxialIMG_h18          # Job name
#SBATCH --mail-type=ALL                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu       # Where to send mail
#SBATCH --ntasks=1                           # Run a single task
#SBATCH --cpus-per-task=36                   # Number of CPU cores per task
#SBATCH --mem=180gb                          # Job memory request
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec
#SBATCH --output=kallistoAxialIMG_h18.log        # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate bowtie2

for seq in $(cat ${1})
do
echo Now processing...${seq}
kallisto quant -i axial_data_concatenated_annotated_metagenomic_index.idx -o ${seq}_AxialIMG_mapped_kallisto -b
100 --single -l 151 -s 41 ${seq}_trim_paired_merged_mRNA.fastq 
done 

# -i is the index file input 
# -o is the directory to write output to 
#-b is the number of bootstrap samples (Amy did 100, so I’ll do 100 and can change it later and rerun if necessary
#with fewer bootstraps)
#-fastq input files go at the end. 

#- The files and index are quite large, so will submit in 3 different batches. 

sbatch kallisto_loop.sh head_37_cut_filenames_sorted.txt
Submitted batch job 1752631
sbatch kallisto_loop.sh tail_37_cut_filenames_sorted.txt
Submitted batch job 1752632
sbatch kallisto_loop.sh middle_37_cut_filenames_sorted.txt 
Submitted batch job 1752633

Failed (timed out) so will instead try it in 6 batches. 
	•	split head_37_cut_filenames_sorted.txt into head_18_cut_filenames_sorted.txt and head_19_cut_filenames_sorted.txt
	⁃	sbatch kallisto_looph_18.sh head_18_cut_filenames_sorted.txt
	⁃	Submitted batch job 1765980
	⁃	sbatch kallisto_loop_h19.sh head_19_cut_filenames_sorted.txt
	⁃	Submitted batch job 1765983
	•	split middle_37_cut_filenames_sorted.txt into middle_18_cut_filenames_sorted.txt and middle_19_cut_filenames_sorted.txt
	⁃	sbatch kallisto_loop_m18.sh middle_18_cut_filenames_sorted.txt
	⁃	Submitted batch job 1765991
	⁃	sbatch kallisto_loop_m19.sh middle_19_cut_filenames_sorted.txt
	⁃	Submitted batch job 1765993
	•	split tail_37_cut_filenames_sorted.txt into tail_18_cut_filenames_sorted.txt tail_19_cut_filenames_sorted.txt
	⁃	sbatch kallisto_loop_t18.sh tail_18_cut_filenames_sorted.txt 
	⁃	Submitted batch job 1765999
	⁃	sbatch kallisto_loop_t19.sh tail_19_cut_filenames_sorted.txt
	⁃	Submitted batch job 1766002

All took ~15 hours and ~65gb to run!

```

That is all of the progress I made so far. Will continue updating as I move through the pipeline. 
