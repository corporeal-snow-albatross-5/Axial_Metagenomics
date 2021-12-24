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
