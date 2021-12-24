#!/bin/bash
#SBATCH --partition=compute                  # Queue selection
#SBATCH --job-name=kallistoAxialIMG_h19          # Job name
#SBATCH --mail-type=ALL                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu       # Where to send mail
#SBATCH --ntasks=1                           # Run a single task
#SBATCH --cpus-per-task=36                   # Number of CPU cores per task
#SBATCH --mem=180gb                          # Job memory request
#SBATCH --time=24:00:00                      # Time limit hrs:min:sec
#SBATCH --output=kallistoAxialIMG_h19.log        # Standard output/error
#export OMP_NUM_THREADS=36

module load anaconda/5.1
source activate bowtie2

for seq in $(cat ${1})
do
echo Now processing...${seq}
kallisto quant -i axial_data_concatenated_annotated_metagenomic_index.idx -o ${seq}_AxialIMG_mapped_kallisto -b 100 --single -l 151 -s 41 ${seq}_trim_paired_merged_mRNA.fastq 
done 

