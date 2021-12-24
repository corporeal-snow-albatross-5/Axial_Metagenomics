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
