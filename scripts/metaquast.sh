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
