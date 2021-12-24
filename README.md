# Axial_Metagenomics
Metagenomics pipeline for processing Axial Seamount diffuse fluid samples

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
-gff_conversion.sh  
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
