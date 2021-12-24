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
