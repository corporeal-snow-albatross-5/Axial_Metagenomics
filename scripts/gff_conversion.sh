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
```{r, eval=FALSE}
1. mkdir concatenation_station
2. Move all of your .ffn files to the concatenation_station
3. use the cat command to concatenate all of your files into one. See my actual code below:

cat assembled_AnemonePlume_1500.ffn assembled_CTDBack.ffn assembled_FS906.ffn assembled_FS915.ffn
assembled_FS891.ffn assembled_FS907.ffn assembled_FS917.ffn assembled_FS903.ffn assembled_FS908.ffn
assembled_FS904.ffn assembled_FS914.ffn > axial_data_concatenated_annotated_metagenomic_file.fasta

#could potentially use srun if you find that it is taking too long. 
