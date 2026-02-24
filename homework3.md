## Summarize Genome Assembly

The script we use is located at /code/scripts/hw3_genome_summary.sh
The way it works is I took a look at the flybase website md5sum.txt file
and located the hash for the all chromosomes fasta file. I then used the md5sum command on
our download and checked if the expected and computed hashes were equal and printed out a statement if they were.
We also use faSize to summarize the fasta file we downloaded.

The file printed out the following when we run it using ./hw3_genome_summary.sh:

```
Integrity Verified
143726002 bases (1152978 N's 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
N count: mean 616.6 sd 6960.7
U count: mean 76242.3 sd 1379508.4
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real
```

As you can see our file is properly downloaded because our script prints out 'Integrity Verified' and
we can see that there are:

143,726,002 bases, 1,152,978 N's, and 1870 sequences. 



## Summarize Annotation File

The script we use is located at /code/scripts/hw3_annotation_summary.sh
The way it works is I took a look at the flybase website md5sum.txt file
and located the hash for the gtf file. I then used the md5sum command on
our download and checked if the expected and computed hashes were equal and printed out a statement if they were.
We also use bioawk to summarize the gtf file we downloaded (sorted features).

The file printed out the following when we run it using ./hw3_annotation_summary.sh:

```
Integrity Verified
 190176 exon
 163377 CDS
  46856 5UTR
  33778 3UTR
  30922 start_codon
  30862 stop_codon
  30836 mRNA
  17872 gene
   3059 ncRNA
    485 miRNA
    365 pseudogene
    312 tRNA
    270 snoRNA
    262 pre_miRNA
    115 rRNA
     32 snRNA
   3508 2L
   3649 2R
   3481 3L
   4226 3R
    114 4
   2704 X
    113 Y
```

The output shows that the file was downloaded properly. Additionally we get the number of feathres. For exons there are 190176, CDS is 163377, 5' UTR and 3'UTR are 46856 and 
33778, start and stop codon are 30922 and 30862, mRNA is 30836, gene is 17872, ncRNA is 3059, miRNA is 485, pseudogene is 365, tRNA is 312, snoRNA is 270, pre-miRNA is 262,
rRNA is 115, snRNA is 32. 

For the genes per chromosome arm 2L is 2508, 2R is 3649, 3L is 3481, 3R is 4226, 4 is 115, X is 2704, and Y is 113. 


NOTE: we downloaded using 
```
# Fasta and GTF respectively
wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/fasta/dmel-all-chromosome-r6.66.fasta.gz
wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/gtf/dmel-all-r6.66.gtf.gz
```
