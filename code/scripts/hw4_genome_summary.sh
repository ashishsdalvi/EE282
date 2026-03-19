#!/bin/bash
GENOME="/pub/asdalvi/informatics_class/ee282_hw4/data/raw/dmel-all-chromosome-r6.66.fasta.gz"

bioawk -c fastx '{ 
    if(length($seq) > 100000) { 
        n_count = gsub(/[Nn]/, "", $seq); 
        print length($seq), n_count 
    } 
}' $GENOME | awk '{len+=$1; ns+=$2; count++} END {print len, ns, count}'

bioawk -c fastx '{ 
    if(length($seq) <= 100000) { 
        n_count = gsub(/[Nn]/, "", $seq); 
        print length($seq), n_count 
    } 
}' $GENOME | awk '{len+=$1; ns+=$2; count++} END {print len, ns, count}'
