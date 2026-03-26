#!/bin/bash

ASSEMBLY="/pub/asdalvi/informatics_class/ee282_hw4/data/processed/dmel_hifi.fasta"
LENGTHS="/pub/asdalvi/informatics_class/ee282_hw4/data/processed/sorted_lengths.txt"

bioawk -c fastx '{ print length($seq) }' $ASSEMBLY | sort -rn > $LENGTHS

TOTAL=$(awk '{sum+=$1} END {print sum}' $LENGTHS)
TARGET=$(echo "$TOTAL / 2" | bc)

awk -v target="$TARGET" '{
    sum+=$1; 
    if (sum >= target) {
        print "N50: " $1; 
        exit;
    }
}' $LENGTHS
