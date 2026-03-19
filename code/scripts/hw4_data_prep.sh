#!/bin/bash
GENOME="/pub/asdalvi/informatics_class/ee282_hw4/data/raw/dmel-all-chromosome-r6.66.fasta.gz"

bioawk -c fastx '{ if(length($seq) > 100000) print length($seq), gc($seq) }' $GENOME > /pub/asdalvi/informatics_class/ee282_hw4/data/data_large.txt
bioawk -c fastx '{ if(length($seq) <= 100000) print length($seq), gc($seq) }' $GENOME > /pub/asdalvi/informatics_class/ee282_hw4/data/data_smalll.txt
