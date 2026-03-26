#!/bin/bash

export PERL5LIB=/data/homezvol3/asdalvi/miniconda3/envs/ee282/lib/perl5/site_perl/5.22.0:$PERL5LIB

BIN="/pub/jje/ee282/bin"
MY_ASM="/pub/asdalvi/informatics_class/ee282_hw4/data/processed/dmel_hifi.fasta"
FB_SCAF="/pub/asdalvi/informatics_class/ee282_hw4/data/raw/dmel-all-chromosome-r6.66.fasta.gz"
PROCESSED_DIR="/pub/asdalvi/informatics_class/ee282_hw4/data/processed"
OUT_DIR="/pub/asdalvi/informatics_class/ee282_hw4/output/figures"

bioawk -c fastx '{print length($seq)}' $MY_ASM | sort -rn > $PROCESSED_DIR/my_hifi.sizes.txt

bioawk -c fastx '{print length($seq)}' $FB_SCAF | sort -rn > $PROCESSED_DIR/fb_scaff.sizes.txt

$BIN/faSplitByN $FB_SCAF /dev/stdout 10 | bioawk -c fastx '{print length($seq)}' | sort -rn > $PROCESSED_DIR/fb_ctg.sizes.txt

$BIN/plotCDF $PROCESSED_DIR/*.sizes.txt $OUT_DIR/contiguity_comparison.png
