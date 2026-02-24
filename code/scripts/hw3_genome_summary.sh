#!/bin/bash

EXPECTED="ccb86e94117eb4eeaaf70efb6be1b6b9"
ACTUAL=$(md5sum ../../data/raw/dmel-all-chromosome-r6.66.fasta.gz | cut -d' ' -f1)

if [ "$EXPECTED" == "$ACTUAL" ]; then echo "Integrity Verified"; else echo "Integrity Failed"; fi

faSize ../../data/raw/dmel-all-chromosome-r6.66.fasta.gz
