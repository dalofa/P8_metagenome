#!/bin/bash
INPUT=$1
OUTPUT=$2

# Run seqkit stats
seqkit stats --all $1 > $2

# C
