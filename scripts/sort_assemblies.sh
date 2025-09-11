#!/bin/bash
set -euo pipefail

# make sure the output folder exists
mkdir -p 05_assemblies

# Go through 02_metaMDBG
for f in 02_metaflye/*/assembly.fasta; do
    if [[ -f "$f" ]]; then
        subdir=$(basename "$(dirname "$f")")
        echo $subdir
        cp "$f" "05_assemblies/${subdir}_metaflye.fa"
    fi
done


# Go through 03_metaMDBG
for f in 03_metaMDBG/*/contigs.fasta.gz; do
    if [[ -f "$f" ]]; then
        subdir=$(basename "$(dirname "$f")")
        echo $subdir
        zcat "$f" > "05_assemblies/${subdir}_metaMDBG.fa"
    fi
done

# Go through 04_myloasm
for f in 04_myloasm/*/assembly_primary.fa; do
    if [[ -f "$f" ]]; then
        subdir=$(basename "$(dirname "$f")")
	echo $subdir
        cp "$f" "05_assemblies/${subdir}_myloasm.fa"
    fi
done
