#!/bin/bash

# Move files
#cp ~/net/o-drive/CFB-S-NewBioactiveCompounds/00_Personal-FoldersStaff/dalofa_David_Faurdal/dalofa_20250704_soil/dalofa_20250704_soil/20250704_1302_P2S-01248-A_PBE11238_c46f0c04/fastq_pass data/. -r

# Bundle barcodes
# In this case barcodes 1-12 were used
mkdir data/np_data
for barcode in barcode{01..12}; do
	cat data/fastq_pass/$barcode/*.gz > data/np_data/$barcode.fq.gz
done

rm data/fastq_pass
