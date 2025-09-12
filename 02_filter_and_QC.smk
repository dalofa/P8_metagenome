# Workflow for filtering, classification and QC of metagenome assemblies
import os
from glob import glob

input_files = glob("05_assemblies/*.fa")
ASSEMBLIES = glob_wildcards("05_assemblies/{sample}_{assembler}.fa")

#-----------------RULE ALL: REQUEST ALL OUTPUT FILES
rule all:
    input:
        expand(
            "06_tmp/{sample}_{assembler}_filt.fa",
            sample=ASSEMBLIES.sample,
            assembler=ASSEMBLIES.assembler
        ),
        expand(
            "06_tiara/{sample}_{assembler}_classification.tsv",
            sample=ASSEMBLIES.sample,
            assembler=ASSEMBLIES.assembler
        )


##-----------------Filter and classify contigs
# The filtered reads are used for assembly with myloasm
rule filter_contigs:
    input:
        "05_assemblies/{sample}_{assembler}.fa"
    output:
        "06_tmp/{sample}_{assembler}_filt.fa"
    threads:
        32
    conda:
        "envs/seqkit_tiara_env.yaml"
    shell:
        """
          seqkit seq -m 3000 {input} > {output} -j {threads}
        """
##-----------------Filter and classify contigs
# The filtered reads are used for assembly with myloasm
rule tiara_filt:
    input:
        "06_tmp/{sample}_{assembler}_filt.fa"
    output:
        "06_tiara/{sample}_{assembler}_classification.tsv"
    threads:
        32
    conda:
        "envs/seqkit_tiara_env.yaml"
    shell:
        """
          tiara -i {input} -t {threads} -m 3000 -o {output}
          grep -e "prokarya" -e "bacteria" -e "archaea" -e "unknown" {output} | cut -f1 > 06_tiara/bacterial_contigs.txt
          grep -e "eukarya" {output} | cut -f1 > 06_tiara/eukaryote_contigs.txt
          seqkit grep -f 06_tiara/bacterial_contigs.txt {input} > 06_tiara/bacterial_contigs.fasta
          seqkit grep -f 06_tiara/eukaryote_contigs.txt {input} > 06_tiara/eukaryote_contigs.fasta
        """





