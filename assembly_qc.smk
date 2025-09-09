# Workflow for Quality control of metagenome assemblies
import os
from glob import glob

input_files = glob("data/np_data/*.fq.gz")
SAMPLES = [os.path.basename(f).replace(".fastq.gz", "").replace(".fq.gz", "") for f in input_files]

print(SAMPLES)
#-----------------RULE ALL: REQUEST ALL OUTPUT FILES
rule all:
    input:


#-----------------QC OF RAW READS
# QC-plots and stats are generated for the each readset
rule nanoplot:
    input:
        "data/np_data/{sample}.fq.gz"
    output:
        html="00_nanoplot/{sample}/NanoPlot-report.html"
    conda:
        "envs/nanoplot_env.yaml"  # Path to your environment YAML file
    shell:
        """
        NanoPlot --fastq {input} --outdir 00_nanoplot/{wildcards.sample} --loglength --plots dot kde
        """

##-----------------Filtlong
# Reads smaller than a 1 kb are excluded for assembly
rule filtlong:
    input:
        "data/np_data/{sample}.fq.gz"
    output:
        "01_filtered_reads/{sample}.fastq.gz"
    params:
        min_length=1000
    conda:
        "envs/filtlong_env.yaml"
    shell:
        """
        filtlong --min_length {params.min_length} {input} | gzip > {output}
        """
##-----------------Metaflye assembly
# The filtered reads are used for assembly with metaflye
rule metaflye:
    input:
        "01_filtered_reads/{sample}.fastq.gz"
    output:
        "02_metaflye/{sample}/assembly.fasta"
    params:
        directory=lambda wildcards, output: os.path.dirname(output[0])
    threads:
         32
    conda:
         "envs/metaflye_env.yaml"
    shell:
        """
          flye --nano-hq {input} --meta --out-dir {params.directory} --threads {threads}
        """


##-----------------metaMDBG
# The filtered reads are used for assembly with metaMDBG
rule metaMDBG:
    input:
        "01_filtered_reads/{sample}.fastq.gz"
    output:
        "03_metaMDBG/{sample}/contigs.fasta.gz"
    params:
        directory=lambda wildcards, output: os.path.dirname(output[0])
    threads:
         32
    conda:
         "envs/metaMDBG_env.yaml"
    shell:
        """
          metaMDBG asm --out-dir {params.directory} --in-ont {input} --threads {threads}
        """
##-----------------myloasm
# The filtered reads are used for assembly with myloasm
rule myloasm:
    input:
        "01_filtered_reads/{sample}.fastq.gz"
    output:
        "04_myloasm/{sample}/assembly_primary.fa"
    params:
        directory=lambda wildcards, output: os.path.dirname(output[0])
    threads:
         32
    conda:
         "envs/myloasm_env.yaml"
    shell:
        """
          myloasm {input} -o {params.directory} -t {threads}
        """
