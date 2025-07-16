import os
from glob import glob

input_files = glob("data/np_data/*.fq.gz")
SAMPLES = [os.path.basename(f).replace(".fastq.gz", "").replace(".fq.gz", "") for f in input_files]

#-----------------RULE ALL: REQUEST ALL OUTPUT FILES
rule all:
    input:
        expand("00_nanoplot/{sample}/NanoPlot-report.html", sample=SAMPLES),
	expand("01_filtered_reads/{sample}.fastq.gz", sample=SAMPLES)

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
        "envs/nanoplot_env.yaml"
    shell:
        """
        filtlong --min_length {params.min_length} {input} | gzip > {output}
        """
