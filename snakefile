import os
from glob import glob

input_files = glob("data/np_data/*.fq.gz")
SAMPLES = [os.path.basename(f).replace(".fastq.gz", "").replace(".fq.gz", "") for f in input_files]

#-----------------RULE ALL: REQUEST ALL OUTPUT FILES
rule all:
    input:
        expand("00_nanoplot/{sample}/NanoPlot-report.html", sample=SAMPLES)
	expand


#-----------------QC OF RAW READS
rule nanoplot:
    input:
        "data/np_data/{sample}.fq.gz"
    output:
        html="NanoPlot_results/{sample}/NanoPlot-report.html"
    conda:
        "envs/nanoplot_env.yaml"  # Path to your environment YAML file
    shell:
        """
        NanoPlot --fastq {input} --outdir 00_nanoplot/{wildcards.sample} --loglength --plots dot kde
        """

##-----------------QC OF RAW READS



# rule nanoplot_qc:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"
