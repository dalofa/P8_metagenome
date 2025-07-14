import os
from glob import glob

input_files = glob("data/np_data/*.fq.gz")
SAMPLES = [os.path.basename(f).replace(".fastq.gz", "").replace(".fq.gz", "") for f in input_files]

#-----------------RULE ALL: REQUEST ALL OUTPUT FILES
rule all:
    input:
        expand("NanoPlot_results/{sample}/NanoPlot-report.html", sample=SAMPLES)


#-----------------QC RAW READS
rule nanoplot:
    input:
        "data/np_data/{sample}.fq.gz"
    output:
        html="NanoPlot_results/{sample}/NanoPlot-report.html"
    shell:
        """
        NanoPlot --fastq {input} --outdir NanoPlot_results/{wildcards.sample} --loglength --loggc
        """



# rule nanoplot_qc:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"