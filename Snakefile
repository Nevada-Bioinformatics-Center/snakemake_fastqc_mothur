import glob

wildcard_constraints:
    sample="\w+\d+_\w+_\w+\d+_.+_\d",
    runid="fastq_runid_.+_0_0",


SAMPLES, RUNID, = glob_wildcards("../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq.gz", followlinks=True)
print(SAMPLES)
print(RUNID)

##### target rules #####
rule all:
    input: 
       "qc/multiqc.html"


rule fastqc_pretrim:
    input:
        "../snakemake_guppy_basecall_dev/basecall/{sample}/pass/{runid}.fastq.gz"
    output:
        html="qc/fastqc_pretrim/{sample}_{runid}.html",
        zip="qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}_{runid}.log"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip", sample=SAMPLES, runid=RUNID)
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.77.0/bio/multiqc"



