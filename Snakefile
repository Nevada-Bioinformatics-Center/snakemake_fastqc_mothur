import os
import glob


configfile: "config.yaml"

inputdirectory=config["directory"]
outdirectory=config["directory"]
database=config["database"]
lineageremove=config["lineageremove"]
#barcodefile=config["barcodefile"]

wildcard_constraints:
    sample="barcode\d\d|unclassified",
    runid="fastq_runid_[A-Za-z0-9]+_\d",

#For demultiplex
SAMPLES_demultiplex, RUNIDS_demultiplex, = glob_wildcards(inputdirectory+"/basecall/demultiplex/{sample}/{runid}.fastq", followlinks=True)
print("Demultiplex samples")
print(SAMPLES_demultiplex)
print("Demultiplex runids")
print(RUNIDS_demultiplex)

ALL_SAMPLES = SAMPLES_demultiplex
ALL_RUNIDS = RUNIDS_demultiplex

def get_symlink_results_demultiplex_input(wildcards):
    fastq_dir = os.path.join(inputdirectory + "/basecall/demultiplex/", wildcards.sample)
    fastq_file = glob.glob(fastq_dir + "/*.fastq")[0] # this assumes there is only ever one fastq file in a directory
    return os.path.join(fastq_dir, fastq_file)

##### target rules #####
rule all:
    input: 
       expand(outdirectory+"/mothur/{sample}.fastq", sample=ALL_SAMPLES),
       expand(outdirectory+"/qc/fastqc_pretrim/{sample}_fastqc.zip", sample=ALL_SAMPLES),
       outdirectory+"/qc/multiqc.html",
       expand(outdirectory+"/mothur/{sample}.fasta", sample=ALL_SAMPLES),
       expand(outdirectory+"/mothur/{sample}.trim.fasta", sample=ALL_SAMPLES),
       outdirectory+"/mothur/work_dir/merged_results.fasta"



rule symlink_results_demultiplex:
    input:
        get_symlink_results_demultiplex_input
    output:
        outdirectory+"/mothur/{sample}.fastq"
    threads: 1
    shell:
        "ln -s {input} {output}"


rule fastqc_pretrim:
    input:
        outdirectory+"/mothur/{sample}.fastq"
    output:
        html=outdirectory+"/qc/fastqc_pretrim/{sample}.html",
        zip=outdirectory+"/qc/fastqc_pretrim/{sample}_fastqc.zip" 
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}.log"
    threads: 4
    wrapper:
        "v1.3.2/bio/fastqc"



rule multiqc:
    input:
        expand(outdirectory+"/qc/fastqc_pretrim/{sample}_fastqc.zip", sample=ALL_SAMPLES),
    output:
        outdirectory+"/qc/multiqc.html"
    params:
        "--cl_config \"{read_count_multiplier: 0.001, read_count_prefix: \"K\", read_count_desc: \"thousands\"}\""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.3.2/bio/multiqc"

###Mothur rules
rule mothur_split_qual_fasta:
    input:
        infile=outdirectory+"/mothur/{sample}.fastq",
    output:
        outfile=outdirectory+"/mothur/{sample}.fasta",
        qual=outdirectory+"/mothur/{sample}.qual",
    params:
        fq="./{sample}.fastq",
        #oligos=barcodefile,
        indir=outdirectory+"/mothur/",
        outdir=outdirectory+"/mothur/",
    log:
        "logs/mothur_split/{sample}.log"
    threads: 16
    conda:
        "mothur.yaml"
    shell:
        """
        cd {params.indir}
        mothur "#set.dir(output={params.outdir}); fastq.info(fastq={params.fq})"
        """
        #mothur "#set.dir(output={params.outdir}); fastq.info(fastq={params.fq}, oligos={params.oligos})"


rule mothur_trim:
    input:
        fasta=outdirectory+"/mothur/{sample}.fasta",
        qual=outdirectory+"/mothur/{sample}.qual",
    output:
        outdirectory+"/mothur/{sample}.trim.fasta"
    params:
        fasta="./{sample}.fasta",
        indir=outdirectory+"/mothur/",
        qual="./{sample}.qual",
    log:
        "logs/mothur_trim/{sample}.log"
    threads: 16
    conda:
        "mothur.yaml"
    shell:
        """
        cd {params.indir}
        mothur "#trim.seqs(fasta={params.fasta}, qfile={params.qual}, qaverage=10, processors=16)"
        """

rule mothur_main:
    input:
        fasta=expand(outdirectory+"/mothur/{sample}.trim.fasta", sample=ALL_SAMPLES),
        refbac=database+"/silva.bacteria.fasta",
        trainsetfasta=database+"/trainset16_022016.pds.fasta",
        trainsettax=database+"/trainset16_022016.pds.tax",
    output:
        mergefasta=outdirectory+"/mothur/work_dir/merged_results.fasta",
        taxon=outdirectory+"/mothur/work_dir/merged_results.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy",
        shared=outdirectory+"/mothur/work_dir/merged_results.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.shared",
    params:
        fasta="-".join(expand("./{sample}.trim.fasta", sample=ALL_SAMPLES)),
        workingdir=outdirectory+"/mothur/work_dir",
        mothurdir=outdirectory+"/mothur",
        lineageremove=config["lineageremove"],
        groups="-".join([ele.removesuffix(".trim.fasta") for ele in expand("{sample}.trim.fasta", sample=ALL_SAMPLES)])
    log:
        "logs/mothur_main/all.log"
    threads: 16
    conda:
        "mothur.yaml"
    shell:
        """
        cd {params.mothurdir}
        mothur "#set.dir(output={params.workingdir});
	merge.files(input={params.fasta}, output=merged_results.fasta);
	make.group(fasta={params.fasta}, groups={params.groups});
	screen.seqs(fasta=merged_results.fasta, group=current, maxambig=0, maxlength=1700, maxhomop=8);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference={input.refbac});
	filter.seqs(fasta=current, vertical=T);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference={input.trainsetfasta}, taxonomy={input.trainsettax}, cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon={params.lineageremove});
	phylotype(taxonomy=current);
	make.shared(list=current, count=current, label=1);
	classify.otu(list=current, count=current, taxonomy=current, label=1) > {log}  2>&1"
        """



