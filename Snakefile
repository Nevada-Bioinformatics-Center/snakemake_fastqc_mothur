import glob

configfile: "config.yaml"

ref_data=config["ref_data"]


#wildcard_constraints:
#    sample="[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+\d+_[A-Za-z0-9]+_\d",
#    runid="[A-Za-z0-9]_[A-Za-z0-9]_\w+_0_\d",
wildcard_constraints:
    runid="[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_\d",


#SAMPLES, RUNIDS, = glob_wildcards("../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq.gz", followlinks=True)
SAMPLES, RUNIDS, = glob_wildcards("../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq", followlinks=True)
print(SAMPLES)
print(RUNIDS)

##### target rules #####
rule all:
    input: 
       #expand("symlink/{sample}_{runid}.fastq.gz", zip, sample=SAMPLES, runid=RUNIDS),
       expand("mothur/{sample}_{runid}.fastq", zip, sample=SAMPLES, runid=RUNIDS),
       expand("mothur/{sample}_{runid}.fasta", zip, sample=SAMPLES, runid=RUNIDS),
       expand("mothur/{sample}_{runid}.trim.fasta", zip, sample=SAMPLES, runid=RUNIDS),
       "qc/multiqc.html",
       "mothur/work_dir/merged_frese.fasta"

rule symlink_results:
    input:
        #"../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq.gz"
        "../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq"
    output:
        #"symlink/{sample}_{runid}.fastq.gz"
        "mothur/{sample}_{runid}.fastq"
    params:
        #"../../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq.gz"
        "../../snakemake_guppy_basecall/basecall/{sample}/pass/{runid}.fastq"
    threads: 1
    shell:
        "ln -s {params} {output}"

rule fastqc_pretrim:
    input:
        #"symlink/{sample}_{runid}.fastq.gz"
        "mothur/{sample}_{runid}.fastq"
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
        expand("qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip", zip, sample=SAMPLES, runid=RUNIDS)
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.77.0/bio/multiqc"

###Mothur rules
rule mothur_split_qual_fasta:
    input:
        #expand("symlink/{sample}_{runid}.fastq", zip, sample=SAMPLES, runid=RUNIDS)
        #expand("symlink/{sample}_{runid}.fastq", zip, sample=SAMPLES, runid=RUNIDS)
        "mothur/{sample}_{runid}.fastq"
    output:
        "mothur/{sample}_{runid}.fasta"
    params:
        "./{sample}_{runid}.fastq"
    log:
        "logs/mothur_split/{sample}_{runid}.log"
    threads: 1
    shell:
        """
        cd mothur
        mothur "#set.dir(output=mothur/split/); fastq.info(fastq={params})"
        """


rule mothur_trim:
    input:
        fasta="mothur/{sample}_{runid}.fasta",
        qual="mothur/{sample}_{runid}.qual",
    output:
        "mothur/{sample}_{runid}.trim.fasta"
    params:
        fasta="./{sample}_{runid}.fasta",
        qual="./{sample}_{runid}.qual",
    log:
        "logs/mothur_trim/{sample}_{runid}.log"
    threads: 1
    shell:
        """
        cd mothur
        mothur "#trim.seqs(fasta={params.fasta}, qfile={params.qual}, qaverage=10, processors=8)"
        """

rule mothur_main:
    input:
        fasta=expand("mothur/{sample}_{runid}.trim.fasta", zip, sample=SAMPLES, runid=RUNIDS),
        refbac=ref_data+"/silva.bacteria.fasta",
        trainsetfasta="/trainset16_022016.pds.fasta",
        trainsettax="/trainset16_022016.pds.tax",
    output:
        "mothur/work_dir/merged_frese.fasta",
        "mothur/work_dir/merged_frese.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy",
        "mothur/work_dir/merged_frese.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.shared",
    params:
        fasta="-".join(expand("./{sample}_{runid}.trim.fasta", zip, sample=SAMPLES, runid=RUNIDS)),
        groups="-".join([ele.removesuffix(".trim.fasta") for ele in expand("{sample}_{runid}.trim.fasta", zip, sample=SAMPLES, runid=RUNIDS)]),
        #groups="sample1-sample2-sample3-sample4-sample5-sample6-sample7"
    log:
        "logs/mothur_main/all.log"
    threads: 1
    shell:
        """
        cd mothur
        mothur "#set.dir(output=work_dir);
	merge.files(input={params.fasta}, output=merged_frese.fasta);
	make.group(fasta={params.fasta}, groups={params.groups});
	screen.seqs(fasta=merged_frese.fasta, group=current, maxambig=0, maxlength=1700, maxhomop=8);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference={input.refbac});
	filter.seqs(fasta=current, vertical=T);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference={input.trainsetfasta}, taxonomy={input.trainsettax}, cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	phylotype(taxonomy=current);
	make.shared(list=current, count=current, label=1);
	classify.otu(list=current, count=current, taxonomy=current, label=1)"
        """



