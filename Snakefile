import glob

wildcard_constraints:
    runid="[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_\d",

configfile: "config.yaml"

inputdirectory=config["directory"]
outdirectory=config["directory"]
#lineageremove=config["lineageremove"]


SAMPLES, RUNIDS, = glob_wildcards(inputdirectory+"/basecall/pass/{sample}/pass/{runid}.fastq", followlinks=True)
print("Sequencer passed fastq sample files")
print(SAMPLES)
print(RUNIDS)

SAMPLES_skip, RUNIDS_skip, = glob_wildcards(inputdirectory+"/basecall/skip/{sample_skip}/pass/{runid_skip}.fastq", followlinks=True)
print("Sequencer skipped fastq sample files")
print(SAMPLES_skip)
print(RUNIDS_skip)
ALL_SAMPLES = SAMPLES + SAMPLES_skip
ALL_RUNIDS = RUNIDS + RUNIDS_skip

##### target rules #####
rule all:
    input: 
       expand(outdirectory+"/mothur/{sample}_{runid}.fastq", zip, sample=SAMPLES, runid=RUNIDS),
       expand(outdirectory+"/mothur/{sample_skip}_{runid_skip}.fastq", zip, sample_skip=SAMPLES_skip, runid_skip=RUNIDS_skip),
       expand(outdirectory+"/qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
       outdirectory+"/qc/multiqc.html",
       expand(outdirectory+"/mothur/{sample}_{runid}.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
       #expand(outdirectory+"/mothur/{sample_skip}_{runid_skip}.fasta", zip, sample_skip=SAMPLES_skip, runid_skip=RUNIDS_skip),
       expand(outdirectory+"/mothur/{sample}_{runid}.trim.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
       outdirectory+"/mothur/work_dir/merged_results.fasta"

rule symlink_results_pass:
    input:
        inputdirectory+"/basecall/pass/{sample}/pass/{runid}.fastq"
    output:
        outdirectory+"/mothur/{sample}_{runid}.fastq"
    threads: 1
    shell:
        "ln -s {input} {output}"

rule symlink_results_skip:
    input:
        inputdirectory+"/basecall/skip/{sample}/pass/{runid}.fastq"
    output:
        outdirectory+"/mothur/{sample}_{runid}.fastq"
    threads: 1
    shell:
        "ln -s {input} {output}"


rule fastqc_pretrim:
    input:
        outdirectory+"/mothur/{sample}_{runid}.fastq"
    output:
        html=outdirectory+"/qc/fastqc_pretrim/{sample}_{runid}.html",
        zip=outdirectory+"/qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}_{runid}.log"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"

rule multiqc:
    input:
        expand(outdirectory+"/qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
    output:
        outdirectory+"/qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.77.0/bio/multiqc"

###Mothur rules
rule mothur_split_qual_fasta:
    input:
        infile=outdirectory+"/mothur/{sample}_{runid}.fastq",
    output:
        outfile=outdirectory+"/mothur/{sample}_{runid}.fasta",
    params:
        fq="./{sample}_{runid}.fastq",
        indir=outdirectory+"/mothur/",
        outdir=outdirectory+"/mothur/",
    log:
        "logs/mothur_split/{sample}_{runid}.log"
    threads: 1
    conda:
        "mothur.yaml"
    shell:
        """
        cd {params.indir}
        mothur "#set.dir(output={params.outdir}); fastq.info(fastq={params.fq})"
        """


rule mothur_trim:
    input:
        fasta=outdirectory+"/mothur/{sample}_{runid}.fasta",
        qual=outdirectory+"/mothur/{sample}_{runid}.qual",
    output:
        outdirectory+"/mothur/{sample}_{runid}.trim.fasta"
    params:
        fasta="./{sample}_{runid}.fasta",
        indir=outdirectory+"/mothur/",
        qual="./{sample}_{runid}.qual",
    log:
        "logs/mothur_trim/{sample}_{runid}.log"
    threads: 8
    conda:
        "mothur.yaml"
    shell:
        """
        cd {params.indir}
        mothur "#trim.seqs(fasta={params.fasta}, qfile={params.qual}, qaverage=10, processors=8)"
        """

rule mothur_main:
    input:
        fasta=expand(outdirectory+"/mothur/{sample}_{runid}.trim.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
        refbac="/home/nmhd/databases/nanopore/silva.bacteria.fasta",
        trainsetfasta="/home/nmhd/databases/nanopore/trainset16_022016.pds.fasta",
        trainsettax="/home/nmhd/databases/nanopore/trainset16_022016.pds.tax",
    output:
        mergefasta=outdirectory+"/mothur/work_dir/merged_results.fasta",
        taxon=outdirectory+"/mothur/work_dir/merged_results.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy",
        shared=outdirectory+"/mothur/work_dir/merged_results.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.shared",
    params:
        fasta="-".join(expand("./{sample}_{runid}.trim.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS))
        workingdir=outdirectory+"/mothur/work_dir",
        mothurdir=outdirectory+"/mothur",
        lineageremove=config["lineageremove"],
        groups="-".join([ele.removesuffix(".trim.fasta") for ele in expand("{sample}_{runid}.trim.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS)])
    log:
        "logs/mothur_main/all.log"
    threads: 1
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



