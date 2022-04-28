import glob


configfile: "config.yaml"

inputdirectory=config["directory"]
outdirectory=config["directory"]
database=config["database"]
lineageremove=config["lineageremove"]
#barcodefile=config["barcodefile"]

wildcard_constraints:
    sample_demultiplex="barcode\d\d|unclassified",
#    #runid="[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_\d",
#    #runid_demultiplex="[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_\d",
#    runid="fastq_runid_[A-Za-z0-9]+_\d",
#    runid_demultiplex="fastq_runid_[A-Za-z0-9]+_\d",

#For demultiplex
SAMPLES_demultiplex, RUNIDS_demultiplex, = glob_wildcards(inputdirectory+"/basecall/demultiplex/{sample_demultiplex}/{runid_demultiplex}.fastq", followlinks=True)
print("Demultiplex samples")
print(SAMPLES_demultiplex)
print("Demultiplex runids")
print(RUNIDS_demultiplex)

ALL_SAMPLES = SAMPLES_demultiplex
ALL_RUNIDS = RUNIDS_demultiplex


##### target rules #####
rule all:
    input: 
       expand(outdirectory+"/mothur/{sample_demultiplex}_{runid_demultiplex}.fastq", zip, sample_demultiplex=SAMPLES_demultiplex, runid_demultiplex=RUNIDS_demultiplex),
       expand(outdirectory+"/qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
       outdirectory+"/qc/multiqc.html",
       expand(outdirectory+"/mothur/{sample}_{runid}.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
       expand(outdirectory+"/mothur/{sample}_{runid}.trim.fasta", zip, sample=ALL_SAMPLES, runid=ALL_RUNIDS),
       outdirectory+"/mothur/work_dir/merged_results.fasta"



rule symlink_results_demultiplex:
    input:
        inputdirectory+"/basecall/demultiplex/{sample_demultiplex}/{runid_demultiplex}.fastq"
    output:
        outdirectory+"/mothur/{sample_demultiplex}_{runid_demultiplex}.fastq"
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
    wildcard_constraints:
        sample="barcode\d\d|unclassified"
    log:
        "logs/fastqc_pretrim/{sample}_{runid}.log"
    threads: 3
    wrapper:
        "v1.3.2/bio/fastqc"



rule multiqc:
    input:
        expand(outdirectory+"/qc/fastqc_pretrim/{sample}_{runid}_fastqc.zip", zip, sample=SAMPLES_demultiplex, runid=RUNIDS_demultiplex),
    output:
        outdirectory+"/qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.3.2/bio/multiqc"

###Mothur rules
rule mothur_split_qual_fasta:
    input:
        infile=outdirectory+"/mothur/{sample}_{runid}.fastq",
    output:
        outfile=outdirectory+"/mothur/{sample}_{runid}.fasta",
    params:
        fq="./{sample}_{runid}.fastq",
        #oligos=barcodefile,
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
        #mothur "#set.dir(output={params.outdir}); fastq.info(fastq={params.fq}, oligos={params.oligos})"


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
        fasta=expand(outdirectory+"/mothur/{sample}_{runid}.trim.fasta", zip, sample=SAMPLES_demultiplex, runid=RUNIDS_demultiplex),
        refbac=database+"/silva.bacteria.fasta",
        trainsetfasta=database+"/trainset16_022016.pds.fasta",
        trainsettax=database+"/trainset16_022016.pds.tax",
    output:
        mergefasta=outdirectory+"/mothur/work_dir/merged_results.fasta",
        taxon=outdirectory+"/mothur/work_dir/merged_results.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.1.cons.taxonomy",
        shared=outdirectory+"/mothur/work_dir/merged_results.good.unique.filter.unique.precluster.pick.pds.wang.pick.tx.shared",
    params:
        fasta="-".join(expand("./{sample}_{runid}.trim.fasta", zip, sample=SAMPLES_demultiplex, runid=RUNIDS_demultiplex)),
        workingdir=outdirectory+"/mothur/work_dir",
        mothurdir=outdirectory+"/mothur",
        lineageremove=config["lineageremove"],
        groups="-".join([ele.removesuffix(".trim.fasta") for ele in expand("{sample}_{runid}.trim.fasta", zip, sample=SAMPLES_demultiplex, runid=RUNIDS_demultiplex)])
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



