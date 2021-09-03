# Snakemake workflow: snakemake_fastqc_mothur

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow  takes the output from the GUPPY workflow and runs fastqc, multiqc, and mothur.

## Authors

* Hans Vasquez-Gross
* Lucas Bishop

## Usage

### Simple

#### Step 1: Install workflow

clone this workflow to your local computer


#### Step 2: Configure workflow

Configure the workflow according to your needs by editing the config.yaml to configure your input basespace PROJECT directory.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n
