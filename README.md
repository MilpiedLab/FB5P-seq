![Screen1](https://github.com/Chuang1118/FB5PE/blob/master/screen1.png?raw=true)


# Pipeline Facs Based 5 Prime End ScRNAseq FB5Pseq

<a name="logo"/>
<div align="center">
<img src="https://github.com/Chuang1118/FB5PE/blob/master/LOGOn.png?raw=true" alt="FB5PE Logo"  width="25%" height="25%" ></img>
</a>
</div>

**Copyright 2019: PMlab**
**This work is distributed under the terms of the GNU General Public License. It is free to use for all purposes.**

----------------

## Introduction

**_corrige moi syntaxe_**

FB5PEseq is a computational pipeline to process Single-cell RNA sequencing data (scRNAseq). It is designed to be run on the cluster systems Gauss (at the Ciml) or Mesocentre (at the Aix Marseille université). The pipeline relies on 4 software principles:

*  The Drop-seq software for get a "digital expression matrix" that will contain integer counts of number of transcripts for each gene in each cell.
*  The Trinity software for the efficient and robust de novo reconstruction of transcriptomes from RNAseq data.
*  The kallisto program for quantifying abundances of transcripts.
*  MiGMAP: IgBlast V-(D)-J mapping tool.
*  Blastn : mapper for constant region.

Some other tools are used for other crucial steps. Here are the main steps of the FB5PEseq pipeline calculations:

1. 1_Preprocessing: This step is done with the Drop-seq tools V1.12. Converting FASTQ to BAM by using Picard, Tagging Bam file (cell barcode and molecule umi) and filtering bad barcode read by TagBamWithReadSequenceExtended and FiLterBAM, respectively.  

2. 2_Preprocessing2: Using Picard tool to convert bam to fastq and sort by queryname for the unmapped bam.

3. 3_Alignment: Mapping the reads on the genome(GRCh38): this step is include merge bam alignment with unmapped bam that have produced by precede step 2, also convert bam to sam.

4. 4_counting: Counting reads in features with htseq-count, converting sam to bam, adding right header with samtools reheader and a python step to merge the BC with a mismatch 1 max.

5. 5_dgeSummary: Extracting Digital Gene Expression (DGE) data from an library is done before using the Drop-seq program DigitalExpression. There are two outputs levels available: tanscript and  UMI. Each level with 2 outputs, the primary is the DGE matrix, with each a row for each gene, and a column for each cell, the secondary analysis is a summary of the DGE matrix on a per-cell level, indicating the number of genes ans transcripts observed.


## FB5PE scRNA-seq workflow

![workflow]( https://github.com/Chuang1118/FB5PE/blob/master/Overview_simplify.png?raw=true )

_Figure 1: An overview of the FB5PE scRNA-seq workflow_

## 
![workflow_detail](https://github.com/Chuang1118/FB5PE/blob/master/Overview.png?raw=true)

_Figure 2: the details of the FB5PE scRNA-seq workflow_



## Resource consideration

Single cell expression analysis is data intensive, and requires substantial computing resources. The pipeline uses the STAR aligner for read mapping, so the memory requirements will scale with the size of the genome. Please look at the STAR manual for the concrete number about the memory requirements. For the human/mouse genome it requires ~ 40Gb of RAM. The pipeline produces temporary files which require a substantial amount of disk space. Please ensure that you have at least 30Gb of disk space per 100 milion sequenced reads.

Important: please make sure that the temporary directory has adequate free space



## Getting Started

To run FB5Pseq on your experimental data, first enter the necessary parameters in the spreadsheet file (see Settings File section), and then from the terminal type. To run the pipeline, you will also need the appropriate genome sequence in fasta format, and the genome annotation in a gtf format

### Prerequisites

Minimum hardware requirements:

* input 22G
* 30 G of RAM
* 1 T of drive space 

Recommended hardware requirement:

* HPC Slurm environment 

Software requirements:

* Git
* Singularity
* Snakemake

### Preparing the environment for the development

What things you need to install the software and how to install them

Create an [environnement virtuel with python3/3.3.6]( http://sametmax.com/les-environnement-virtuels-python-virtualenv-et-virtualenvwrapper/)

then,

```
pip install snakemake
```

### Download scripts, singularity images, Reference genome, data .....

A step by step series of examples that tell you how to get a development env running

#### 1. Use git clone for get all scripts, configuration files and reference genome for Blastn
Get them from Github

```
git clone https://github.com/Chuang1118/FB5PE.git
```

If you scucceed, in the folder cdong/, the folder structure like following:
```
|-- References
|   |-- blast   
|       |-- BCR_CstRegion
|           `-- IMGT_IGH_ConstantRegion_nt.fasta
|           `-- IMGT_IGH_ConstantRegion_nt.fasta.nhr
|           `-- IMGT_IGH_ConstantRegion_nt.fasta.nin
|           `-- IMGT_IGH_ConstantRegion_nt.fasta.nsq
|           `-- IMGT_IGLK_ConstantRegion_nt.fasta
|           `-- IMGT_IGLK_ConstantRegion_nt.fasta.nhr
|           `-- IMGT_IGLK_ConstantRegion_nt.fasta.nin
|           `-- IMGT_IGLK_ConstantRegion_nt.fasta.nsq
|       |-- TCR_CstRegion
|           `-- IMGT_TRA_ConstantRegion_nt.fasta
|           `-- IMGT_TRB_ConstantRegion_nt.fasta
|           `-- IMGT_TR_ConstantRegion_nt.fasta
|           `-- IMGT_TR_ConstantRegion_nt.fasta.nhr
|           `-- IMGT_TR_ConstantRegion_nt.fasta.nin
|           `-- IMGT_TR_ConstantRegion_nt.fasta.nsq
|           `-- IMGT_TRD_ConstantRegion_nt.fasta
|           `-- IMGT_TRG_ConstantRegion_nt.fasta
|-- Run_BCR
|   `-- BCRpreRawdata_PlateSampleList.sh
|   `-- BCRlaunch_1_lines_analysis.sh
|   `-- BCRlaunch_1_lines_analysis_unlock.sh
|   `-- custom_180416_h_HuPhysioB_1_metadata.csv
|   |-- snakafile
|       `-- snakefile_all_BCR.yml
|       `-- dropseq_pipeline.yml
|       `-- trinity_pipeline.ym
|       `-- migmapBCR_pipeline.yml
|       `-- kallisto_pipeline.yml
|    |-- config
|       `-- config_BCR.yml
|       `-- cluster.json
|    |-- barcode
|       `-- barcode_seq_2ndSet.txt
|    |-- script
|       `-- ModifyBarcodeInBam_2ndSetBC.py
|       `-- SplitBamByBarcodeID_2ndSetBC.py
```
#### 2. Use wget for get all singularity images, Star reference genome, reference genome .fa and annotation .gtf

Download reference genome from zenodo and putting them at the right places
```
cd FB5PE/cdong
wget https://zenodo.org/record/3403043/files/starGenome_GRCh38_ERCC92.tar.gz
tar zxvf starGenome_GRCh38_ERCC92.tar.gz
wget https://zenodo.org/record/3403043/files/GRCh38_ERCC92.tar.gz
tar zxvf GRCh38_ERCC92.tar.gz
```

Get singularity images from zenodo and create the folder called "simg" under Run_BCR folder.

```
cd Run_BCR
mkdir simg
cd simg
wget https://zenodo.org/record/3403043/files/bcr_pipe_masters_v2.img
wget https://zenodo.org/record/3403043/files/dropseqtools_1.12.img
wget https://zenodo.org/record/3403043/files/htseq.img
wget https://zenodo.org/record/3403043/files/migmap.img
wget https://zenodo.org/record/3403043/files/tools.img
```

#### 3. Download rawdata .fastq.gz format from the GEO site.

```
cd ..
mkdir Rawdata
cd Rawdata
wget [OPTION]... [URL]...
```

#### 4. Prepare metadata: custom_180416_h_HuPhysioB_1_metadata.csv

In this version, just last 2 lines are necessary, Sample and Plate
```
Keyword,Description,Options
Protocol,Protocol that was used to generate the scRNAseq libraries,Custom5prime
Date,Date of scRNAseq data reception,180416
Species,Species of cells analyzed,H
ProjectH,Project for which human cells were analyzed,HuPhysioB_1
ProjectM,Project for which mouse cells were analyzed,
Description,Short description of the experiment,text
ReferenceGenome,Reference genome to be used for alignment,GRCh38
AdditionalReference,Whether additional references should be used for alignment,
BCR,Whether BCR sequence data should be reconstructed,Yes
TCR,Whether TCR sequence data should be reconstructed,No
IndexSorting,Whether the experiment involves index sorting parameters,No
IndexSortingFiles,Link to the Index Sorting Files Folder,
Read1,Length and expected composition of Read1,"16nt, Well barcode (8nt) and UMI (5 nt), and TAT"
Readi7sequences,i7 library index sequences,
Readi7libraries,i7 indexed libraries names,
Read2,Length and expected composition of Read2,"67nt, gene (scRNAseq)"
PlateBarcodes,Whether the experiment involves plate barcodes parameters,Yes
PlateMapFile,Link to the Plates Map description file,
HTO,Whether HTOs were used for multiplexing,No
HTOsequences,HTO sequences,
HTOlibraries,HTO sample names,
SortingDate,,
SequencingDate,,
SequencingPlatform,,HalioDx NextSeq550
Sample,,10633669_S1;10633670_S2;10633671_S3;10633672_S4;10633673_S5
Plate,,5
```

#### 5. Launch BCRpreRawdata_PlateSampleList.sh

This script take custom_180416_h_HuPhysioB_1_metadata.csv as parameter.

```
./BCRpreRawdata_PlateSampleList.sh custom_180416_h_HuPhysioB_1_metadata.csv
```

First, create a automatic folder Input which each plate folder have the paired end read .fastq.gz that have the symbolic link point to .fastq.gz of Rawdata like below.

```
|-- Input
|   |-- plate1
|       `-- 10633669_S1_R1_001.fastq.gz(symlink)
|       `-- 10633669_S1_R2_001.fastq.gz
|   |-- plate2
|       `-- 10633669_S1_R1_001.fastq.gz
|       `-- 10633669_S1_R2_001.fastq.gz
|   |-- etc…
```
Second, generate to the file plates_samples_list.txt in config folder. 

config/plates_samples_list.txt 
```
plate1   10633669_S1
plate2   10633670_S2
plate3   10633671_S3
plate4   10633672_S4
plate5   10633673_S5
```
#### 6. Settings config file

The configuration file is a YAML file which specifies:

- Locations:

	- The location of the fasta file with the reference genome (must be prepared by the user)

	- The location of a GTF file with genome annotations

	- The location of a fasta file with BCR or TCR constant region

- Migmap parameter (i.e. human/mouse )

- barcodes

Example of config_BCR.yml, when you apply the pipeline, you need change the path **/scratch/cdong/References/......**

```
# path of the file containing the list of plates and samples to be analyzed
plates_samples_list: config/plates_samples_list.txt
 
#Picard: GRCh38_ERCC92.fa and HTSeq : GRCh38_ERCC92.gtf
genomeDir: /scratch/cdong/References/GRCh38_ERCC92

stargenomeDir: /scratch/cdong/References/starGenome_GRCh38_ERCC92

#96 barcodes, idem barcodes following
bcfile: files/barcode_seq_2ndSet.txt

#BCR constant region
ref: /scratch/cdong/References/blast/BCR_CstRegion/IMGT_IGH_ConstantRegion_nt.fasta

#BCR Constant region 
ref2: /scratch/cdong/References/blast/BCR_CstRegion/IMGT_IGLK_ConstantRegion_nt.fasta

#migmap parameter
species: human

barcodes:
 - CGTCTAAT
 - AGACTCGT
 - GCACGTCA
 - TCAACGAC
 - ATTTAGCG
 - ATACAGAC
 - TGCGTAGG
 - TGGAGCTC
 - TGAATACC
 - TCTCACAC
 - TACTGGTA
 - etc...
```

#### 7. Settings example of snakemake using HPC cluster via SLURM scheduler

cluster.json, configuration of each Snakemake rule, I have used 1 cpu and time limited 1h30 by default, and my Input .fastq.gz around 3.3G for Read1 and 1.2G for Read2 were performed using a NextSeq. If your dataset have the huge size, you must change the parameter mem-per-cpu and time limited to run correctly. 
```
{
    "documentation": {
        "cmdline": "Use with snakemake --cluster-config cluster.slurm.cheaha.json --cluster 'sbatch --job-name {cluster.job-name} --ntasks {cluster.ntasks} --cpus-per-task {threads} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'"
    },
    "__default__" : {
        "job-name"       : "c.{rule}",
        "project"        : "b098",
        "partition"      : "skylake",
        "time"           : "01:30:00",
        "ntasks"         : 1,
        "cpus-per-task"  : 1,
        "mem-per-cpu"    : 4000,
        "output"         : "log/%N.%j.%a.out",
        "error"          : "log/%N.%j.%a.err",
        "mail-user"      : "dong@ciml.univ-mrs.fr",
        "mail-type"      : "FAIL"
    },
    "star_map" : {
        "time"           : "05:00:00",
        "mem-per-cpu"    : 30000
    },
    "trinity_on_each_bc_rule" : {
        "time"           : "03:00:00",
        "cpus-per-task"  : 2,
        "mem-per-cpu"    : 8000
    }
}


```
#### 8. Setting Launch BCR file

**_BCRlaunch_1_lines_analysis.sh_** have -- jobs parameter, each plate have 96 barcodes and dataset have 5 plates, so 96*5 = 480 jobs parallel in Trinity rules(i.e. bam_to_fastq_rule, trinity_on_each_bc_rule and so on). 

You need change the parameter -B, A user-bind path specification, to run correctly.

```
#!/bin/bash

module load userspace/all
module load python3/3.6.3
source /home/cdong/snakemake/bin/activate

snakemake --jobs 480 --use-singularity --singularity-args "-B /scratch/cdong:/scratch/cdong" --snakefile snakefile/snakefile_all_BCR.yml --configfile config/config_BCR.yml --cluster-config config/cluster.json --cluster 'sbatch -A {cluster.project} --job-name {cluster.job-name} --ntasks {cluster.ntasks} --cpus-per-task {threads} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'
```
#### 9. Launch the pipeline

If you follow all the procedures correctly, now you can launch the pipeline.

Now Job Launch

```
mkdir log
nohup ./BCRlaunch_1_lines_analysis.sh &
```
After severals hours or 1 day, output like Output directory structure(see following section). The dataset custom_180416_h_HuPhysioB_1 have 5 plates, time execution about 6 hours.


### Output directory structure

```
|-- 1_preprocessing
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 2_preprocessing2
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 3_Alignment
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 4_counting
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 5_dgsummary
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 6_trinity_preprocessing
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 7_trinity
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 8_kallisto
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 9_migmap
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
|-- 10_blast
|   |-- plate1
|       ` -- output_name
|   |-- plate2
|       ` -- output_name
|   |-- etc…
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deploymen

Add additional notes about how to deploy this on a live systemPlease read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests to us.

Versioning
We use SemVer for versioning. For the versions available, see the tags on this repository.

Authors


## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
