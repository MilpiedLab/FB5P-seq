# FB5P-seq: FACS-based 5'-end single-cell RNA-seq

<a name="logo"/>
<div align="center">
<img src="https://github.com/Chuang1118/FB5PE/blob/master/LOGOn.png?raw=true" alt="FB5PE Logo"  width="25%" height="25%" ></img>
</a>
</div>

**Copyright 2019: PMlab, Centre d'Immunologie de Marseille-Luminy**
**This work is distributed under the terms of the GNU General Public License. It is free to use for all purposes.**

----------------

## Introduction

FB5P-seq is a computational pipeline to process single-cell RNA sequencing (scRNAseq) data produced with the FB5P-seq protocol designed by the Milpied lab at Centre d'Immunologie de Marseille-Luminy. The pipeline relies on 5 main softwares:

*  The Drop-seq software to get a "digital expression matrix" that will contain integer counts of number of transcripts for each gene in each cell.
*  The Trinity software for the efficient and robust de novo reconstruction of transcriptomes from RNAseq data.
*  The kallisto program for quantifying abundances of transcripts.
*  MiGMAP: IgBlast V-(D)-J mapping tool to identify and filter T cell receptor (TCR) or B cell receptor (BCR) sequences.
*  Blastn : mapper for constant regions of TCR or BCR.

Some other tools are used at different steps. Here are the main steps of the FB5P-seq pipeline:

1. 1_Preprocessing: This step is done with the Drop-seq tools V1.12. Converting FASTQ to BAM by using Picard, Tagging Bam file (cell barcode and molecule umi) and filtering bad barcode reads by TagBamWithReadSequenceExtended and FiLterBAM, respectively.  

2. 2_Preprocessing2: Using Picard tool to convert bam to fastq and sort by queryname for the unmapped bam.

3. 3_Alignment: Mapping the reads on the genome (here, we use the human GRCh38 with additional references for ERCC spike-ins): this step includes the merging of bam alignment with the unmapped bam that was produced in step 2, and also converts bam to sam.

4. 4_counting: Counting reads in features with htseq-count, converting sam to bam, adding right header with samtools reheader and a python step to merge the cell barcodes with 1 mismatch.

5. 5_dgeSummary: Extracting Digital Gene Expression (DGE) data is done using the Drop-seq program DigitalExpression. There are two output types available: transcript and  UMI. Each type with 2 outputs, the primary output is the DGE matrix, with a row for each gene, and a column for each cell; the secondary output is a summary of the DGE matrix on a per-cell level, indicating the number of genes and transcripts observed.


## FB5P-seq pipeline workflow

![workflow]( https://github.com/Chuang1118/FB5PE/blob/master/Overview_simplify.png?raw=true )

_Figure 1: overview of the FB5P-seq pipeline workflow_
**Copyright 2019: PMlab, Centre d'Immunologie de Marseille-Luminy**

## 
![workflow_detail](https://github.com/Chuang1118/FB5PE/blob/master/Overview.png?raw=true)

_Figure 2: detailed view of the FB5P-seq pipeline workflow_
**Copyright 2019: PMlab, Centre d'Immunologie de Marseille-Luminy**


## Resource consideration

Single cell gene expression expression analysis requires substantial computing resources. The pipeline uses the STAR aligner for read mapping, so the memory requirements will scale with the size of the genome. Please look at the STAR manual for specific memory requirements. For the human or mouse genome it requires ~ 40 Gb of RAM. The pipeline produces temporary files which require a substantial amount of disk space. Please ensure that you have at least 30 Gb of disk space per 100 million sequenced reads.

Important: make sure that the temporary directory has adequate free space.



## Getting Started

To run FB5P-seq on your experimental data, first enter the necessary parameters in the spreadsheet file (see Settings File section), and then from the terminal type. To run the pipeline, you will also need the appropriate genome sequence in fasta format, and the genome annotation in a gtf format.

### Prerequisites

Minimum hardware requirements:

* input 29 Gb
* 30 Gb of RAM
* 1.5 Tb of drive space (output folder 1.2 Tb)

Recommended hardware requirement:
We provide two versions of the pipeline depending on the type of hardware:

* HPC Slurm environment 
* Multiple cores Server > 10CPU/20Threads, RAM > 100 GB

Software requirements:

* Git
* Singularity
* Snakemake

### 2 Versions
* Total parallel jobs run in HPC Slurm
```
112 nodes Dell PowerEdge C6420
32 cores/node, processor Intel® Xeon® Gold 6142 (Sky Lake) at 2.6 GHz 
192 Gb RAM/node
```
* Semi parallel jobs run in Multiple cores Server
```
Cores: 20CPU/40Threads
RAM: 190GB
```
## Version 1: HPC Slurm hardware

### Preparing the environment for the development

Things you need to install the software and how to install them.

Create a [virtual environment with python3]( http://sametmax.com/les-environnement-virtuels-python-virtualenv-et-virtualenvwrapper/),within Python Virtual Environment install snakemake 

```
$ pip install --user virtualenv
$ virtualenv fb5p -p /usr/bin/python3.5
$ source fb5p/bin/activate
(fb5p)$ pip3 install snakemake
```

### Download scripts, singularity images, Reference genome, data .....

A step by step series of examples that tell you how to get a development env running.

#### 1. Use git clone to get all scripts, configuration files and reference genome for Blastn
Get them from Github

```
git clone https://github.com/MilpiedLab/FB5P-seq.git
```

If you succeed, in the folder FB5P-seq/, the folder structure will be the following:
```
├── cdong
│   ├── References
│   │   └── blast
│   │       ├── BCR_CstRegion
│   │       │   ├── IMGT_IGH_ConstantRegion_aa.fasta
│   │       │   ├── IMGT_IGH_ConstantRegion_nt.fasta
│   │       │   ├── IMGT_IGH_ConstantRegion_nt.fasta.nhr
│   │       │   ├── IMGT_IGH_ConstantRegion_nt.fasta.nin
│   │       │   ├── IMGT_IGH_ConstantRegion_nt.fasta.nsq
│   │       │   ├── IMGT_IGLK_ConstantRegion_aa.fasta
│   │       │   ├── IMGT_IGLK_ConstantRegion_nt.fasta
│   │       │   ├── IMGT_IGLK_ConstantRegion_nt.fasta.nhr
│   │       │   ├── IMGT_IGLK_ConstantRegion_nt.fasta.nin
│   │       │   └── IMGT_IGLK_ConstantRegion_nt.fasta.nsq
│   │       └── TCR_CstRegion
│   │           ├── IMGT_TRA_ConstantRegion_nt.fasta
│   │           ├── IMGT_TRB_ConstantRegion_nt.fasta
│   │           ├── IMGT_TR_ConstantRegion_nt.fasta
│   │           ├── IMGT_TR_ConstantRegion_nt.fasta.nhr
│   │           ├── IMGT_TR_ConstantRegion_nt.fasta.nin
│   │           ├── IMGT_TR_ConstantRegion_nt.fasta.nsq
│   │           ├── IMGT_TRD_ConstantRegion_nt.fasta
│   │           └── IMGT_TRG_ConstantRegion_nt.fasta
│   ├── Run_BCR
│   │   ├── barcode
│   │   │   └── barcode_seq_2ndSet.txt
│   │   ├── BCRlaunch_1_lines_analysis.sh
│   │   ├── BCRlaunch_1_lines_analysis_unlock.sh
│   │   ├── BCRpreRawdata_PlateSampleList.sh
│   │   ├── config
│   │   │   ├── cluster.json
│   │   │   └── config_BCR.yml
│   │   ├── custom_180416_h_HuPhysioB_2_metadata.csv
│   │   ├── script
│   │   │   ├── ModifyBarcodeInBam_2ndSetBC.py
│   │   │   └── SplitBamByBarcodeID_2ndSetBC.py
│   │   └── snakefile
│   │       ├── dropseq_pipeline.yml
│   │       ├── kallisto_pipeline.yml
│   │       ├── migmapBCR_pipeline.yml
│   │       ├── snakefile_all_BCR.yml
│   │       └── trinity_pipeline.yml
│   └── Run_TCR
│       ├── barcode
│       │   └── barcode_seq_2ndSet.txt
│       ├── config
│       │   ├── cluster.json
│       │   └── config_TCR.yml
│       ├── custom_190220_h_HuTcells-AR_metadata.csv
│       ├── preRawdata_PlateSampleList_TCR.sh
│       ├── script
│       │   ├── ModifyBarcodeInBam_2ndSetBC.py
│       │   └── SplitBamByBarcodeID_2ndSetBC.py
│       ├── snakefile
│       │   ├── dropseq_pipeline.yml
│       │   ├── kallisto_pipeline.yml
│       │   ├── migmapTCR_pipeline.yml
│       │   ├── snakefile_all_TCR.yml
│       │   └── trinity_pipeline.yml
│       ├── TCRlaunch_1_lines_analysis.sh
│       └── TCRlaunch_1_lines_analysis_unlock.sh
├── LOGOn.png
├── Overview.png
├── Overview_simplify.png
├── README.md
└── screen1.png
```
#### 2. Use wget for get all singularity images, Star reference genome, reference genome .fa and annotation .gtf

Download reference genome (human GRCh38 plus additional references for ERCC spike-ins) from zenodo and put them at the right place
```
cd FB5P-seq/cdong/References/
wget https://zenodo.org/record/3403043/files/starGenome_GRCh38_ERCC92.tar.gz
tar zxvf starGenome_GRCh38_ERCC92.tar.gz
wget https://zenodo.org/record/3403043/files/GRCh38_ERCC92.tar.gz
tar zxvf GRCh38_ERCC92.tar.gz
```

Get singularity images from zenodo and create the folder called "simg" under Run_BCR folder.

```
cd ../Run_BCR/
mkdir simg
cd simg
wget https://zenodo.org/record/3403043/files/bcr_pipe_masters_v2.img
wget https://zenodo.org/record/3403043/files/dropseqtools_1.12.img
wget https://zenodo.org/record/3403043/files/htseq.img
wget https://zenodo.org/record/3403043/files/migmap.img
wget https://zenodo.org/record/3403043/files/tools.img
```

#### 3. Download or copy rawdata in .fastq.gz format.

```
cd ..
mkdir Rawdata
cd Rawdata
wget [OPTION]... [URL]...
```

#### 4. Prepare metadata: 

example file in Run_BCR: custom_180416_h_HuPhysioB_2_metadata.csv
In this version, just the last 2 lines are necessary, Sample and Plate
```
Keyword,Description,Options
Protocol,Protocol that was used to generate the scRNAseq libraries,Custom5prime
Date,Date of scRNAseq data reception,180416
Species,Species of cells analyzed,H
ProjectH,Project for which human cells were analyzed,HuPhysioB_2
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
Sample,,10633663_S1;10633664_S2;10633665_S3;10633666_S4;10633667_S5;10633668_S6
Plate,,6
```

#### 5. Launch BCRpreRawdata_PlateSampleList.sh

For this example, the script takes custom_180416_h_HuPhysioB_2_metadata.csv as input parameter.

```
cd ..
./BCRpreRawdata_PlateSampleList.sh custom_180416_h_HuPhysioB_2_metadata.csv
```

First, automatically create a folder Input in which each plate folder contains the paired end read files .fastq.gz that have the symbolic link point to .fastq.gz of Rawdata like below.

```
|-- Input
|   |-- plate1
|       `-- 10633663_S1_R1_001.fastq.gz(symlink)
|       `-- 10633663_S1_R2_001.fastq.gz
|   |-- plate2
|       `-- 10633664_S2_R1_001.fastq.gz
|       `-- 10633664_S2_R2_001.fastq.gz
|   |-- etc…
```
Second, generate the file plates_samples_list.txt in config folder. 

config/plates_samples_list.txt 
```
plate1   10633663_S1
plate2   10633664_S2
plate3   10633665_S3
plate4   10633666_S4
plate5   10633667_S5
plate6   10633668_S6
```
#### 6. Settings config file

The configuration file is a YAML file which specifies:

- Locations:

	- The location of the fasta file with the reference genome (must be prepared by the user)

	- The location of a GTF file with genome annotations

	- The location of a fasta file with BCR or TCR constant region

- Migmap parameter (i.e. human/mouse )

- barcodes

Example of config_BCR.yml, when you apply the pipeline, you need change the path **/scratch/cdong/References/......** to your specified path where you have installed the FB5P-seq folders.

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

#### 7. Example settings of snakemake for running on HPC cluster via SLURM scheduler

_cluster.json_, configuration of each Snakemake rule, I have used 2 cpu ("mem-per-cpu":16000) and time limited to 1h30 by default, for my Input .fastq.gz around 3.8 Gb for Read1 and 1.4 Gb for Read2 (performed on NextSeq550 with HighOutput 75-cycle flow cell). Depending on the size of your dataset, you will have to change the parameters mem-per-cpu and time limited to run correctly. 
```
{
    "documentation": {
        "cmdline": "Use with snakemake --cluster-config cluster.slurm.cheaha.json --cluster 'sbatch --job-name {cluster.job-name} --ntasks {cluster.ntasks} --cpus-per-task {threads} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'"
    },
    "__default__" : {
        "job-name"       : "e.{rule}",
        "project"        : "b098",
        "partition"      : "skylake",
        "time"           : "01:30:00",
        "ntasks"         : 1,
        "cpus-per-task"  : 1,
        "mem-per-cpu"    : 8000,
        "output"         : "log/%N.%j.%a.out",
        "error"          : "log/%N.%j.%a.err",
        "mail-user"      : "dong@ciml.univ-mrs.fr",
        "mail-type"      : "FAIL"
    },
    "picard_sort_unmapped" : {
        "mem-per-cpu"    : 16000
    },
    "picard_sort_mapped" : {
        "mem-per-cpu"    : 16000
    },
    "star_map" : {
        "time"           : "05:00:00",
        "mem-per-cpu"    : 30000
    },
    "merged_alignment"   : {
        "mem-per-cpu"    : 16000
    },
    "trinity_on_each_bc_rule" : {
        "time"           : "03:00:00",
        "cpus-per-task"  : 2,
        "mem-per-cpu"    : 8000
    }
}
~      


```
#### 8. Setting Launch BCR file

**_BCRlaunch_1_lines_analysis.sh_** have -- jobs parameter, each plate have 96 barcodes and dataset have 6 plates, so 96*6 = 576(around 600) jobs in parallel for Trinity rules (i.e. bam_to_fastq_rule, trinity_on_each_bc_rule and so on). 

You will need to change the parameter -B to specify the path to your FB5P-seq run folder, to run correctly.

```
#!/bin/bash

module load userspace/all
module load python3/3.6.3
source /home/cdong/snakemake/bin/activate

snakemake --jobs 600 --use-singularity --singularity-args "-B /scratch/cdong:/scratch/cdong" --snakefile snakefile/snakefile_all_BCR.yml --configfile config/config_BCR.yml --cluster-config config/cluster.json --cluster 'sbatch -A {cluster.project} --job-name {cluster.job-name} --ntasks {cluster.ntasks} --cpus-per-task {threads} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'
```
#### 9. Launch the pipeline

If you follow all the procedures correctly, now you can launch the pipeline.

Now Job Launch

```
mkdir log
nohup ./BCRlaunch_1_lines_analysis.sh &
```
After severals hours or 1 day, the output will be found in the Output directory (see following section). My dataset custom_180416_h_HuPhysioB_2 has scRNAseq data for 6 plates of 96 cells, time execution is 05:07:33.


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

### Log file

nohup

```
[cdong@login02 Run_2019_09_11]$ head -n 50 nohup.out 
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 600
Job counts:
	count	jobs
	576	add_column_to_tsv_rule
	6	add_right_header
	1	all
	6	bam_to_fastq
	576	bam_to_fastq_rule
	6	bam_to_sam
	6	blastn_csr_HeavyCstRegion_rule
	6	blastn_csr_LightCstRegion_rule
	6	detect_errors
	576	fastq_modification_for_trinity_rule
	6	fastq_to_bam
	6	filter_bad_bc_read
	6	htseq_count
	576	index_transcriptome_for_kallisto_rule
	6	merge_kallisto_abundance_rule
	6	merge_trinity_results_rule
	6	merged_alignment
	6	migmap_rule
	6	modify_barcode
	6	modify_barcode_in_bam
	576	modify_trinity_results_rule
	6	picard_sort_mapped
	6	picard_sort_unmapped
	576	quantify_tpm_kallisto_rule
	6	sam_to_bam
	6	split_bam_by_barcode
	6	star_map
	6	summary_matrix_transcript
	6	summary_matrix_umi
	6	tag_cell_barcode
	6	tag_molecular_umi
	576	trinity_on_each_bc_rule
	4177

[Thu Sep 12 15:41:43 2019]
rule fastq_to_bam:
    input: Input/plate1/10633663_S1_R1_001.fastq.gz, Input/plate1/10633663_S1_R2_001.fastq.gz
    output: output/1_preprocessing/plate1/10633663_S1.bam
    jobid: 4171
    wildcards: plate=plate1, sample=10633663_S1

Submitted job 4171 with external jobid 'Submitted batch job 1425241'.

[Thu Sep 12 15:41:43 2019]
rule fastq_to_bam:
[cdong@login02 Run_2019_09_11]$ tail -n 50 nohup.out 
Finished job 35.
4172 of 4177 steps (100%) done
[Thu Sep 12 20:35:55 2019]
Finished job 1200.
4173 of 4177 steps (100%) done

[Thu Sep 12 20:35:55 2019]
rule summary_matrix_transcript:
    input: output/5_dgsummary/plate6/10633668_S6_merged_clean_HTSeqCount_modifiedBC.bam
    output: output/5_dgsummary/plate6/10633668_S6_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate6/10633668_S6_dge_results_HTSeqCount_modifiedBC_transcript.txt
    jobid: 42
    wildcards: plate=plate6, sample=10633668_S6

Submitted job 42 with external jobid 'Submitted batch job 1429776'.

[Thu Sep 12 20:35:55 2019]
rule summary_matrix_umi:
    input: output/5_dgsummary/plate6/10633668_S6_merged_clean_HTSeqCount_modifiedBC.bam
    output: output/5_dgsummary/plate6/10633668_S6_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate6/10633668_S6_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
    jobid: 36
    wildcards: plate=plate6, sample=10633668_S6

Submitted job 36 with external jobid 'Submitted batch job 1429777'.
[Thu Sep 12 20:46:15 2019]
Finished job 5.
4174 of 4177 steps (100%) done
[Thu Sep 12 20:48:36 2019]
Finished job 42.
4175 of 4177 steps (100%) done
[Thu Sep 12 20:49:16 2019]
Finished job 36.
4176 of 4177 steps (100%) done

[Thu Sep 12 20:49:16 2019]
localrule all:
    input: output/9_migmap/plate1/10633663_S1_migmap_output_filtered.csv, output/9_migmap/plate2/10633664_S2_migmap_output_filtered.csv, output/9_migmap/plate3/10633665_S3_migmap_output_filtered.csv, output/9_migmap/plate4/10633666_S4_migmap_output_filtered.csv, output/9_migmap/plate5/10633667_S5_migmap_output_filtered.csv, output/9_migmap/plate6/10633668_S6_migmap_output_filtered.csv, output/10_blast/plate1/10633663_S1_blastn_csr_HeavyCstRegion.out, output/10_blast/plate2/10633664_S2_blastn_csr_HeavyCstRegion.out, output/10_blast/plate3/10633665_S3_blastn_csr_HeavyCstRegion.out, output/10_blast/plate4/10633666_S4_blastn_csr_HeavyCstRegion.out, output/10_blast/plate5/10633667_S5_blastn_csr_HeavyCstRegion.out, output/10_blast/plate6/10633668_S6_blastn_csr_HeavyCstRegion.out, output/10_blast/plate1/10633663_S1_blastn_csr_LightCstRegion.out, output/10_blast/plate2/10633664_S2_blastn_csr_LightCstRegion.out, output/10_blast/plate3/10633665_S3_blastn_csr_LightCstRegion.out, output/10_blast/plate4/10633666_S4_blastn_csr_LightCstRegion.out, output/10_blast/plate5/10633667_S5_blastn_csr_LightCstRegion.out, output/10_blast/plate6/10633668_S6_blastn_csr_LightCstRegion.out, output/8_kallisto/plate1/10633663_S1_kallisto_merged_abundance.tsv, output/8_kallisto/plate2/10633664_S2_kallisto_merged_abundance.tsv, output/8_kallisto/plate3/10633665_S3_kallisto_merged_abundance.tsv, output/8_kallisto/plate4/10633666_S4_kallisto_merged_abundance.tsv, output/8_kallisto/plate5/10633667_S5_kallisto_merged_abundance.tsv, output/8_kallisto/plate6/10633668_S6_kallisto_merged_abundance.tsv, output/7_trinity/plate1/10633663_S1_all_trinity/trinity_results.fasta, output/7_trinity/plate2/10633664_S2_all_trinity/trinity_results.fasta, output/7_trinity/plate3/10633665_S3_all_trinity/trinity_results.fasta, output/7_trinity/plate4/10633666_S4_all_trinity/trinity_results.fasta, output/7_trinity/plate5/10633667_S5_all_trinity/trinity_results.fasta, output/7_trinity/plate6/10633668_S6_all_trinity/trinity_results.fasta, output/7_trinity/plate1/10633663_S1_all_trinity/merge_trinity_results_rule.log, output/7_trinity/plate2/10633664_S2_all_trinity/merge_trinity_results_rule.log, output/7_trinity/plate3/10633665_S3_all_trinity/merge_trinity_results_rule.log, output/7_trinity/plate4/10633666_S4_all_trinity/merge_trinity_results_rule.log, output/7_trinity/plate5/10633667_S5_all_trinity/merge_trinity_results_rule.log, output/7_trinity/plate6/10633668_S6_all_trinity/merge_trinity_results_rule.log, output/5_dgsummary/plate1/10633663_S1_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate2/10633664_S2_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate3/10633665_S3_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate4/10633666_S4_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate5/10633667_S5_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate6/10633668_S6_dge_summary_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate1/10633663_S1_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate2/10633664_S2_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate3/10633665_S3_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate4/10633666_S4_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate5/10633667_S5_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate6/10633668_S6_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt, output/5_dgsummary/plate1/10633663_S1_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate2/10633664_S2_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate3/10633665_S3_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate4/10633666_S4_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate5/10633667_S5_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate6/10633668_S6_dge_summary_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate1/10633663_S1_dge_results_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate2/10633664_S2_dge_results_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate3/10633665_S3_dge_results_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate4/10633666_S4_dge_results_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate5/10633667_S5_dge_results_HTSeqCount_modifiedBC_transcript.txt, output/5_dgsummary/plate6/10633668_S6_dge_results_HTSeqCount_modifiedBC_transcript.txt
    jobid: 0

[Thu Sep 12 20:49:16 2019]
Finished job 0.
4177 of 4177 steps (100%) done
Complete log: /scratch/cdong/Projet_FB5PE/custom_180416_h_HuPhysioB_2/Run_2019_09_11/.snakemake/log/2019-09-12T154137.140191.snakemake.log
Duration dropseq: 0:00:00.017415
Duration trinity: 0:00:00.026882
Duration kallisto: 0:00:00.032312
Duration migmap_BCR: 0:00:00.036807
***** Reading file with plates and samples definition
  List of plates to analyse: ['plate1', 'plate2', 'plate3', 'plate4', 'plate5', 'plate6']
  List of samples to analyse: ['10633663_S1', '10633664_S2', '10633665_S3', '10633666_S4', '10633667_S5', '10633668_S6']
***** End reading file with plates and samples definition
```
### Output important for downstream analyses
```
├── 10_blast
│   ├── plate1
│   │   ├── 10633663_S1_blastn_csr_HeavyCstRegion.out
│   │   └── 10633663_S1_blastn_csr_LightCstRegion.out
│   ├── plate2
│   │   ├── 10633664_S2_blastn_csr_HeavyCstRegion.out
│   │   └── 10633664_S2_blastn_csr_LightCstRegion.out
│   ├── plate3
│   │   ├── 10633665_S3_blastn_csr_HeavyCstRegion.out
│   │   └── 10633665_S3_blastn_csr_LightCstRegion.out
│   ├── plate4
│   │   ├── 10633666_S4_blastn_csr_HeavyCstRegion.out
│   │   └── 10633666_S4_blastn_csr_LightCstRegion.out
│   ├── plate5
│   │   ├── 10633667_S5_blastn_csr_HeavyCstRegion.out
│   │   └── 10633667_S5_blastn_csr_LightCstRegion.out
│   └── plate6
│       ├── 10633668_S6_blastn_csr_HeavyCstRegion.out
│       └── 10633668_S6_blastn_csr_LightCstRegion.out
├── 5_dgsummary
│   ├── plate1
│   │   ├── 10633663_S1_dge_results_HTSeqCount_modifiedBC_transcript.txt
│   │   └── 10633663_S1_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
│   ├── plate2
│   │   ├── 10633664_S2_dge_results_HTSeqCount_modifiedBC_transcript.txt
│   │   └── 10633664_S2_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
│   ├── plate3
│   │   ├── 10633665_S3_dge_results_HTSeqCount_modifiedBC_transcript.txt
│   │   └── 10633665_S3_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
│   ├── plate4
│   │   ├── 10633666_S4_dge_results_HTSeqCount_modifiedBC_transcript.txt
│   │   └── 10633666_S4_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
│   ├── plate5
│   │   ├── 10633667_S5_dge_results_HTSeqCount_modifiedBC_transcript.txt
│   │   └── 10633667_S5_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
│   └── plate6
│       ├── 10633668_S6_dge_results_HTSeqCount_modifiedBC_transcript.txt
│       └── 10633668_S6_dge_results_HTSeqCount_modifiedBC_UMI_EditDistance0_IslamFilter.txt
├── 8_kallisto
│   ├── plate1
│   │   └── 10633663_S1_kallisto_merged_abundance.tsv
│   ├── plate2
│   │   └── 10633664_S2_kallisto_merged_abundance.tsv
│   ├── plate3
│   │   └── 10633665_S3_kallisto_merged_abundance.tsv
│   ├── plate4
│   │   └── 10633666_S4_kallisto_merged_abundance.tsv
│   ├── plate5
│   │   └── 10633667_S5_kallisto_merged_abundance.tsv
│   └── plate6
│       └── 10633668_S6_kallisto_merged_abundance.tsv
└── 9_migmap
    ├── plate1
    │   └── 10633663_S1_migmap_output_filtered.csv
    ├── plate2
    │   └── 10633664_S2_migmap_output_filtered.csv
    ├── plate3
    │   └── 10633665_S3_migmap_output_filtered.csv
    ├── plate4
    │   └── 10633666_S4_migmap_output_filtered.csv
    ├── plate5
    │   └── 10633667_S5_migmap_output_filtered.csv
    └── plate6
        └── 10633668_S6_migmap_output_filtered.csv
```
## Version 2: Multiple cores Server

This version is similar to Version 1:  HPC Slurm hardware, just change the config_TCR.yml file, including Path of genome reference. Here the example is for running on a human T cell FB5P-seq dataset.
```
#Path of the file containing the list of plates and samples to be analyzed
plates_samples_list: config/plates_samples_list.txt

genomeDir: /scratch/cdong/References/GRCh38_ERCC92

stargenomeDir: /scratch/cdong/References/starGenome_GRCh38_ERCC92

bcfile: barcode/barcode_seq_2ndSet.txt

ref_TCR: /scratch/cdong/References/blast/TCR_CstRegion/IMGT_TR_ConstantRegion_nt.fasta

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

```
#### Launch the pipeline
snakemake -rn : dryrun print a summary of the DAG of jobs.
if you don't get any error, you want to launch the jobs delete "n".
```
$ pwd
/FB5P-seq/cdong/Run_TCR
$ source fb5p/bin/activate
(fb5p)$ nohup snakemake -j 30 -rn --use-singularity -s snakefile/snakefile_all_TCR.yml &
(fb5p)$ nohup snakemake -j 30 -r --use-singularity -s snakefile/snakefile_all_TCR.yml &
```

## Cleaning up intermediate files [Optional]
In the current workflow, we have not implemented the method to mark output files as temporary or protected files in snakemake.
To free up disk space after executing the pipeline, you can use the following command lines and save only the important outputs for downstream analyses (see previous section for details).
```
cd output
rm -r 1_preprocessing
rm -r 2_preprocessing2
rm -r 3_Alignment
rm -r 4_counting
rm 5_dgsummary/plate*/*.bam
rm -r 6_trinity_preprocessing
rm -r 7_trinity/plate*/BAMbyBC
rm -r 8_kallisto/plate*/BAMbyBC
```

## Authors

**Chuang Dong, Pierre Milpied, Iñaki Cervera-Marzal**

Contact: dong@ciml.univ-mrs.fr or milpied@ciml.univ-mrs.fr 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
