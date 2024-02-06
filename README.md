# ChIPseq-Pipeline
Standard ChIPseq pipeline in the Turck lab

The ChIP-seq peipeline uses `BWA` for mapping `Fastq` reads to the `Arabidopsis thaliana` genome version `TAIR10`, `EPIC2` to identifiy peaks; `IDR` to overlap replicates and `HOMER` to annotate the peaks to genes.

# Before starting the analysis

## Install micromamba

Install `micromamba` (or mamba, conda, anaconda) on a linux based server following instructions provided on the micromamba website <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>.

## Clone this GitHub repository

Move to the location where you want to set up your working_directory and type

``` bash
gh repo clone Turckie/ChIPseq-Pipeline
```

Preserve the folder structure and call all scripts from the working directory to preserve file path definitions.

## Install all programs and dependencies

Install everything required for the analysis in a specialized `micromamba` environment by running the commands below. It installs an environment called Turck_ChIPseq that contains all programs needed for the analysis.

If the installation goes smoothly, the analyses should be run through once the environment is activated.

``` bash
micromamba activate
micromamba create -n Turck_ChIPseq -f Turck_ChIPseq_env.yml
micromamba activate Turck_ChIPseq
```

# EPIC and IDR based ChIP seq analysis

Scripts starting with `EPIC-ChIP-seq_N*` can be run consecutively in the order indicated by the number to recapitulate the analysis of Krause et al. Alternatively, it would also be possible to create a custom `bash` script to run one block after the other.

The pipeline includes the generation of a custom `blacklist` of regions that are over or under-represented in the mapping and should therefore be excluded.

The annotation is based on the generation of a custom version of the `Araport11` genome annotation that is not excluded in the `HOMER` suite.


1.  **Create a custom annotation for Araport11 and place it in the correct folder in HOMER**

The `Homer` depository of genome annotations only contains the older `TAIR10` annotation, not the more recent `Araport11` annotation. Furthermore, the default annotation is not quite suited for the small *Arabidopsis thaliana* genome since the TSS and TTS region are extended by 1000bp, which often includes the entire intergenic region. It is better to use a more fine-grained approach with Araport11. With the following helper script, Araport11 annotation is added to the Homer repository with a custom annotation of promoters.

With this custom Araport11 annotation, the strategy of assigning genes will be in this order:

1.  -100 to +100 around the TSS as promoter-TSS
2.  genebody (exon or intron)
3.  -100 to +100 at the TTS
4.  1kb promoter
5.  3kb promoter.
6.  Everything else will be assigned as intergenic region.

Provided that you have first installed the TRB_Krause_et_al environment in `micromamba`, this should get everything in place if launched from your working directory:

``` bash
./scripts/Custom_annotation_Homer.bash
```

2. **Get the read data and create a sample table**
*In case of published data:*
With the study name on hand, the EBI website <https://www.ebi.ac.uk/ena/browser/view/PRJNA329443> the column selector allows the creation of a text file following the example below. Replace the study name with the relevant study. 
The first columns are downloaded as `.csv` file from EBI, the `sample title` needed a bit of by-hand reformatting, and the last column with `type` is added by hand. Safe the table as Linux-formated (end-of-line) text file. The file `Samples-ChIPseq.txt` contains an example. 

```text
study_accession	sample_accession	run_accession	fastq_ftp	sra_md5	sample_title	type

```
*In case of data provided from a sequencing centre*
Data from the Max Planck Genome Centre Cologne are downloaded from a WebSafe server. Prepare a `.txt` file with the following headers (the file `Samples.txt` contains an example:

```text
link	file	md5sum	sample	type QC

```


3. **Download the read data**
The data downloaded from the `ftp` link provided in the `Samples_ChIPseq.txt` will be saved in `data/fastq`. The script also checks if the download was complete using the provided md5 checksums.

``` shell
#for published data
./scripts/EPIC_ChIP-seq_1.get.reads.bash Samples.ChIP-seq.txt
#for genome center data
./scripts/EPIC_ChIP-seq_1c.get.reads.wget.bash Samples.txt
```
Alternatively to all this scripting, the data can be downloaded "by Hand" and stored in the directory `data/fastq`. This can be done via a browser or by using the terminal:

``` shell
#first create the folder in which the data will be placed

working_directory=$(pwd)
path_reads=${working_directory}/data/fastq
mkdir -p ${path_reads} 

#use wget to download data into a specificed folder
wget -P ${path_reads} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR216/001/ERR2167261/ERR2167261.fastq.gz #example for one link
```


4. **Align reads to genome**
The script expects to find the `Samples_ChIPseq.txt` folder and the base name of the reads under the header `run_accession`. It uses `bwa` for read alignment. If no index for the genome is found, the script will generate the index. The alignments will end up in a `data/alignmentsChIP-seq`. The reads are automatically sorted and indexed. A second version of alignments filtered to contain only uniquely mapping reads is also stored.

```bash
./scripts/EPIC_ChIP-seq_2.alignment.bash #aligns ChIP-seq reads
echo "alignment completed"
#if the script is submitted using the lsf queue, use this command
#bsub -q multicore20 -n 8 -R "span[hosts=1] rusage[mem=40000]" -M 48000 ./scripts/EPIC_ChIPseq_2.alignment.bash Samples_ChIPseq.txt


```
Before continuing with the rest of the analysis it is good practice to have a look at the alignments, e.g. in the IGV browser. 
* Can some peaks, e.g. at already known target genes be detected?
* Does the background look even, with fluctuations that look slightly like nucleosome profiles or *spikey*?
* Do the controls/Inputs look more evenly distributed than the ChIP samples or do they have even more spikes?
* Are there tracks that look like outliers?

5. **Pool controls and inputs**
Usually the replicates of ChIP, input and control data are not linked. Since the controls can sometimes also show problems, such as random spikes, it is a better strategy to pool the controls
6. 


2.  **ChIP-seq analysis**

This should run through, but it is probably better to write another shell script that can be submitted to a queue. Make sure that the Turck_ChIPseq environment is activated in the submission system and that all scripts are executable to the system (depending on the system set-up).

``` shell


./scripts/EPIC_ChIP-seq_3.control_pooling.bash # creates pooled control samples since there is no obvious pairing between control and samples. Expects to find randomsplitbam.py in the scripts folder
echo "control samples are pooled"
./scripts/
./scripts/EPIC_ChIP-seq_4.make.blacklist.bash #Creates a blacklist of regions that are over- or underrepresented in reads. Expects to find blacklist.r in the scripts folder
echo "final blacklist generated".
./scripts/EPIC_ChIP-seq_5.epic2-IDR.bash #Binding site enrichment analysis with epic2 combined with replicate pooling using IDR
echo "bindig site enrichment completed"
./scripts/EPIC_ChIP-seq_6.PseudoIDR.bash #as above but on Psuedoreplicates and subsamples as quality check
echo "Pseudoreplication quality check completed"
./scripts/EPIC_ChIP-seq_7.Homer_Annotation.bash #annotates the target sites to genes and location
echo "binding sites annotated to genes"
./scripts/EPIC_ChIP-seq_8.makebigwig.bash #uses deeptools to generate CPM normalized bigwig coverage files, expects "Samples.ChIP-seq.txt" in the working directory
echo "Bigwig files for visualization generated"
./scripts/EPIC_ChIP-seq_9.Heatmap.bash #using deeptools, draws diverse heatmaps to visualize the data
echo "Heatmaps drawn"
./scripts/EPIC_ChIP-seq_10.makeGoodShuffle.bash # creates a shuffled control data set that keeps the relative frequency of peaks to genes as the original dataset.
echo "good shuffle set created"
./scripts/EPIC_ChIP-seq_11.peakcategories.bash #uses a merged peak dataset to annotate many features in a complete data table. Expects to find Peak_categories.r in the scirpts folder
echo "final master table for analysis generated"
```
