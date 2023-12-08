# ChIPseq-Pipeline
Standard ChIPseq pipeline in the Turck lab

The ChIP-seq peipeline uses `BWA` for mapping `Fastq` reads to the `Arabidopsis thaliana` genome version `TAIR10`, `EPIC2` to identifiy peaks; `IDR` to overlap replicates and `HOMER` to annotate the peaks to genes.

# Before starting the analysis

## Install micromamba

Install `micromamba` (or mamba, conda, anaconda) on a linux based server following instructions provided on the micromamba website <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>.

## Clone this GitHub repository

Move to the location that you want to set up your working_directory and type

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

Scripts starting with `EPIC-ChIP-seq_N*` can be run consecutively in the order indicated by the number to recapitulate the analysis of Krause et al. Alternatively, it would also be possible to create a `snakefile` or a custom `bash` script to run one block after the other.

The pipeline includes the generation of a custom `blacklist` of regions that are over or under-represented in the mapping and should therefore be excluded.

The annotation is based on the generation of a custom version of the `Araport11` genome annotation that is not excluded in the `HOMER` suite.


1.  **Create a custom annotation for Araport11 and place it in the correct folder in HOMER**

The `Homer` depository of genome annotations only contains the older `TAIR10` annotation, not the more recent `Arapor11` annotation. Furthermore, the default annotation is not quite suited for the small *Arabidopsis thaliana* genome since the TSS and TTS region are extended by 1000bp, which often includes the entire intergenic region. It is better to use a more fine-grained approach with Araport11. With the following helper script, Araport11 annotation is added to the Homer repository with a custom annotation of promoters.

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

With the study name on hand, the EBI website <https://www.ebi.ac.uk/ena/browser/view/PRJNA329443> the column selector allows the creation of a text file following the example below. Replace the study name with the relevant study. 
The first columns are downloaded as `.csv` file from EBI, the `sample title` neede a bit fo by hand reformating and a last column with `Type` was added by hand. 

```text
study_accession	sample_accession	run_accession	fastq_ftp	sra_md5	sample_title	Type
PRJNA329443	SAMN05412815	SRR3928027	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/007/SRR3928027/SRR3928027.fastq.gz	b1040c4f369409b43a93bcfd56d100d8	H3-ChIP-Col-bioRep2	H3
PRJNA329443	SAMN05412819	SRR3928031	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/001/SRR3928031/SRR3928031.fastq.gz	ee4c505db7944472bcb856a8a9592f40	H3K27me3-ChIP-Col-bioRep3	H3K27me3
PRJNA329443	SAMN05412820	SRR3928032	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/002/SRR3928032/SRR3928032.fastq.gz	d98fe7b0ec6d62fb58fbe3945f16336c	Input-C24-bioRep-1	Input-C24
PRJNA329443	SAMN05412821	SRR3928033	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/003/SRR3928033/SRR3928033.fastq.gz	b15b44a3d3f8f02d3e507f2fd493033b	Input-chromatin-C24-bioRep2	Input-C24
PRJNA329443	SAMN05412823	SRR3928035	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/005/SRR3928035/SRR3928035.fastq.gz	1e8b5476f61a1f43f51ef72a260ec773	FIE-ChIP-C24-bioRep1	FIE
PRJNA329443	SAMN05412824	SRR3928036	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/006/SRR3928036/SRR3928036.fastq.gz	f9efcbb361d195612ac51a88e43c06cd	FIE-ChIP-C24-bioRep2	FIE
PRJNA329443	SAMN05412825	SRR3928037	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/007/SRR3928037/SRR3928037.fastq.gz	3f75ea4ba64f6abe9beb23644fbfe3f3	FIE-ChIP-C24-bioRep3	FIE
PRJNA329443	SAMN05412831	SRR3928043	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/003/SRR3928043/SRR3928043.fastq.gz	cb22caac4d2ab22082712cc19b8a7c88	AZF1-ChIP-Col-bioRep3	AZF
PRJNA329443	SAMN05412814	SRR3928026	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/006/SRR3928026/SRR3928026.fastq.gz	ae926f0df55cf4008cab7769aeb1354f	H3-ChIP-Col-bioRep1	H3
PRJNA329443	SAMN05412817	SRR3928029	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/009/SRR3928029/SRR3928029.fastq.gz	aad35a8b51bf22ad30046cb813674eea	H3K27me3-ChIP-Col-bioRep1	H3K27me3
PRJNA329443	SAMN05412828	SRR3928040	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/000/SRR3928040/SRR3928040.fastq.gz	5e5bd58e148ed2facf6c435ced8755e6	Input-Col-bioRep3	Input
PRJNA329443	SAMN05412830	SRR3928042	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/002/SRR3928042/SRR3928042.fastq.gz	75df4b72e40858943740eefe8d37b8f3	AZF1-ChIP-Col-bioRep2	AZF1
PRJNA329443	SAMN05412834	SRR3928046	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/006/SRR3928046/SRR3928046.fastq.gz	b2016569c55dc2d36d083b7bdf8f26ad	BPC1-ChIP-Col-bioRep3	BPC1
PRJNA329443	SAMN05412816	SRR3928028	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/008/SRR3928028/SRR3928028.fastq.gz	d8593d7239217c4ffe61c1a7e1712eaf	H3-ChIP-Col-bioRep3	H3
PRJNA329443	SAMN05412818	SRR3928030	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/000/SRR3928030/SRR3928030.fastq.gz	eb6cbd916c0fc6d3925570429451f93a	H3K27me3-ChIP-Col-bioRep2	H3K27me3
PRJNA329443	SAMN05412822	SRR3928034	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/004/SRR3928034/SRR3928034.fastq.gz	dbfe39ea87385667aae249df5de30a98	Input-C24-bioRep3	Input-C24
PRJNA329443	SAMN05412826	SRR3928038	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/008/SRR3928038/SRR3928038.fastq.gz	7ca717afa2614a66fb334986510675a9	Input-Col-bioRep1	Input
PRJNA329443	SAMN05412827	SRR3928039	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/009/SRR3928039/SRR3928039.fastq.gz	5ccd3c657f6fcfa17bf11513cb201530	Input-Col-bioRep2	Input
PRJNA329443	SAMN05412829	SRR3928041	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/001/SRR3928041/SRR3928041.fastq.gz	df970632445c0d3f04a700fe7081bf1c	AZF1-ChIP-Col-bioRep1	AZF1
PRJNA329443	SAMN05412832	SRR3928044	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/004/SRR3928044/SRR3928044.fastq.gz	0022af5eebed98cf0a0f75ebc835be0b	BPC1-ChIP-Col-bioRep1	BPC1
PRJNA329443	SAMN05412833	SRR3928045	ftp.sra.ebi.ac.uk/vol1/fastq/SRR392/005/SRR3928045/SRR3928045.fastq.gz	e5643cfa69c692af3d03ef9a92e20cff	BPC1-ChIP-Col-bioRep2	BPC1
```

3. **Download the read data**
The data downloaded from the `ftp` link provided in the `Samples_ChIPseq.txt` will be saved in `data/fastq`. The script also checks if the download was complete using the provided md5 checksums.

``` shell
./scripts/EPIC_ChIP-seq_1.get.reads.bash Samples.ChIP-seq.txt #downloads reads and annotations, expects "Samples.ChIP-seq.txt" in the working directory
echo "download completed"
```

5. 


2.  **ChIP-seq analysis**

This should run through, but it is probably better to write another shell script that can be submitted to a queue. Make sure that the Turck_ChIPseq environment is activated in the submission system and that all scripts are executable to the system (depending on the system set-up).

``` shell
./scripts/EPIC_ChIP-seq_1.get.reads.bash Samples.ChIP-seq.txt #downloads reads and annotations, expects "Samples.ChIP-seq.txt" in the working directory
echo "download completed"
./scripts/EPIC_ChIP-seq_2.alignment.bash #aligns ChIP-seq reads
echo "alignment completed"
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
