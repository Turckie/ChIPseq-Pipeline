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
