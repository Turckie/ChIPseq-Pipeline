#!/bin/bash

#############################################################################################################################################################
#ChIP-seq will compare two replicates of IP samples for TRB1/2/3-GFP to two controls that were carried out by precipitating
#an untransformed wild-type Col-0 plant with the same anti-GFP polyconal antibody. We find that this is an important control 
#since we have realized that there is some `chromatin-stickness` that is particular to the specific lot of an antibody preparation. Since there is no obvious pairing
#of the controls to the replicates, it makes sense to create a pool of the libraries and compare each ChIP-seq library to this pooled
#control. We also have input samples, which can be used as control. Since they are not paired to a specific sample, they will also be subsampled and merged.
#############################################################################################################################################################

working_directory=$(pwd)
path_alignments=${working_directory}/data/alignments_ChIP-seq/unique
control=(ERR11556839 ERR11556840 ERR833053 ERR833054)
input=(ERR11556845 ERR11556846 ERR11556847)

cd ${path_alignments}

for i in ${control[@]}
do
samtools view -bhs 42.25 ${i}.uniq.bam -o ${i}.control.subsample.bam #using 42 as random number initiator this subsamples 25% of the reads
done

samtools cat *.control.subsample.bam | samtools sort - -o controls.uniq.bam
samtools index controls.uniq.bam

for i in ${input[@]}
do
samtools view -bhs 42.33 ${i}.uniq.bam -o ${i}.input.subsample.bam #using 42 as random number initiator this subsamples 33% of the reads
done

samtools cat *.input.subsample.bam | samtools sort - -o inputs.uniq.bam
samtools index inputs.uniq.bam


#clean up
rm *.subsample.bam

cd ${working_directory}
