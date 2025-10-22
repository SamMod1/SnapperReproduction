#!/bin/bash

# Note: this is an example script. The this code cannot run without access to the genomic data which cannot be shared publically.
# The samplesheets can be found in the samplesheets directory.

OUTDIR=output
#mkdir -p output

SAMPLESHEET=/workspace/cfnsjm/my_files/C_auratus/transcriptomics/rnaseq/samplesheets/gonad_samplesheet.csv
INPUT_GENOME=/workspace/cfnsjm/my_files/C_auratus/transcriptomics/rnaseq/ChrAur2_renamed_fixed.fasta
GENOME_ANNOTATION=/workspace/cfnsjm/my_files/C_auratus/transcriptomics/rnaseq/ChrAurV1_fixed.gff

module load nextflow/24.04.4

nextflow pull nf-core/rnaseq
nextflow run nf-core/rnaseq --input $SAMPLESHEET --outdir $OUTDIR --fasta $INPUT_GENOME --gff $GENOME_ANNOTATION --extra_trimgalore_args="--clip_r1=1 --clip_r2=1" --skip_biotype_qc -profile singularity