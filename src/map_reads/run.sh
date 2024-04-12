#!/bin/bash

nextflow run nf-core/rnaseq \
  --input samplesheet.csv \
  --fasta /home/biostats_share/mancchri/genomes/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz \
  --gtf /home/biostats_share/mancchri/genomes/Sus_scrofa.Sscrofa11.1.108.gtf.gz \
  --pseudo_aligner salmon \
  --skip_markduplicates \
  --skip_bigwig \
  --skip_qualimap \
  --skip_stringtie \
  --extra_salmon_quant_args '--seqBias --gcBias' \
  --outdir . \
  -c config.conf \
  -resume \
  -bg
   


