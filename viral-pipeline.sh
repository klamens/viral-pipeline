#!/bin/bash
INPUT="./fastq_directories_example.csv";
WORKDIR="workDirectory";
STORAGEDIR="storageDirectory";
perl ./GetUnmappedReads.pl $INPUT $STORAGEDIR $WORKDIR;
perl ./map_virusses_frhit_v1.4.pl ${STORAGEDIR}/1st_unmapped $WORKDIR viral-genomes.fasta;
perl ./make_fr-hit_count_report_paired_v1.4.pl $WORKDIR ${WORKDIR}/fr-hit_results ./viral_ebi-ncbi_12122013_unique.csv;