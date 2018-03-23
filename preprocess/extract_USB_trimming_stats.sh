#!/bin/bash

# call: ./extract_USB_trimming_stats.sh ./trimming max_length

dir=$1
ml=$3


#collect quality metrics
stats=$dir/"USB_trimming_stats.tsv"
echo -e "sample\treads\tfrac_adapter\tremoved_shorter_16\tremoved_longer_${ml}\tfrac_quality\treads_trimmed" > $stats

for f in $dir/*fastq.gz_trimming_report.txt; do
  sample=$( basename $( grep "Input filename" $f | cut -d " " -f3 ) | cut -d "_" -f1 )
  rps=$( grep "Total reads processed" $f | rev | cut -d " " -f1 | rev )
  far=$( grep "Reads with adapters" $f | rev | cut -d " " -f1 | rev | tr -d '()' )
  rm16=$( grep "Sequences removed because they became shorter than the length cutoff of 16 bp" $f | cut -f2 | cut -d " " -f1 )
  rm49=$( grep "Sequences removed because after trimming they were longer than the maximum length cutoff of $ml bp" $f | cut -f2 | cut -d " " -f1 )
  # open second report
  f2=${f/.fastq.gz/_trimmed.fq.gz}
  fqr=$( grep "Quality-trimmed: " $f2 | rev | cut -d " " -f1 | rev | tr -d '()' )
  #rm1=$( grep "Sequences removed because they became shorter than the length cutoff of 16 bp" $f2 | cut -f2 | cut -d " " -f1 )
  proc2=$( grep "sequences processed in total" $f2 | cut -d " " -f1 )
  reads_trimmed=$((proc2-rm1))
  printf "%s\t" $sample $rps $far >> $stats
  printf "%'d\t" $rm16 $rm49 >> $stats
  printf "%s\t" $fqr >> $stats
  printf "%'d\n" $reads_trimmed >> $stats 
done 
