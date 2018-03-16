#!/bin/bash

# call ./USB_trimming.sh ./fastqs <fq|fastq> ./trimmed trim_galore_bin max_length

DIR=$1
fq=$2
OUT=$3
TRIM_GALORE=$4
ml=$5

# first round
for f in $DIR/*{fq}.gz; do
	echo "first round: trim $f ..."
	$TRIM_GALORE --length 16 -q 0 --trim-n --max_length $ml $f -o $OUT
done

# second round
for f in $OUT/*trimmed.fq.gz;do 
	echo "second round: trim $f ..."
	$TRIM_GALORE -q 20 --length 1 -a CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC --stringency 50 $f -o $OUT
done
