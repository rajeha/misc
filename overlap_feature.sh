#!/usr/bin/env bash
set -e

usage='overlap_feature.sh <svs.bedpe> <features.bed>'
# features.bed must only have 3 fields

if [ $# -ne 2 ]; then
	echo $usage	
	exit 1
fi

# 1) PROCESS EVERYTHING BUT TRANSLOCATIONS
cat $1 | grep -v TRA | cut -f1,3,5,7 | perl -nale 'if ($F[1]>$F[2]){print join "\t", @F[0,2,1,3]} else {print}' | sort -k1,1 -k2,2n -k3,3n > nontra.oftmp

bedtools closest -t first -D ref -a nontra.oftmp -b $2 | awk '$8==0' | cut -f4 > nontra_ids.oftmp

LC_ALL=C grep -w -F -f nontra_ids.oftmp $1  > nontra_overlap.oftmp

# 2) PROCESS TRANSLOCATIONS
# 2.1) Find overlaps with first end. 
cat $1 | grep TRA | cut -f1,2,3,7 | perl -nale 'if ($F[1]>$F[2]){print join "\t", @F[0,2,1,3]} else {print}' | sort -k1,1 -k2,2n -k3,3n > tra1.oftmp

bedtools closest -t first -D ref -a tra1.oftmp -b $2 | awk '$8==0' | cut -f4 > tra1_ids.oftmp

# 2.2) Find overlaps with second end.
cat $1 | grep TRA | cut -f4,5,6,7 | perl -nale 'if ($F[1]>$F[2]){print join "\t", @F[0,2,1,3]} else {print}'| sort -k1,1 -k2,2n -k3,3n > tra2.oftmp

bedtools closest -t first -D ref -a tra2.oftmp -b $2 | awk '$8==0' | cut -f4 > tra2_ids.oftmp

cat tra1_ids.oftmp tra2_ids.oftmp | sort | uniq > tra_ids.oftmp
 
LC_ALL=C grep -w -F -f tra_ids.oftmp $1  > tra_overlap.oftmp || true

cat nontra_overlap.oftmp tra_overlap.oftmp

rm -f *.oftmp


