#!/bin/bash
echo $1 $(bcftools view  -h $1 | grep '^#CHROM' | tr '\t' '\n' | wc -l)
