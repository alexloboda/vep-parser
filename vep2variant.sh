#!/bin/bash

# Convert vep.tsv to file containing list of variant in the format 
# chr:pos\tref\talt

# Usage: bash vep2variant.sh vep.tsv > variant.txt

awk -F"\t" 'NR>1 {print "chr"$1":"$2"\t"$3"\t"$4}' $1