#!/bin/bash
set -euo pipefail
# Description: This script validates that the number of features in a feature table matches the number of variants in a VCF file.
# Usage: ./validator.sh <feature_table> <vcf_file>

table=$1
vcf=$2

# count feature table rows excluding header
if [[ -z ${table} ]]; then
    >&2  echo 'ERROR: No file argument provided'
    exit 1
elif [[ -s ${table} ]]; then
    line_count=$( grep -v "^#" ${table} | wc -l )
    feature_count=$( echo ${line_count} - 1 | bc )
else
    >&2  echo 'ERROR: file does not exist or is empty. Script requires at least a file with a header line.'
    exit 1
fi
# count vcf
if [[ -e ${vcf} ]]; then
    if [[ $vcf == *.vcf.gz ]]; then
        count=$(gunzip -c ${vcf} | grep -v "^#" | grep "TYPE=SNV" | wc -l) # successfully validated but not empty
    elif [[ $vcf == *.vcf ]]; then
        count=$(grep -v "^#" ${vcf} | grep "TYPE=SNV" | wc -l) # successfully validated but not empty
    else
        >&2  echo "ERROR: VCF file in unexpected format :" ${vcf}
        echo "False"
        exit 1
    fi
else
    count=0
fi

if [[ "$feature_count" -ne "$count" ]]; then
    echo "${count}"
    echo ${feature_count}
    echo "The VCF and feature table counts are not equal."
    exit 1
fi

