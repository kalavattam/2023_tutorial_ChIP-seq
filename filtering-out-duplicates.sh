#!/bin/bash

#  filtering-out-duplicates.sh
#  KA

infile="${1}"
outfile="${2}"

if [[
         -f "${infile}" \
    && ! -f "${outfile}"
]]; then
    samtools view -b -F 1024 "${infile}" -o "${outfile}"
fi

if [[ 
         -f "${outfile}" \
    && ! -f "${outfile}.bai"
]]; then
    samtools index "${outfile}"
fi

if [[ 
         -f "${outfile}" \
    &&   -f "${outfile}.bai" \
    && ! -f
]]; then