#!/bin/bash

#########################################---Step 1: Define Functions---###########################################################
log_path="${1}"

###########################################---Step 7: Find failed logs---############################################################
cd "$log_path"

if [ ! -f  "$log_path/status.txt" ]; then
    touch  "$log_path/status.txt"
fi  

if [ -d "bam_or_cram_to_fastq" ]; then
    cd "bam_or_cram_to_fastq"

    if [ ! -f "SUCCESS_bam_cram_to_fastq.txt" ]; then
        touch "SUCCESS_bam_cram_to_fastq.txt"
    fi

    grep -rl "1 of 1 steps (100%) done" >> "SUCCESS_bam_cram_to_fastq.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tPASSED_UNALIGN"
    done < "SUCCESS_bam_cram_to_fastq.txt" >> "$log_path/status.txt"

    if [ ! -f "FAIL_bam_cram_to_fastq.txt" ]; then
        touch "FAIL_bam_cram_to_fastq.txt"
    fi

    grep -rL "1 of 1 steps (100%) done" >> "FAIL_bam_cram_to_fastq.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tFAILED_UNALIGN"
    done < "FAIL_bam_cram_to_fastq.txt" >> "$log_path/status.txt"

    cd "$log_path"
fi

if [ -d "cutadapt" ]; then
    cd "cutadapt"

    if [ ! -f "SUCCESS_cutadapt.txt" ]; then
        touch "SUCCESS_cutadapt.txt"
    fi

    grep -rl "1 of 1 steps (100%) done" >> "SUCCESS_cutadapt.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tPASSED_CUTADAPT"
    done < "SUCCESS_cutadapt.txt" >> "$log_path/status.txt"

    if [ ! -f "FAIL_cutadapt.txt" ]; then
        touch "FAIL_cutadapt.txt"
    fi

    grep -rL "1 of 1 steps (100%) done" >> "FAIL_cutadapt.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tFAILED_CUTADAPT"
    done < "FAIL_cutadapt.txt" >> "$log_path/status.txt"

    cd "$log_path"
fi

if [ -d "shuffle" ]; then
    cd "shuffle"  

    if [ ! -f "SUCCESS_shuffle.txt" ]; then
        touch "SUCCESS_shuffle.txt"
    fi

    grep -rl "1 of 1 steps (100%) done" >> "SUCCESS_shuffle.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tPASSED_SHUFFLE"
    done < "SUCCESS_shuffle.txt" >> "$log_path/status.txt"

    if [ ! -f "FAIL_shuffle.txt" ]; then
        touch "FAIL_shuffle.txt"
    fi

    grep -rL "1 of 1 steps (100%) done" >> "FAIL_shuffle.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tFAILED_SHUFFLE"
    done < "FAIL_shuffle.txt" >> "$log_path/status.txt"

    cd "$log_path"
fi

if [ -d "star_align_trim" ]; then
    cd "star_align_trim"

    if [ ! -f "SUCCESS_star_align.txt" ]; then
        touch "SUCCESS_star_align.txt"
    fi

    grep -rl "1 of 1 steps (100%) done" >> "SUCCESS_star_align.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tPASSED_ALIGN"
    done < "SUCCESS_star_align.txt" >> "$log_path/status.txt"

    if [ ! -f "FAIL_star_align.txt" ]; then
        touch "FAIL_star_align.txt"
    fi

    grep -rL "1 of 1 steps (100%) done" >> "FAIL_star_align.txt"
    while IFS= read -r file; do
          cleaned=$(echo "$file" | sed 's/.*ID=//; s/,type=.*//')
          ts=$(stat -c "%y" "$file" 2>/dev/null)
          echo -e "$file\t$cleaned\t$ts\tFAILED_ALIGN"
    done < "FAIL_star_align.txt" >> "$log_path/status.txt"

fi

cd "../$(dirname "$0")"

./scripts/make_status_csv.r "$log_path/status.txt" "$log_path/overall_status.csv"
