#!/bin/bash
READ_ID="@NS500418:485:HC7L5BGXY:2:21202:6534:10233"
FILE_R1="/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/fastq/RD046.R1.fastq.gz"
FILE_R2="/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/fastq/RD046.R2.fastq.gz"

# Function to check a single FASTQ file for a given read ID
check_fastq_record() {
  local file=$1
  local id=$2
  zcat "$file" | awk -v id="$id" '
    BEGIN { found=0 }
    {
      if (found==1 && NR==record_line+1) { seq=$0; seq_len=length(seq) }
      if (found==2 && NR==record_line+3) { qual=$0; qual_len=length(qual) }
      if ($0 == id) { found=1; record_line=NR; found_seq=1 }
      else if (found==1 && $0 == "+") { found=2 }
    }
    END {
      if (seq_len == "" || qual_len == "") {
        print "Record not found or incomplete for", id, "in", file
      } else {
        if (seq_len == qual_len) {
          print "OK:", id, "in", file, "-> seq_len=", seq_len, "qual_len=", qual_len
        } else {
          print "MISMATCH:", id, "in", file, "-> seq_len=", seq_len, "qual_len=", qual_len
        }
      }
    }
  '
}

check_fastq_record "$FILE_R1" "$READ_ID"
check_fastq_record "$FILE_R2" "$READ_ID"