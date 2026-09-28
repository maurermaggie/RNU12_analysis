#isolates the header 
samtools view -H /oak/stanford/groups/euan/UDN/gateway/data/Other_sites/Harvard/RNAseq/UDN029422/UDN029422-3247057-MGML0110-BLR1-R1.sorted.md.cram > /oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/cram_to_fastq_workflow/header_nonsense/header.sam

#edits the header
sed -i 's#/mnt/share/database/refgenome/refgenome-hg38_rnaseq_gtex_v10/fasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta#/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/reference/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta#g' /oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/cram_to_fastq_workflow/header_nonsense/header.sam 
#this edits in place ^

#cp the bam file cause of course I can't access the file and need to do all of this on a copy
cp /oak/stanford/groups/euan/UDN/gateway/data/Other_sites/Harvard/RNAseq/UDN029422/UDN029422-3247057-MGML0110-BLR1-R1.sorted.md.cram /oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/converted_fastq_from_cram/copied_cram/UDN029422.cram
samtools reheader -i /oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/cram_to_fastq_workflow/header_nonsense/header.sam /oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/converted_fastq_from_cram/copied_cram/UDN029422.cram

        samtools view -@ 8 -b -T "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/reference/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta" "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/converted_fastq_from_cram/copied_cram_UR_reheadered/UDN029422.cram" | \
        samtools collate -@ 8 -O -u {input.input_file} | \
        samtools fastq -@ 8 \
            -1 {output.fastq_R1} \
            -2 {output.fastq_R2}