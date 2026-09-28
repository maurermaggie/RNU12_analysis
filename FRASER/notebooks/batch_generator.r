library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(broom)
library(rsample)
library(tidyverse)

################################################################
########-----Find if all bam and bai files created-----#########
################################################################

####-----Get list of files for nygc runs and UDN029422-----#####
nygc_samples <- list.files("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_fastqs") %>% as.data.frame() %>%
                filter(grepl(".fastq.gz", .))
colnames(nygc_samples) <- c("file_name")

nygc_samples <- nygc_samples %>% mutate(sample_id_og = sub("_R1\\.fastq\\.gz$", "", file_name)) %>% 
          mutate(sample_id = sub("_R2\\.fastq\\.gz$", "", sample_id_og)) %>% select(-sample_id_og) %>%
          select(-file_name) %>% group_by(sample_id) %>% tally %>% arrange(desc(n)) %>% select(-n) %>% ungroup

nygc_samples$batch <- "nygc_ONLY"

symlinked_bams <- list.files("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams") %>%
                as.data.frame %>% rename("file_name" = ".") %>%
                filter(grepl("\\.bam$", file_name)) %>%
                mutate(sample_id = gsub("\\.bam", "", file_name))

symlinked_bams$batch <- "nygc_and_UDN029422"

background_samples <- bind_rows(nygc_samples, symlinked_bams) %>% select(-file_name) %>% unique
background_samples <- background_samples %>% filter(!sample_id %in% c("UDN029422", "UDN112958", "UDN907359", "UDN629307", "UDN664300")) %>% unique

background_samples$type <- "non-spliceopathy"

RNU12_samples <- data.frame(sample_id=c("UDN029422", "UDN112958", "UDN907359", "UDN629307", "UDN664300"), type = c("RNU12_proband", "RNU12_proband", "RNU12_proband", "related_to_RNU12_proband", "related_to_RNU12_proband"))
RNU12_samples$batch <- c("UDN029422_only", "nygc_ONLY", "nygc_ONLY", "nygc_ONLY", "nygc_ONLY")

background_samples %>% filter(sample_id %in% RNU12_samples$sample_id)
# A tibble: 0 × 3

##########-----Join all sample type (nygc, ect.)-----###########
og_samples <- bind_rows(RNU12_samples, background_samples)

############-----Find if all bam files created-----#############
final_bams <- list.files("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47") %>% 
        as.data.frame %>% filter(grepl("\\.bam$", .)) %>%
        rename("file_name" = ".") %>%
        mutate(sample_id = gsub("\\_Aligned\\.sortedByCoord\\.out\\.bam", "", file_name)) %>% 
        select(-file_name)

intersect(og_samples$sample_id, final_bams$sample_id) %>% as.data.frame %>% arrange(desc(.))
setdiff(og_samples$sample_id, final_bams$sample_id) %>% as.data.frame %>% arrange(desc(.))
setdiff(final_bams$sample_id, og_samples$sample_id) %>% as.data.frame %>% arrange(desc(.))

############-----Find if all bai files created-----#############
final_bais <- list.files("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47") %>% 
        as.data.frame %>% filter(grepl("\\.bai$", .)) %>%
        rename("file_name" = ".") %>%
        mutate(sample_id = gsub("\\_Aligned\\.sortedByCoord\\.out\\.bam\\.bai", "", file_name)) %>% 
        select(-file_name)

intersect(og_samples$sample_id, final_bais$sample_id) %>% as.data.frame %>% arrange(desc(.))
setdiff(og_samples$sample_id, final_bais$sample_id) %>% as.data.frame %>% arrange(desc(.))
setdiff(final_bais$sample_id, og_samples$sample_id) %>% as.data.frame %>% arrange(desc(.))

#missing_bai <- setdiff(og_samples$sample_id, final_bais$sample_id) %>% as.data.frame %>% arrange(desc(.)) %>%
                #rename("sample_id" = ".")
#missing_bai$missing <- "Y"

#left_join(og_samples, missing_bai)
#all nygc_samples need to be run again for indexing

################################################################
###########-----Make run csvs for sample batches-----###########
################################################################
final_bams <- final_bams %>% mutate(batch = case_when(sample_id %in% nygc_samples$sample_id ~ "nygc_only",
                                        !sample_id %in% nygc_samples$sample_id ~ "nygc_and_UDN029422")) %>%
                        mutate(status = case_when(sample_id %in% c("UDN029422", "UDN112958", "UDN907359", "UDN629307", "UDN664300") ~ "RNU12_associated",
                                        ! sample_id %in% c("UDN029422", "UDN112958", "UDN907359", "UDN629307", "UDN664300") ~ "non_spliceopathy"))

#####-----Ensure none of the files are spliceopathies-----######
spliceopathies <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/Ultima_spliceopathy/metadata/intermediary/spliceopathy_cases.csv")

spliceopathies_gss_id <- spliceopathies %>% filter(!is.na(GSS_ID)) %>% pull(GSS_ID)
spliceopathies_udn_id <- spliceopathies %>% filter(!is.na(UDN_ID)) %>% pull(UDN_ID)
spliceopathies_rdid <- spliceopathies %>% filter(!is.na(RDID)) %>% pull(RDID)
spliceopathies_all_ids <- spliceopathies_gss_id %>% append(spliceopathies_udn_id) %>% append(spliceopathies_rdid)

#spliceopathies_ids <- intersect(final_bams$sample_id, spliceopathies_all_ids) %>% unique

final_bams <- final_bams %>% filter(! sample_id %in% spliceopathies_all_ids)

#######-----Ensure none of the files are watchmaker-----######
watchmaker <- read_excel("/oak/stanford/groups/smontgom/shared/UDN/Analysis/Ultima_spliceopathy/metadata/intermediary/Watchmaker.xlsx", col_names=TRUE) %>% tail(-1)
colnames(watchmaker) <- c("sampleID", "barcode")

watchmaker_ids <- intersect(final_bams$sample_id, watchmaker$sampleID)

final_bams <- final_bams %>% filter(! sample_id %in% watchmaker$sampleID)

####################-----UDN029422-----#########################
UDN029422 <- final_bams %>% filter(batch == "nygc_and_UDN029422")

UDN029422 <- data.frame(sampleID = UDN029422$sample_id, 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", UDN029422$sample_id, "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN029422 %>% filter(sampleID == "UDN029422")

write_csv(UDN029422, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN029422_non_nygc.csv")

####################-----UDN112958 nygc-----####################
UDN112958 <- final_bams %>% filter(!sample_id %in% c("UDN029422", "UDN907359", "UDN629307", "UDN664300"))

UDN112958 <- data.frame(sampleID = UDN112958$sample_id, 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", UDN112958$sample_id, "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN112958_ind <- data.frame(sampleID = "UDN112958", 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", "UDN029422", "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN112958 <- bind_rows(UDN112958,UDN112958_ind)

UDN112958 %>% filter(sampleID == "UDN112958")

write_csv(UDN112958, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN112958_nygc.csv")

####################-----UDN907359 nygc-----####################
UDN907359 <- final_bams %>% filter(!sample_id %in% c("UDN029422", "UDN112958", "UDN629307", "UDN664300"))

UDN907359 <- data.frame(sampleID = UDN907359$sample_id, 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", UDN907359$sample_id, "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN907359_ind <- data.frame(sampleID = "UDN907359", 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", "UDN907359", "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN907359 <- bind_rows(UDN907359,UDN907359_ind)

UDN907359 %>% filter(sampleID == "UDN907359")

write_csv(UDN907359, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN907359_nygc.csv")

####################-----UDN629307 nygc-----####################
UDN629307 <- final_bams %>% filter(!sample_id %in% c("UDN029422", "UDN112958", "UDN907359", "UDN664300"))

UDN629307 <- data.frame(sampleID = UDN629307$sample_id, 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", UDN629307$sample_id, "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN629307_ind <- data.frame(sampleID = "UDN629307", 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", "UDN629307", "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN629307 <- bind_rows(UDN629307,UDN629307_ind)

UDN629307 %>% filter(sampleID == "UDN629307")

write_csv(UDN629307, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN629307_nygc.csv")

####################-----UDN664300 nygc-----####################
UDN664300 <- final_bams %>% filter(!sample_id %in% c("UDN029422", "UDN112958", "UDN907359", "UDN629307"))

UDN664300 <- data.frame(sampleID = UDN664300$sample_id, 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", UDN664300$sample_id, "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN664300_ind <- data.frame(sampleID = "UDN664300", 
                bamFile = paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/star_output_trimmed_hg38_v47/", "UDN664300", "_Aligned.sortedByCoord.out.bam", sep=""),
                pairedEnd=TRUE) %>% as.tibble()

UDN664300 <- bind_rows(UDN664300,UDN664300_ind)

UDN664300 %>% filter(sampleID == "UDN664300")

write_csv(UDN664300, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN664300_nygc.csv")