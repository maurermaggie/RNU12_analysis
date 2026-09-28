library(readxl)
require(data.table)
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(ComplexUpset)
library(broom)
library(rsample)
library(tidyverse)

check_ids <- function(df, lst) {
    lst_join <- paste(lst, collapse = "|")
    in_udn <- df %>% filter(grepl(lst_join, UDN_ID)) %>% pull(ID)
    in_gss <- df %>% filter(grepl(lst_join, GSS_ID)) %>% pull(ID)
    in_rd <- df %>% filter(grepl(lst_join, RDID)) %>% pull(ID)

    in_all <- in_udn %>% append(in_gss) %>% append(in_rd)

    df <- df %>% filter(!ID %in% in_all)
    df   
}

################################################################
####################-----Metadata load-----#####################
################################################################

####################-----ID information-----####################
#ids of files from the same batch as the four fastqs
nygc_polyA <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/metadata/nygc_polyA_batch.csv")

#list of sampleIDs assoicated with files processed via watchmaker (at Stanford)
watchmaker <- read_excel("/oak/stanford/groups/smontgom/shared/UDN/Analysis/Ultima_spliceopathy/metadata/intermediary/Watchmaker.xlsx", col_names=TRUE) %>% tail(-1)
colnames(watchmaker) <- c("sampleID", "barcode")

#information (including sampleIDs) of individuals with spliceopathies + their family
spliceopathies <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/Ultima_spliceopathy/metadata/intermediary/spliceopathy_cases.csv")

##################-----Actual filepaths-----###################

#nygc_files
nygc_files <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/metadata/nygc_files.csv", col_names=FALSE)
colnames(nygc_files) <- c("rel_path")
nygc_files$rel_path <- sub("\\.", "", nygc_files$rel_path)
nygc_files$full_path <- paste("/oak/stanford/scg/prj_gregor/transfer/nygc/Project_SMI_20033_B01_NAN_Lane.2025-04-04", nygc_files$rel_path, sep="")
nygc_files <- nygc_files %>% filter(grepl("GSS|UDN", full_path)) %>% select(-rel_path)
write_csv(nygc_files, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/metadata/nygc_all.csv")

#filepaths (no preprocessing) from GREGoR U12
#GREGoR_bams <- list.files("/oak/stanford/scg/prj_gregor/consortium-lr-rna/input_bams") %>% as.data.frame %>% filter(grepl("UDN|RD|GSS", .))
#colnames(GREGoR_bams) <- c("sampleID")
#GREGoR_bams$filepath <- paste("/oak/stanford/scg/prj_gregor/consortium-lr-rna/input_bams/", GREGoR_bams$sampleID, sep= "")

#filepaths of files used in the last FRASER run of my phenotype paper (plus some??)
UDN <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/FRASER/FRASER_snakemake_cleaned_github/config/RNU4ATAC_paper_rm_12_6.csv", col_names=FALSE)
colnames(UDN) <- c("filepath")

######################-----Metadata-----######################
metadata <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/metadata/FRASER_metadata_updated_batches.csv")
#GSS_meta <- read_csv("/oak/stanford/scg/prj_gregor/consortium-lr-rna/participant_information_mar_8.csv")

################################################################
###############-----Remove 'bad' samples'-----##################
################################################################

###############-----Remove all watchmaker-----##################
#this must be done with Tecan
gss_id <- metadata %>% filter(!is.na(GSS_ID)) %>% pull(GSS_ID)
udn_id <- metadata %>% filter(!is.na(UDN_ID)) %>% pull(UDN_ID)
rdid <- metadata %>% filter(!is.na(RDID)) %>% pull(RDID)
all_ids <- gss_id %>% append(udn_id) %>% append(rdid)

watchmaker_ids <- intersect(all_ids, watchmaker$sampleID)

metadata_polyA <- check_ids(metadata, watchmaker_ids)

#############-----Remove all spliceopathies-----################
#the samples must not be spliceopathies (or parents of spliceopathy pts)
spliceopathies_gss_id <- spliceopathies %>% filter(!is.na(GSS_ID)) %>% pull(GSS_ID)
spliceopathies_udn_id <- spliceopathies %>% filter(!is.na(UDN_ID)) %>% pull(UDN_ID)
spliceopathies_rdid <- spliceopathies %>% filter(!is.na(RDID)) %>% pull(RDID)
spliceopathies_all_ids <- spliceopathies_gss_id %>% append(spliceopathies_udn_id) %>% append(spliceopathies_rdid)

spliceopathies_ids <- intersect(all_ids, spliceopathies_all_ids)

metadata_polyA_no_spl <- check_ids(metadata_polyA, spliceopathies_ids)

##############-----Remove w/ missing data-----##################
metadata_polyA_no_spl$RIN %>% unique
metadata_polyA_no_spl$batch %>% unique
metadata_polyA_no_spl$sex %>% unique
metadata_polyA_no_spl$age %>% unique

metadata_polyA_no_spl$affected_status %>% unique

metadata_polyA_no_spl_no_unknown <- metadata_polyA_no_spl %>% filter(affected_status != "Unknown")

##################-----Remove batch 8-----#####################
metadata_polyA_no_spl_no_unknown_no8 <- metadata_polyA_no_spl_no_unknown %>% filter(batch != "8")

########-----Remove spliceopathies from  nygc-----############
nygc_no_rnu12 <- nygc_files %>% filter(!grepl("UDN112958|UDN907359|UDN629307|UDN664300", full_path))

nygc_no_rnu12 <- nygc_no_rnu12 %>%
  mutate(sample_id = sub("_.*", "", basename(full_path)))

spl_ids <- intersect(spliceopathies_all_ids, nygc_no_rnu12$sample_id)

nygc_no_spliceopathies <- nygc_no_rnu12 %>% filter(!sample_id %in% spl_ids) #50 samples

################################################################
##############-----Stratify sample for nygc-----################
################################################################

######################-----get fold-----########################
metadata_polyA_no_spl_no_unknown_no8$RIN <- as.numeric(metadata_polyA_no_spl_no_unknown_no8$RIN)

metadata_polyA_no_spl_no_unknown_no8 <- metadata_polyA_no_spl_no_unknown_no8 %>%
  mutate(
    #age_bin = cut(age, breaks = 9),
    #RIN_bin = cut(RIN, breaks = 7),
    strata = interaction(batch, affected_status, drop = TRUE)
  )

folds <- vfold_cv(metadata_polyA_no_spl_no_unknown_no8, v = 7, strata = strata)

fold_data <- map(folds$splits, ~ assessment(.x))
fold_data <- setNames(fold_data, folds$id)  # Fold_1, Fold_2, ..., Fold_7

subsampled_data <- fold_data[[4]] %>% select(UDN_ID, GSS_ID, RDID)

###################-----Get sampleID-----#######################
UDN <- UDN %>%
  mutate(sample_id = sub("_.*", "", basename(filepath)))

###############-----Collapse subsampled data-----###############
udn_sub <- subsampled_data %>% filter(!is.na(UDN_ID)) %>% pull(UDN_ID)
gss_sub <- subsampled_data %>% filter(!is.na(GSS_ID)) %>% pull(GSS_ID)
rd_sub <- subsampled_data %>% filter(!is.na(RDID)) %>% pull(RDID)

sub_all <- udn_sub %>% append(gss_sub) %>% append(rd_sub)

lst_join <- paste(sub_all, collapse = "|")

################-----Find matching files-----###################
UDN_subsampled <- UDN %>% filter(grepl(lst_join, sample_id))

write_csv(UDN_subsampled, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/metadata/background_samples_for nygc_files.csv")

##################-----Symlink the bams-----####################
UDN_subsampled <- UDN_subsampled %>% mutate(sample_id = sub("\\..*", "", sample_id))
UDN_subsampled$destination <- paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams/", UDN_subsampled$sample_id, ".bam", sep="")
UDN_subsampled$bai <- paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams/", UDN_subsampled$sample_id, ".bam.bai", sep="")

all_filepaths_command <- paste("ln -s", UDN_subsampled$filepath, UDN_subsampled$destination, sep = " ")
all_filepaths_command_bai <- paste("ln -s", UDN_subsampled$filepath, UDN_subsampled$bai, sep = " ")

all_filepaths <- append(all_filepaths_command, all_filepaths_command_bai)

names(all_filepaths) <- NULL  # no column names at all

script_lines <- paste0(c("#! /bin/bash", all_filepaths), collapse = "\n")

output_file <- paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams/pileups_symlinking_instructions_for_nygc_control_bams.sh", sep="")

write(script_lines, output_file)

################################################################
###########-----Stratify sample for UDN029422-----##############
################################################################

######################-----get fold-----########################
subsampled_data <- fold_data[[4]] %>% bind_rows(fold_data[[5]]) %>% select(UDN_ID, GSS_ID, RDID)

###############-----Collapse subsampled data-----###############
udn_sub <- subsampled_data %>% filter(!is.na(UDN_ID)) %>% pull(UDN_ID)
gss_sub <- subsampled_data %>% filter(!is.na(GSS_ID)) %>% pull(GSS_ID)
rd_sub <- subsampled_data %>% filter(!is.na(RDID)) %>% pull(RDID)

sub_all <- udn_sub %>% append(gss_sub) %>% append(rd_sub)

lst_join <- paste(sub_all, collapse = "|")

################-----Find matching files-----###################
UDN_subsampled <- UDN %>% filter(grepl(lst_join, sample_id))

write_csv(UDN_subsampled, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/metadata/background_samples_for_UDN029422.csv")

##################-----Symlink the bams-----####################
#only need to do this for fold 5 since fold 4 was already symlinked

fold5 <- fold_data[[5]] %>% select(UDN_ID, GSS_ID, RDID)

udn_5 <- fold5 %>% filter(!is.na(UDN_ID)) %>% pull(UDN_ID)
gss_5 <- fold5 %>% filter(!is.na(GSS_ID)) %>% pull(GSS_ID)
rd_5 <- fold5 %>% filter(!is.na(RDID)) %>% pull(RDID)

ids_fold5 <- udn_5 %>% append(gss_5) %>% append(rd_5)

lst_join <- paste(ids_fold5, collapse = "|")

UDN_subsampled5 <- UDN %>% filter(grepl(lst_join, sample_id))

UDN_subsampled5 <- UDN_subsampled5 %>% mutate(sample_id = sub("\\..*", "", sample_id))
UDN_subsampled5$destination <- paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams/", UDN_subsampled5$sample_id, ".bam", sep="")
UDN_subsampled5$bai <- paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams/", UDN_subsampled5$sample_id, ".bam.bai", sep="")

all_filepaths_command <- paste("ln -s", UDN_subsampled5$filepath, UDN_subsampled5$destination, sep = " ")
all_filepaths_command_bai <- paste("ln -s", UDN_subsampled5$filepath, UDN_subsampled5$bai, sep = " ")

all_filepaths <- append(all_filepaths_command, all_filepaths_command_bai)

names(all_filepaths) <- NULL  # no column names at all

script_lines <- paste0(c("#! /bin/bash", all_filepaths), collapse = "\n")

output_file <- paste("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_bams/pileups_symlinking_instructions_for_UDN029422_control_bams.sh", sep="")

write(script_lines, output_file)
