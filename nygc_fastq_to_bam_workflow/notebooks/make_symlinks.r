library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

#file_list <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/RNU12_fastq_files.r"
file_list <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/metadata/nygc_all.csv"
input_directory <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/symlinked_fastqs"
 
#all_filepaths <- read_csv(file_list, col_names = FALSE) %>% as.data.frame
all_filepaths <- read_csv(file_list, col_names = TRUE) %>% as.data.frame

colnames(all_filepaths) <- c("filepaths")

all_filepaths$basename <-  basename(all_filepaths$filepaths) #remove everything after UDN***_
all_filepaths$sampleID <- sub("_.*", "", all_filepaths$basename)
all_filepaths <- all_filepaths %>%
  mutate(type_fastq = case_when(
    grepl("R1", basename) ~ "R1",
    grepl("R2", basename) ~ "R2"
    ))
all_filepaths$new_name <- paste(all_filepaths$sampleID, "_", all_filepaths$type, ".fastq.gz", sep="")
all_filepaths$destination <- paste(input_directory, "/", all_filepaths$new_name, sep = "")  

all_filepaths_command <- paste("ln -s", all_filepaths$filepaths, all_filepaths$destination, sep = " ")
names(all_filepaths_command) <- NULL  # no column names at all

script_lines <- paste0(c("#! /bin/bash", all_filepaths_command), collapse = "\n")

#output_file <- paste(input_directory, "/", "pileups_symlinking_instructions.sh", sep="")
output_file <- paste(input_directory, "/", "pileups_symlinking_instructions_no_4_og_files.sh", sep="")

write(script_lines, output_file)
