library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)

args <- commandArgs(TRUE)
file_list <- args[1]
input_directory <- args[2]
 
all_filepaths <- read_csv(file_list, col_names = TRUE) %>% as.data.frame
colnames(all_filepaths) <- c("filepaths")

all_filepaths$basename <-  basename(all_filepaths$filepaths) 
all_filepaths$remove_ending <- gsub("-Z.*", "", all_filepaths$basename)
all_filepaths$new_basename <- sub("^[^-]+-", "", all_filepaths$remove_ending) %>% paste(".cram", sep="")

all_filepaths$destination <- paste(input_directory, "/", all_filepaths$new_basename, sep = "")  

all_filepaths_command <- paste("ln -s", all_filepaths$filepaths, all_filepaths$destination, sep = " ")
names(all_filepaths_command) <- NULL  # no column names at all

script_lines <- paste0(c("#! /bin/bash", all_filepaths_command), collapse = "\n")

output_file <- paste(input_directory, "/", "pileups_symlinking_instructions.sh", sep="")
write(script_lines, output_file)