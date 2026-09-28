library(readxl)
library(magrittr)
require(data.table)
library(R.utils)  
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

txt_fp <- args[1]
txt_file <- fread(txt_fp, header=FALSE) 
output_file_fp <- args[2]

status <- txt_file %>%
  mutate(
    Unaligned = case_when(
      V4 == "PASSED_UNALIGN" ~ 1,
      V4 == "FAILED_UNALIGN" ~ 0
    ),
    Cutadapt = case_when(
      V4 == "FAILED_CUTADAPT" ~ 0,
      V4 == "PASSED_CUTADAPT" ~ 1,
    ),
    Shuffle = case_when(
      V4 == "FAILED_SHUFFLE" ~ 0,
      V4 == "PASSED_SHUFFLE" ~ 1,
    ),
    STAR_Align = case_when(
      V4 == "FAILED_ALIGN" ~ 0,
      V4 == "PASSED_ALIGN" ~ 1,
    )
  )

  write_csv(status, output_file_fp)
