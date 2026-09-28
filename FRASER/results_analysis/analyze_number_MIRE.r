library(readxl)
library(magrittr)
library(biomaRt)
library(tidyverse)
require(data.table)

############################################################################################################
###################################----Read in dataframes----###############################################
############################################################################################################

#you need the input because FRASER output doesn't include samples with 0 outliers

##########################################----FRASER1----###################################################
#UDN029422
UDN029422_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN029422_non_nygc/FRASER1/FRASER_output.csv"
UDN029422 <- read_csv(UDN029422_fp)

UDN029422_input_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN029422_non_nygc.csv"
UDN029422_input <- read_csv(UDN029422_input_fp)

#UDN112958
UDN112958_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN112958_nygc/FRASER1/FRASER_output.csv"
UDN112958 <- read_csv(UDN112958_fp)

UDN112958_input_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN112958_nygc.csv"
UDN112958_input <- read_csv(UDN112958_input_fp)

#UDN629307
UDN629307_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN629307_nygc/FRASER1/FRASER_output.csv"
UDN629307 <- read_csv(UDN629307_fp)

UDN629307_input_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN629307_nygc.csv"
UDN629307_input <- read_csv(UDN629307_input_fp)

#UDN664300
UDN664300_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN664300_nygc/FRASER1/FRASER_output.csv"
UDN664300 <- read_csv(UDN664300_fp)

UDN664300_input_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN664300_nygc.csv"
UDN664300_input <- read_csv(UDN664300_input_fp)

#UDN907359
UDN907359_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN907359_nygc/FRASER1/FRASER_output.csv"
UDN907359 <- read_csv(UDN907359_fp)

UDN907359_input_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/config/UDN907359_nygc.csv"
UDN907359_input <- read_csv(UDN907359_input_fp)

##########################################----FRASER2----###################################################
#UDN029422
UDN029422_fp2 <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN029422_non_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv"
UDN029422_2 <- read_csv(UDN029422_fp2, col_names = FALSE)
colnames(UDN029422_2) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")

#UDN112958
UDN112958_fp2 <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN112958_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv"
UDN112958_2 <- read_csv(UDN112958_fp2, col_names = FALSE)
colnames(UDN112958_2) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")

#UDN629307
UDN629307_fp2 <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN629307_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv"
UDN629307_2 <- read_csv(UDN629307_fp2, col_names = FALSE)
colnames(UDN629307_2) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")
#UDN664300
UDN664300_fp2 <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN664300_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv"
UDN664300_2 <- read_csv(UDN664300_fp2, col_names = FALSE)
colnames(UDN664300_2) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")
#UDN907359
UDN907359_fp2 <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN907359_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv"
UDN907359_2 <- read_csv(UDN907359_fp2, col_names = FALSE)
colnames(UDN907359_2) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")

############################################----MIDB----####################################################
#MIDB file used for filtering on Dec 2nd, 2025
MIDB_all_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Phenotype_Paper/data/Homo_sapiens_intron_new.csv" 
MIDB_all <- read_csv(MIDB_all_fp) 

############################################################################################################
##########################################----Functions----#################################################
############################################################################################################

#####################################----convert_to_gene_names----##########################################
#RATIONAL:
    #the fibroblast FRASER output from the Broad included ENSGs as its output in the hgncSymbol column INSTEAD of hgnc symbols/ gene names
        #so, this function converts ENSGs to hgnc symbols/ gene names
        #it is called in make_outlier_df_gene after checking if the rows in the hgncSymbol column contain the string "ENSG" > 50% of the time
#INPUT: 
    #ensg_df: the filtered FRASER table given to make_outlier_df_gene from change_colnames_bind_rows from the main_function
        #this table is sometimes altered in change_colnames_bind_rows before being passed to make_outlier_df_gene to only include introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
#OUTPUT: the same dataframe given in the INPUT with the values in the "hgncSymbol" column changed from ENSG IDs to hgnc symbols/ gene names

convert_to_gene_names <- function(ensg_df) {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    ensg_ids <- unique(gsub("\\.[0-9]+$", "", ensg_df$gene_symbol))

    # Query BioMart and map ENSG ids to hgnc/ gene symbols
    mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = ensg_ids,
    mart       = ensembl
    )

    # Remove version suffix from the returned ENSG ids so they match my input
    mapping$ensembl_gene_id <- gsub("\\.[0-9]+$", "", mapping$ensembl_gene_id)

    colnames(ensg_df) <- c("ensembl_gene_id")

    # Merge back to original data frame
    filtered_mapped <- left_join(ensg_df, mapping)
    filtered_mapped
}

####################################----make_outlier_df_gene----############################################
#INPUT: 
    #filtered_table: a filtered df created in chage_colnames_bind_rows OR a df given to main_function
        #the df given to main_function is the FRASER output filtered on the padjust, deltaPsi, and type of your choice
        #the df created in change_colnames_bind_rows is the df given to the main function after being filtered to only introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
    #input_sid: a csv file with one column (labeled sampleID) where every row contains the sampleID assocaited with the bam file
        #that was inputted in the FRASER pipeline
    #dir_path: the directory name where the output will be saved under the "unique_genes_per_sample" folder
#OUTPUT:
    #a folder (unique_genes_per_sample) containing one file per sample named [sampleID].csv with two columns
        #ensembl_gene_id: ENSG of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
            #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
            #outliers of only the specified intron type (ex: number of genes with minor intron retention)
        #hgnc_symbol: gene names of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
            #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
            #outliers of only the specified intron type (ex: number of genes with minor intron retention)

make_outlier_df_gene <- function(filtered_tables, filepath) {   
    if (!dir.exists(filepath)) {
        dir.create(filepath, recursive = TRUE)
    }

    filenumber <- length(list.files(filepath))

    ids_w_outliers <- filtered_tables$sampleID %>% unique
    ids_w_outliers_count <- filtered_tables %>% pull(sampleID) %>% unique %>% length()

    if (filenumber < ids_w_outliers_count) {
        for (i in ids_w_outliers) {
            ind_gene_df_sample <- filtered_tables %>% filter(sampleID == i)
            
            ind_gene_df <- ind_gene_df_sample %>% select(gene_symbol) %>% unique()
            
            pro_ensg <- mean(startsWith(ind_gene_df$gene_symbol, "ENSG"), na.rm = TRUE)

            if (pro_ensg > 0.5) {
                ind_gene_df <- convert_to_gene_names(ind_gene_df)
            }
            
            filepath_sample <- paste0(filepath, "/", i, ".csv")
            print(filepath_sample)
            write_csv(ind_gene_df, filepath_sample)
        }
    }
} 

######################################----get_genes_count_df----############################################
#INPUT: 
    #filtered_tables: a filtered df created in chage_colnames_bind_rows OR a df given to main_function
        #the df given to main_function is the FRASER output filtered on the padjust, deltaPsi, and type of your choice
        #the df created in change_colnames_bind_rows is the df given to the main function after being filtered to only introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
    #filepath: the filepath where you would like to save the csv output
#OUTPUT:
    #a csv with three columns:
        #sampleID: the sampleID as found in the FRASER output
        #n: the number of genes containing outliers of the specified intron type (as found in the folder name) impacted by EITHER
            #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
            #outliers of only the specified intron type (ex: number of genes with minor intron retention)
        #Z_score: the Z-score for that sample taken from the mean and sd of the n column

get_genes_count_df <- function(filtered_table, input_sid, level_name, dir_path) {
    count_output_fp <- paste0(level_name, "unique_gene_counts.csv")
    if (!file.exists(count_output_fp)) {
        genes <- filtered_table %>% select(sampleID, gene_symbol) %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
        zero_count_samples <- input_sid %>% filter(! sampleID %in% filtered_table$sampleID)
        zero_count_samples$n <- 0
        count_joined <- bind_rows(genes, zero_count_samples)
        count_joined_sdv <- count_joined %>% mutate(Z_score = (n - mean(n)) / sd(n)) %>%
            arrange(sampleID)
        write_csv(count_joined_sdv, count_output_fp)
        count_joined_sdv
    }
}

####################################----get_junctions_count_df----###########################################
#INPUT: 
    #filtered_tables: a filtered df created in chage_colnames_bind_rows OR a df given to main_function
        #the df given to main_function is the FRASER output filtered on the padjust, deltaPsi, and type of your choice
        #the df created in change_colnames_bind_rows is the df given to the main function after being filtered to only introns
            #with a start=inron_start+/-1 and an end=intron_end+/-1 where
            #start/end are from the FRASER output and
            #intron_start/intron_end are from the MIDB_database
    #filepath: the filepath where you would like to save the csv output
#OUTPUT:
    #a csv with three columns:
        #sampleID: the sampleID as found in the FRASER output
        #n: the number of junctions with outliers of the specified intron type (as found in the folder name) which are EITHER
            #outliers of any intron type (ex: number of intron retention events of any type within minor intron containing genes, or MIGs)
            #outliers of only the specified intron type (ex: number of junctions with minor intron retention)
        #Z_score: the Z-score for that sample taken from the mean and sd of the n column

get_junctions_count_df <- function(filtered_table, input_sid, level_name, dir_path) {
    count_output_fp <- paste0(level_name, "unique_junction_counts.csv")
    if (!file.exists(count_output_fp)) {
        junctions <- filtered_table %>% unique %>% group_by(sampleID) %>% tally %>% arrange(desc(n))
        zero_count_samples <- input_sid %>% filter(! sampleID %in% filtered_table$sampleID)
        zero_count_samples$n <- 0
        count_joined <- bind_rows(junctions, zero_count_samples)
        count_joined_sdv <- count_joined %>% mutate(Z_score = (n - mean(n)) / sd(n)) %>%
            arrange(sampleID)
        write_csv(count_joined_sdv, count_output_fp)
        count_joined_sdv
    } 
}

#################################----change_colnames_bind_rows----###########################################
#INPUT: 
    #intron: the intron type you want to filter the MIDB file by 
        #this will do through every intron type in the MIDB database in for loops in main_function
    #MIDB_file: the file from the MIDB database read in earlier via read_csv
    #FRASER_results_filtered: the results of FRASER given to the main function
        #this should ALREADY be filtered for the padjust, deltaPsi, and type of your choice
    #FRASER_input: a csv file with one column (labeled sampleID) where every row contains the sampleID assocaited with the bam file
        #that was inputted in the FRASER pipeline
    #dir_path: 
        #the path you want the outputs of get_genes_count_df and, if focusing on a per-gene level, make_outlier_df_gene
    #index: the number in the list given to the function by main_funcion
        #here, the odd numbers indicate that we are looking for either only:
            #genes with outliers of ONLY the specified intron type (as found in the folder name) OR
            #junctions with ONLY outliers of the specified intron type (as found in the folder name)
                #as opposed to, for example, junctions with outliers of ANY intron type in a gene containing at least one intron of 
                    #the specified intron type (as found in the folder name)
#OUTPUT:
    #if looking on a per-gene basis:
        #the output of get_genes_count_df AND make_outlier_df_gene
    #if looking on a per-junction basis:
        #the output of get_junctions_count_df

change_colnames_bind_rows <- function(intron, MIDB_file, FRASER_results_filtered, FRASER_input, level, index, dir_path) {
    pro_ensg <- mean(startsWith(FRASER_results_filtered$hgncSymbol, "ENSG"), na.rm = TRUE)

    if (pro_ensg > 0.5) {
        intron_list <- MIDB_file %>% filter(intron_class == intron) %>% pull(ensembl_gene_id) %>% unique

        names(MIDB_file)[names(MIDB_file) == "ensembl_gene_id"] <- "hgncSymbol"
        
        ####Filter FRASER results to MIGs
        FRASER_results_filtered <- FRASER_results_filtered %>% filter(hgncSymbol %in% intron_list) #this is an input for the counts function, ect.

        ####Join to intron file
        names(MIDB_file)[names(MIDB_file) == "gene_symbol"] <- "gene_name"
        MIDB_joined <- left_join(MIDB_file, FRASER_results_filtered, by="hgncSymbol")
        names(MIDB_joined)[names(MIDB_joined) == "hgncSymbol"] <- "gene_symbol"
        names(FRASER_results_filtered)[names(FRASER_results_filtered) == "hgncSymbol"] <- "gene_symbol"

    } else {
        intron_list <- MIDB_file %>% filter(intron_class == intron) %>% pull(gene_symbol) %>% unique

        ####Filter FRASER results to MIGs
        FRASER_results_filtered <- FRASER_results_filtered %>% filter(hgncSymbol %in% intron_list) #this is an input for the counts function, ect.

        ####Join to intron file
        names(FRASER_results_filtered)[names(FRASER_results_filtered) == "hgncSymbol"] <- "gene_symbol"
        MIDB_joined <- left_join(MIDB_file, FRASER_results_filtered, by="gene_symbol")
        names(FRASER_results_filtered)[names(FRASER_results_filtered) == "hgncSymbol"] <- "gene_symbol"
    }

    if (index %% 2 == 1) {
    FRASER_results_filtered <- MIDB_joined %>% filter(intron_class == "minor") %>% 
                                       filter(start == intron_start | start == intron_start-1 | start == intron_start+1 |
                                       end == intron_end|  end == intron_end-1 | end == intron_end+1) %>% 
                                       select(any_of(c("sampleID", "gene_symbol", "intron_name", "intron_start", "intron_end", "start", "end", "intron_class", "pValue", "padjust", 
                                       "zScore", "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts"))) %>% unique #this is an input for the counts function, ect.

    }

    if (grepl("genes_w", level)) {
        get_genes_count_df(FRASER_results_filtered, FRASER_input, level, dir_path)
        make_outlier_df_gene(FRASER_results_filtered, paste0(level, "unique_genes_per_sample"))
    } else {
        get_junctions_count_df(FRASER_results_filtered, FRASER_input, level, dir_path)
    }
    
}

#################################----change_colnames_bind_rows----###########################################
#INPUT: 
    #filtered_FRASER_results: the results of FRASER given to the main function
        #this should ALREADY be filtered for the padjust, deltaPsi, and type of your choice
    #FRASER_input_files: a csv file with one column (labeled sampleID) where every row contains the sampleID assocaited with the bam file
        #that was inputted in the FRASER pipeline
    #minor_intron_file: the file from the MIDB database read in earlier via read_csv
    #dir_path: the directory where you want all results saved
    #FRASER_output_fp: the filepath where filtered_FRASER_results originates from
        #this will be recorded in a readme
    #FRASER_input_fp: the filepath where FRASER_input_files originates from
        #this will be recorded in a readme
    #MIG_fp: the filepath where minor_intron_file originates from
        #this will be recorded in a readme
#OUTPUT:
    #output of get_junctions_count_df for all 6 intron types in the MIDB df (minor, minor_hybrid, major_hybrid, minor-like, major-like, and non-canonical): 
        #csv with three columns:
            #sampleID: the sampleID as found in the FRASER output
            #n: the number of junctions with outliers of the specified intron type (as found in the folder name) which are EITHER
                #outliers of any intron type (ex: number of intron retention events of any type within minor intron containing genes, or MIGs)
                #outliers of only the specified intron type (ex: number of junctions with minor intron retention)
            #Z_score: the Z-score for that sample taken from the mean and sd of the n column
    #output of get_genes_count_df for all 6 intron types in the MIDB df (minor, minor_hybrid, major_hybrid, minor-like, major-like, and non-canonical):
        #a csv with three columns:
            #sampleID: the sampleID as found in the FRASER output
            #n: the number of genes containing outliers of the specified intron type (as found in the folder name) impacted by EITHER
                #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
                #outliers of only the specified intron type (ex: number of genes with minor intron retention)
            #Z_score: the Z-score for that sample taken from the mean and sd of the n column
    #the output of make_outlier_df_gene for all 6 intron types in the MIDB df (minor, minor_hybrid, major_hybrid, minor-like, major-like, and non-canonical):
         #a folder (unique_genes_per_sample) containing one file per sample named [sampleID].csv with two columns
            #ensembl_gene_id: ENSG of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
                #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
                #outliers of only the specified intron type (ex: number of genes with minor intron retention)
            #hgnc_symbol: gene names of the gene containing the specified intron type (as found in the folder name) impacted by EITHER
                #outliers of any intron type (ex: number of minor intron containing genes, or MIGs, with intron retention of any type)
                #outliers of only the specified intron type (ex: number of genes with minor intron retention)
    #a readme containing:
            #the date this was run
            #filepath of MIDB file used
            #filepath of the FRASER output used
            #filepath of the FRASER input used

main_function <- function(filtered_FRASER_results, FRASER_input_files, minor_intron_file, dir_path, FRASER_output_fp, FRASER_input_fp, MIG_fp) {
    ###############----Define all levels----##############
    gene_minor_juncfiltered <- "genes_w_minor_intron_retention"
    gene_minor_NONfiltered <- "genes_w_minor_introns_w_intron_retention"
    junc_minor_juncfiltered <- "juncs_w_minor_intron_retention"
    junc_minor_NONfiltered <- "juncs_w_intron_retention_in_genes_w_minor_intron"

    gene_majorlike_juncfiltered <- "genes_w_major_like_intron_retention"
    gene_majorlike_NONfiltered <- "genes_w_major_like_introns_w_intron_retention"
    junc_majorlike_juncfiltered <- "juncs_w_major_like_intron_retention"
    junc_majorlike_NONfiltered <- "juncs_w_intron_retention_in_genes_w_major_like_intron"

    gene_minorlike_juncfiltered <- "genes_w_minor_like_intron_retention"
    gene_minorlike_NONfiltered <- "genes_w_minor_like_introns_w_intron_retention"
    junc_minorlike_juncfiltered <- "juncs_w_minor_like_intron_retention"
    junc_minorlike_NONfiltered <- "juncs_w_intron_retention_in_genes_w_minor_like_intron"

    gene_noncanonical_juncfiltered <- "genes_w_non_canonical_intron_retention"
    gene_noncanonical_NONfiltered <- "genes_w_non_canonical_introns_w_intron_retention"
    junc_noncanonical_juncfiltered <- "juncs_w_non_canonical_intron_retention"
    junc_noncanonical_NONfiltered <- "juncs_w_intron_retention_in_genes_w_non_canonical_genes"

    gene_majorhybrid_juncfiltered <- "genes_w_major_hybrid_intron_retention"
    gene_majorhybrid_NONfiltered <- "genes_w_major_hybrid_introns_w_intron_retention"
    junc_majorhybrid_juncfiltered <- "juncs_w_major_hybrid_intron_retention"
    junc_majorhybrid_NONfiltered <- "juncs_w_major_hybrid_in_genes_w_major_hybrid_genes"

    gene_minorhybrid_juncfiltered <- "genes_w_minor_hybrid_intron_retention"
    gene_minorhybrid_NONfiltered <- "genes_w_minor_hybrid_introns_w_intron_retention"
    junc_minorhybrid_juncfiltered <- "juncs_w_minor_hybrid_intron_retention"
    junc_minorhybrid_NONfiltered <- "juncs_w_minor_hybrid_in_genes_w_minor_hybrid_genes"

    level <- list(gene_minor_juncfiltered, gene_minor_NONfiltered, 
                        junc_minor_juncfiltered, junc_minor_NONfiltered, 
                        gene_majorlike_juncfiltered, gene_majorlike_NONfiltered,
                        junc_majorlike_juncfiltered, junc_majorlike_NONfiltered,
                        gene_minorlike_juncfiltered, gene_minorlike_NONfiltered,
                        junc_minorlike_juncfiltered, junc_minorlike_NONfiltered,
                        gene_noncanonical_juncfiltered, gene_noncanonical_NONfiltered,
                        junc_noncanonical_juncfiltered, junc_noncanonical_NONfiltered,
                        gene_majorhybrid_juncfiltered, gene_majorhybrid_NONfiltered,
                        junc_majorhybrid_juncfiltered, junc_majorhybrid_NONfiltered,
                        gene_minorhybrid_juncfiltered, gene_minorhybrid_NONfiltered,
                        junc_minorhybrid_juncfiltered,junc_minorhybrid_NONfiltered)
    
    ###############----Make output directories----##############
    level_value <- list()
    for (i in level) {
         if (grepl("genes_w", i)) {
            type <- "gene_level/"         
        }

        if (grepl("juncs_w", i)) {
            type <- "junction_level/"
        }
                
        if (grepl("minor_intron|minor_gene", i)) {
            intron <- "minor/"
        }

        if (grepl("minor_hybrid", i)) {
            intron <- "minor_hybrid/"
        }

        if (grepl("major_hybrid", i)) {
            intron <- "major_hybrid/"
        }

        if (grepl("minor_like", i)) {
            intron <- "minor_like/"
        }

        if (grepl("major_like", i)) {
            intron <- "major_like/"
        }

        if (grepl("non_canonical", i)) {
            intron <- "non_canonical/"
        }

        dir_path_level <- paste0(dir_path, intron, type, i, "/")

        if (!dir.exists(dir_path_level)) {
                print(dir_path_level)
                dir.create(dir_path_level, recursive = TRUE)
        }

        level_value[[length(level_value) + 1]] <- dir_path_level
        level_value <- unlist(level_value)
    }

    ###############----Join and filter dataframes----##############
    for (i in 1:4) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("minor", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 5:8) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("major-like", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 9:12) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("minor-like", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 13:16) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("non-canonical", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 17:20) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("major_hybrid", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)
    }

    for (i in 21:24) {
        lvl <- level_value[[i]]
        change_colnames_bind_rows("minor_hybrid", minor_intron_file, filtered_FRASER_results,
            FRASER_input = FRASER_input_files, level = lvl,index = i, dir_path)

    }

    readme_simple <- c(
        "# Run completed on:",
            paste("**", format(Sys.Date()), "**"),
        "## Filepath MIG file used:",
            paste("`", MIG_fp, "`"),
        "## Filepath of FRASER output used",
            paste("`", FRASER_output_fp, "`"),
        "## Filepath of FRASER input used",
            paste("`", FRASER_input_fp, "`")
    )
    
    writeLines(readme_simple, paste0(dir_path, "README.md"))
}

############################################################################################################
##################################----Get results FRASER 1----##############################################
############################################################################################################
UDN029422_filtered_05_abs3 <- UDN029422 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")
UDN029422_filtered_05_abs3_psi3 <- UDN029422 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi3")
UDN029422_filtered_05_abs3_psi5 <- UDN029422 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi5")

UDN112958_filtered_05_abs3 <- UDN112958 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")
UDN112958_filtered_05_abs3_psi3 <- UDN112958 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi3")
UDN112958_filtered_05_abs3_psi5 <- UDN112958 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi5")

UDN629307_filtered_05_abs3 <- UDN629307 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")
UDN629307_filtered_05_abs3_psi3 <- UDN629307 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi3")
UDN629307_filtered_05_abs3_psi5 <- UDN629307 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi5")

UDN664300_filtered_05_abs3 <- UDN664300 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")
UDN664300_filtered_05_abs3_psi3 <- UDN664300 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi3")
UDN664300_filtered_05_abs3_psi5 <- UDN664300 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi5")

UDN907359_filtered_05_abs3 <- UDN907359 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "theta")
UDN907359_filtered_05_abs3_psi3 <- UDN907359 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi3")
UDN907359_filtered_05_abs3_psi5 <- UDN907359 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3) %>% filter(type == "psi5")

########################################----UDN029422----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN029422/"
main_function(UDN029422_filtered_05_abs3, UDN029422_input, MIDB_all, dirpath, UDN029422_fp, UDN029422_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN029422/"
main_function(UDN029422_filtered_05_abs3_psi3, UDN029422_input, MIDB_all, dirpath, UDN029422_fp, UDN029422_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN029422/"
main_function(UDN029422_filtered_05_abs3_psi5, UDN029422_input, MIDB_all, dirpath, UDN029422_fp, UDN029422_input_fp, MIDB_all_fp)

########################################----UDN112958----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN112958/"
main_function(UDN112958_filtered_05_abs3, UDN112958_input, MIDB_all, dirpath, UDN112958_fp, UDN112958_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN112958/"
main_function(UDN112958_filtered_05_abs3_psi3, UDN112958_input, MIDB_all, dirpath, UDN112958_fp, UDN112958_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN112958/"
main_function(UDN112958_filtered_05_abs3_psi5, UDN112958_input, MIDB_all, dirpath, UDN112958_fp, UDN112958_input_fp, MIDB_all_fp)

########################################----UDN629307----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN629307/"
main_function(UDN629307_filtered_05_abs3, UDN629307_input, MIDB_all, dirpath, UDN629307_fp, UDN629307_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN629307/"
main_function(UDN629307_filtered_05_abs3_psi3, UDN629307_input, MIDB_all, dirpath, UDN629307_fp, UDN629307_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN629307/"
main_function(UDN629307_filtered_05_abs3_psi5, UDN629307_input, MIDB_all, dirpath, UDN629307_fp, UDN629307_input_fp, MIDB_all_fp)

########################################----UDN664300----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN664300/"
main_function(UDN664300_filtered_05_abs3, UDN664300_input, MIDB_all, dirpath, UDN664300_fp, UDN664300_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN664300/"
main_function(UDN664300_filtered_05_abs3_psi3, UDN664300_input, MIDB_all, dirpath, UDN664300_fp, UDN664300_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN664300/"
main_function(UDN664300_filtered_05_abs3_psi5, UDN664300_input, MIDB_all, dirpath, UDN664300_fp, UDN664300_input_fp, MIDB_all_fp)

########################################----UDN907359----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN907359/"
main_function(UDN907359_filtered_05_abs3, UDN907359_input, MIDB_all, dirpath, UDN907359_fp, UDN907359_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN907359/"
main_function(UDN907359_filtered_05_abs3_psi3, UDN907359_input, MIDB_all, dirpath, UDN907359_fp, UDN907359_input_fp, MIDB_all_fp)

dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN907359/"
main_function(UDN907359_filtered_05_abs3_psi5, UDN907359_input, MIDB_all, dirpath, UDN907359_fp, UDN907359_input_fp, MIDB_all_fp)

############################################################################################################
################################----Combine results FRASER 1----############################################
############################################################################################################

#############################----Combine results FRASER 1 Theta----#########################################
UDN029422_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN029422/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN029422' = 'n')

UDN112958_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN112958/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN112958' = 'n')

UDN629307_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN629307/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN629307' = 'n')

UDN664300_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN664300/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN664300' = 'n')

UDN907359_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN907359/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN907359' = 'n')

MIRE <- full_join(UDN029422_MIRE, UDN112958_MIRE) %>% full_join(UDN629307_MIRE) %>% full_join(UDN664300_MIRE) %>% full_join(UDN907359_MIRE)

##############################----Combine results FRASER 1 Psi3----#########################################
UDN029422_P3 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN029422/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN029422' = 'n')

UDN112958_P3 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN112958/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN112958' = 'n')

UDN629307_P3 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN629307/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN629307' = 'n')

UDN664300_P3 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN664300/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN664300' = 'n')

UDN907359_P3 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN907359/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN907359' = 'n')

P3 <- full_join(UDN029422_P3, UDN112958_P3) %>% full_join(UDN629307_P3) %>% full_join(UDN664300_P3) %>% full_join(UDN907359_P3)

##############################----Combine results FRASER 1 Psi3----#########################################
UDN029422_P5 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN029422/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN029422' = 'n')

UDN112958_P5 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN112958/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN112958' = 'n')

UDN629307_P5 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN629307/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN629307' = 'n')

UDN664300_P5 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN664300/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN664300' = 'n')

UDN907359_P5 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN907359/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN907359' = 'n')

P5 <- full_join(UDN029422_P5, UDN112958_P5) %>% full_join(UDN629307_P5) %>% full_join(UDN664300_P5) %>% full_join(UDN907359_P5)

##################################----Avg FRASER1 results----###############################################
MIRE <- MIRE %>% mutate(avg_MIRE = rowMeans(across(n_UDN029422:n_UDN907359), na.rm = TRUE))
P3 <- P3 %>% mutate(avg_Psi3 = rowMeans(across(n_UDN029422:n_UDN907359), na.rm = TRUE))
P5 <- P5 %>% mutate(avg_Psi5 = rowMeans(across(n_UDN029422:n_UDN907359), na.rm = TRUE))

#####################################----Get Z-scores----###################################################
MIRE_mean <- mean(MIRE$avg_MIRE)
MIRE_sd <- sd(MIRE$avg_MIRE)
MIRE$Z_score <- (MIRE$avg_MIRE - MIRE_mean)/MIRE_sd

write_csv(MIRE, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/MIRE_FRASER1.csv")

P3_mean <- mean(P3$avg_Psi3)
P3_sd <- sd(P3$avg_Psi3)
P3$Z_score <- (P3$avg_Psi3 - P3_mean)/P3_sd

write_csv(P3, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/Psi3_FRASER1.csv")

P5_mean <- mean(P5$avg_Psi5)
P5_sd <- sd(P5$avg_Psi5)
P5$Z_score <- (P5$avg_Psi5 - P5_mean)/P5_sd

write_csv(P5, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/Psi5_FRASER1.csv")

############################################################################################################
##################################----Get results FRASER 2----##############################################
############################################################################################################
UDN029422_filtered_05_abs3_2 <- UDN029422_2 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)
UDN112958_filtered_05_abs3_2 <- UDN112958_2 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)
UDN629307_filtered_05_abs3_2 <- UDN629307_2 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)
UDN664300_filtered_05_abs3_2 <- UDN664300_2 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)
UDN907359_filtered_05_abs3_2 <- UDN907359_2 %>% filter(padjust <= 0.05) %>% filter(abs(deltaPsi) >= 0.3)

########################################----UDN029422----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN029422/"
main_function(UDN029422_filtered_05_abs3_2, UDN029422_input, MIDB_all, dirpath, UDN029422_fp2, UDN029422_input_fp, MIDB_all_fp)

########################################----UDN112958----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN112958/"
main_function(UDN112958_filtered_05_abs3_2, UDN112958_input, MIDB_all, dirpath, UDN112958_fp2, UDN112958_input_fp, MIDB_all_fp)

########################################----UDN629307----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN629307/"
main_function(UDN629307_filtered_05_abs3_2, UDN629307_input, MIDB_all, dirpath, UDN629307_fp2, UDN629307_input_fp, MIDB_all_fp)

########################################----UDN664300----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN664300/"
main_function(UDN664300_filtered_05_abs3_2, UDN664300_input, MIDB_all, dirpath, UDN664300_fp2, UDN664300_input_fp, MIDB_all_fp)

########################################----UDN907359----##################################################
dirpath <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN907359/"
main_function(UDN907359_filtered_05_abs3_2, UDN907359_input, MIDB_all, dirpath, UDN907359_fp2, UDN907359_input_fp, MIDB_all_fp)

############################################################################################################
################################----Combine results FRASER 2----############################################
############################################################################################################
UDN029422_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN029422/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN029422' = 'n')

UDN112958_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN112958/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN112958' = 'n')

UDN629307_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN629307/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN629307' = 'n')

UDN664300_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN664300/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN664300' = 'n')

UDN907359_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN907359/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_UDN907359' = 'n')

Jaccard_FRASER2 <- full_join(UDN029422_J, UDN112958_J) %>% full_join(UDN629307_J) %>% full_join(UDN664300_J) %>% full_join(UDN907359_J)

##################################----Avg FRASER1 results----###############################################
Jaccard_FRASER2 <- Jaccard_FRASER2 %>% mutate(avg_Jaccard = rowMeans(across(n_UDN029422:n_UDN907359), na.rm = TRUE))

Jaccard_FRASER2 %>% arrange(desc(avg_Jaccard))

Jaccard_FRASER2 %>% filter(sampleID == "UDN907359")

#####################################----Get Z-scores----###################################################
J_mean <- mean(Jaccard_FRASER2$avg_Jaccard)
J_sd <- sd(Jaccard_FRASER2$avg_Jaccard)
Jaccard_FRASER2$Z_score <- (Jaccard_FRASER2$avg_Jaccard - J_mean)/J_sd

write_csv(Jaccard_FRASER2, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/MI_Jaccard_FRASER2.csv")

############################################################################################################
################################----Join FRASER1 and FRASER2----############################################
############################################################################################################
J_results <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/MI_Jaccard_FRASER2.csv") %>%
                select(sampleID, avg_Jaccard, Z_score) %>%
                rename("Z_score_Jaccard" = "Z_score")
MIRE_results <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/MIRE_FRASER1.csv") %>%
                select(sampleID, avg_MIRE, Z_score) %>%
                rename("Z_score_MIRE" = "Z_score")
P3_results <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/Psi3_FRASER1.csv") %>% 
                select(sampleID, avg_Psi3, Z_score) %>%
                rename("Z_score_Psi3" = "Z_score")
P5_results <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/Psi5_FRASER1.csv") %>%
                select(sampleID, avg_Psi5, Z_score)  %>%
                rename("Z_score_Psi5" = "Z_score")

joined_results <- left_join(J_results, MIRE_results) %>% left_join(P3_results) %>% left_join(P5_results)

joined_results <- joined_results %>% mutate(type = case_when(sampleID %in% c("UDN029422", "UDN112958", "UDN907359") ~ "RNU12-opathy",
                                                        sampleID %in% c("UDN629307", "UDN664300") ~ "RNU12_unaffected_family",
                                                        !sampleID %in% c("UDN029422", "UDN112958", "UDN907359", "UDN629307", "UDN664300") ~ NA))

write_csv(joined_results, "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/joined_results.csv")

joined_results %>% filter(type == "RNU12-opathy")

############################################################################################################
########################################----Analyze Results----#############################################
############################################################################################################
results <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/joined_results.csv")

results %>% arrange(desc(Z_score_MIRE))

