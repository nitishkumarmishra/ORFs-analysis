# ##############################################################################
# Author: Nitish K. Mishra
# Copyright (c) Nitish K. Mishra, 2023
# Email:  nitishimtech@gmail.com
# Date: 15 February, 2023
# Script Name: run_homer_analysis.R
################################################################################
# ########################  Script Description: ################################
# This code will perform analysis by using R code "function_homer_by_nitish.R"
################################## Notes #######################################
# This is the R code for running RBP binding motif analysis by unig RBP motifs.
# We used motifs from 
################################################################################
########################## SET WORKING DIRECTORY ###############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
###### change "wd" accordingly based on project current working directory ######
wd <- "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/data/mgi"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
################################################################################

################################################################################
suppressMessages(suppressWarnings(source("function_homer_by_nitish.R")))
# ******************************************************************************
# ******************** make *****************************

species <- "m39"
myData.species  <-  species
myData.species
myData.transcriptome.gtf <- load_transcriptome_gtf(species = myData.species)
myData.transcriptome.gtf <- prepare_transcriptome_gtf(gtf = myData.transcriptome.gtf)
myData.transcriptome.info <- get_all_transciptome_classification_info(gtf = myData.transcriptome.gtf)

gene_set <- build_gene_set_with_rnaseq_riboseq()
verbose <- TRUE


################################################################################
################################################################################
# Make motif dataframe from Bioconductor package "trnasite": https://transite.mit.edu/ 
# https://www.sciencedirect.com/science/article/pii/S2211124720310494?via%3Dihub
# Add three extra motif manually which used by HK i.e. "TISU","GGAGG", "GGC-rich"

suppressMessages(suppressWarnings(library(transite)))
# ******************************************************************************

motif_all <- get_motifs_meta_info()
motiffName <- paste0(motif_all$id,"_", gsub(",.*$", "", motif_all$rbps))
motifSeq <- motif_all$iupac
# transite motif database is for the RNA while HOMER create PWM matrix for DNA motif. 
# So replace "U" with "T" to make RNA motif sequences in DNA sequences.
motifSeq <- gsub("U", "T", motifSeq)
motifSeq <- motif_list <- as.data.frame(cbind(motiffName, motifSeq))


EXTRA <- data.frame(motiffName=c("TISU","GGAGG", "GGC-rich"),
           motifSeq=c("SAAGATGGCGGC", "GGAGG","GCGGCGGCGGCG"))
motif_list <- rbind(EXTRA,motif_list)


################################################################################
################################################################################
# Motif number 136 and 137 have some issue. Below code stop running after 135.
# Below I just rerun manually from 138 
# ******************************************************************************
# ******************************************************************************

# for (idx in seq(nrow(motif_list))){
#   
#   #set variables to loop through
#   name <- motif_list[idx,1]
#   print(name)
#   motif <- motif_list[idx,2]
#   print(motif)
#   
#   list_tests <- perform_tests(c(motif_name = name, motif = motif), 
#                               c("five_prime_utr", "CDS", "5UTR-CDS", "three_prime_utr"), 
#                               "rna-binding_protein_mrna_motif", 
#                               gene_set = gene_set, 
#                               dir_out = paste0("2023.emt.rprof/", name), 
#                               fname_prefix = paste0("emt.rprof.", name))
#   list_tests$tdf
#   
#   
#   list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
#                           row_pattern="^(?!.*translationONLY.*)", 
#                           col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
#                           cluster_rows = T, column_title=name, 
#                           fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
#   fname  <-  sprintf("%s_%s", name, list_out$type_col_short )
#   print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)
# }

# Either I have to remove 136-137 motif from dataframe or run 1:135 & 138:nrow(motif_list)  
for (idx in 1:135){
  
  #set variables to loop through
  name <- motif_list[idx,1]
  print(name)
  motif <- motif_list[idx,2]
  print(motif)
  
  list_tests <- perform_tests(c(motif_name = name, motif = motif), 
                              c("five_prime_utr", "CDS", "5UTR-CDS", "three_prime_utr"), 
                              "rna-binding_protein_mrna_motif", 
                              gene_set = gene_set, 
                              dir_out = paste0("2023.emt.rprof/", name), 
                              fname_prefix = paste0("emt.rprof.", name))
  list_tests$tdf
  
  
  list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
                          row_pattern="^(?!.*translationONLY.*)", 
                          col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
                          cluster_rows = T, column_title=name, 
                          fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
  fname  <-  sprintf("%s_%s", name, list_out$type_col_short )
  print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)
}

#for (idx in 138:140){
for (idx in 138:nrow(motif_list)){
   
  #set variables to loop through
  name <- motif_list[idx,1]
  print(name)
  motif <- motif_list[idx,2]
  print(motif)
  
  list_tests <- perform_tests(c(motif_name = name, motif = motif), 
                              c("five_prime_utr", "CDS", "5UTR-CDS", "three_prime_utr"), 
                              "rna-binding_protein_mrna_motif", 
                              gene_set = gene_set, 
                              dir_out = paste0("2023.emt.rprof/", name), 
                              fname_prefix = paste0("emt.rprof.", name))
  list_tests$tdf
  
  
  list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
                          row_pattern="^(?!.*translationONLY.*)", 
                          col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
                          cluster_rows = T, column_title=name, 
                          fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
  fname  <-  sprintf("%s_%s", name, list_out$type_col_short )
  print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)
}

################################################################################
################################################################################
suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(suppressWarnings(library("data.table")))
suppressMessages(suppressWarnings(library("dplyr")))
# ******************************************************************************
# ******************************************************************************
# Read all motif analysis result files.
# Select rows which have abs(dlog10p) >= 1.30103), and make subset by using select

data_dir <- "2023.emt.rprof"

csv_files <- fs::dir_ls(data_dir, regexp = "\\.tests.plot.tdf.txt", recurse = TRUE)


my_data <- csv_files %>% 
  purrr::map_dfr(fread, .id = "source") %>%
  mutate(source=sapply(strsplit(source, "/"),  "[", 2)) %>%
  filter(abs(dlog10p) >= 1.30103) %>%
  select(c(source, class, group, direc, dlog10p, geneset)) %>%
  dplyr::rename("RBP_motif_name"="source") %>%
  arrange(RBP_motif_name, desc(dlog10p)) %>%
  mutate(RBF_class=paste0(geneset, "_", group)) #%>%
  #write.xlsx(my_data, "RBP_motif_analysis.xlsx")
# If I will use xlsx::write.xlsx in %>% then it will not make my_data. 
# So, I will make my_data then save it by using xlsx::write.xlsx
xlsx::write.xlsx(my_data, "RBP_motif_analysis.xlsx")



################################################################################
################################################################################
# Further filtering based on RBP motif which have inverse relation
# Select which have opposite dlog10p in Up/Down or UpDown/DownUp
# RBP motifs in different RBF classes (geneset)
# Adapted https://stackoverflow.com/questions/57171967
# filter(n_distinct(sign(dlog10p)) == 2) check both in up/down it's not over/underrepresented
# In unt48.tgfb48.DEtranslation_up/unt48.tgfb48.DEtranslation_down both should not be under/over-represented

# ******************************************************************************
# ******************************************************************************

# -------- "three_prime_utr" -----------
RBP.motifs.three_prime_utr.DEtranslation <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="three_prime_utr" & geneset=="unt48.tgfb48.DEtranslation") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% # in both unt48.tgfb48.DEtranslation_up & unt48.tgfb48.DEtranslation_down
  #head() # Total 35 RBPs are available in 3'UTR
  

RBP.motifs.three_prime_utr.translation.CX <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="three_prime_utr" & geneset=="reversible.translation.CX") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% 
  #head(n=30) # Total 35 RBPs are available in 3'UTR


# -------- "five_prime_utr" -----------
RBP.motifs.five_prime_utr.DEtranslation <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="five_prime_utr" & geneset=="unt48.tgfb48.DEtranslation") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% 
  #head(n=10) ## Total 5 RBPs are available in 5'UTR


RBP.motifs.five_prime_utr.translation.CX <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="five_prime_utr" & geneset=="reversible.translation.CX") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% 
  #head(n=30) # Total 35 RBPs are available in 3'UTR

# -------------- "CDS" ----------------
RBP.motifs.CDS.DEtranslation <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="CDS" & geneset=="unt48.tgfb48.DEtranslation") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% # in both unt48.tgfb48.DEtranslation_up & unt48.tgfb48.DEtranslation_down
#head() # Total 35 RBPs are available in 3'UTR


RBP.motifs.CDS.translation.CX <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="CDS" & geneset=="reversible.translation.CX") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% 
#head(n=30) # Total 35 RBPs are available in 3'UTR


# ------------ "5UTR-CDS" -------------
RBP.motifs.5UTR_CDS.DEtranslation <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="5UTR-CDS" & geneset=="unt48.tgfb48.DEtranslation") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% # in both unt48.tgfb48.DEtranslation_up & unt48.tgfb48.DEtranslation_down
#head() # Total 35 RBPs are available in 3'UTR


RBP.motifs.5UTR_CDS.translation.CX <- 
  my_data %>% 
  group_by(RBP_motif_name) %>%
  filter(class=="5UTR-CDS" & geneset=="reversible.translation.CX") %>%
  filter(n_distinct(sign(dlog10p)) == 2) #%>% 
#head(n=30) # Total 35 RBPs are available in 3'UTR


################################################################################
################################################################################
# Don't need this part as I already added these three motif in motif_all
## TISU notif
list_tests <- perform_tests(c(motif_name = "TISU", motif = "SAAGATGGCGGC"), c("five_prime_utr", "CDS", "5UTR-CDS", "three_prime_utr"), "rna-binding_protein_mrna_motif", 
                            gene_set = gene_set, dir_out = "2023.emt.rprof/TISU", fname_prefix = "emt.rprof.TISU")
list_tests$tdf

list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
                        row_pattern="^(?!.*translationONLY.*)", 
                        col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
                        cluster_rows = T, column_title='TISU', 
                        fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
fname  <-  sprintf("TISU_%s", list_out$type_col_short )
print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)


## CGAGG motif
list_tests <- perform_tests(c(motif_name = "GGAGG", motif = "GGAGG"), c("five_prime_utr", "CDS", "three_prime_utr"), "rna-binding_protein_mrna_motif", 
                            gene_set = gene_set, dir_out = "2023.emt.rprof/GGAGG", fname_prefix = "emt.rprof.GGAGG")
list_tests$tdf

list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
                        row_pattern="^(?!.*translationONLY.*)", 
                        col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
                        cluster_rows = T, column_title='GGAGG', 
                        fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
fname  <-  sprintf("GGAGG_%s", list_out$type_col_short )
print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)


## eIF4A1 long 5'UTR GC Rich motif
list_tests <- perform_tests(c(motif_name = "GGC-rich", motif = "GCGGCGGCGGCG"), c("five_prime_utr", 
                                                                                  "CDS", "three_prime_utr"), "rna-binding_protein_mrna_motif", gene_set = gene_set, 
                            dir_out = "2023.emt.rprof/GGC-rich", fname_prefix = "emt.rprof.GGC-rich")
list_tests$tdf


list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
                        row_pattern="^(?!.*translationONLY.*)", 
                        col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
                        cluster_rows = T, column_title='GGC-rich', 
                        fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
fname  <-  sprintf("GGC-rich_%s", list_out$type_col_short )
print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)


## 5'UTR TOP motif
# https://www.mdpi.com/2218-273X/10/7/969
# https://www.pnas.org/doi/10.1073/pnas.1912864117
list_tests <- perform_tests(c(motif_name = "TOP-motif", motif = "CTCTTCC"), c("five_prime_utr", 
                                                                                  "CDS", "three_prime_utr"), "rna-binding_protein_mrna_motif", gene_set = gene_set, 
                            dir_out = "2023.emt.rprof/TOP-motif", fname_prefix = "emt.rprof.TOP-motif")
list_tests$tdf


list_out <- get_heatmap(list_tests$tmat, vec_color_group = NULL, type_col = "dlog10p", 
                        row_pattern="^(?!.*translationONLY.*)", 
                        col = colorRamp2(c(-5,-1,0,1,5),c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000")),
                        cluster_rows = T, column_title='TOP-motif', 
                        fontsize_row_names=9, fontsize_column_names=9, column_name_angle=50)
fname  <-  sprintf("TOP-motif_%s", list_out$type_col_short )
print_figure(list_out$list_ht, width = 4.25, height = 2.5, file = fname)
