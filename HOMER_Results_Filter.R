# ##############################################################################
# Author: Nitish K. Mishra
# Copyright (c) Nitish K. Mishra, 2023
# Email:  nitishimtech@gmail.com
# Date: 14 March, 2023
# Script Name: HOMER_Results_Filter.R
################################################################################
# ########################  Script Description: ################################
# R code for selecting RBPs HOMER motif from run_homer_analysis.R analysis result .txt files. 
################################## Notes #######################################
# This code will make a dataframe from HOMER results .txt file (from run_homer_analysis.R)
# Further we filter dataframe based on logP (-dlog10p)
################################################################################
########################## SET WORKING DIRECTORY ###############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
###### change "wd" accordingly based on project current working directory ######
wd <- "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/data/mgi"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
################################################################################
# load global R library
suppressMessages(suppressWarnings(library("tidyverse")))

################################################################################
# Make list of HOMER motif result text file and make dataframe for the analysis
data_dir <- "2023.emt.rprof"
csv_files <- fs::dir_ls(data_dir, regexp = "\\.tests.plot.tdf.txt", recurse = TRUE)

################################################################################
# Filter result- select RBP motif which have HOMER p-vale <= 0.001 (dlog10p >=3)
# Although P-value cutoff is very low, but it will avaoid random motif

my_data = csv_files |> 
  purrr::map_dfr(fread, .id = "source") |>
  mutate(source=sapply(strsplit(source, "/"),  "[", 2)) |>
  #filter(abs(dlog10p) >= 1.30103) |> # P-value <= 0.05
  filter(abs(dlog10p) >= 3) |> # P-value <= 0.001
  select(c(source, class, direc, dlog10p, geneset, group)) |>
  dplyr::rename("RBP_motif_name"="source") |>
  arrange(RBP_motif_name, desc(dlog10p)) |>
  mutate(Extra=paste0(RBP_motif_name, "-", class, "-", geneset)) |>
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  mgsub::mgsub(pattern = c("GGAGG", "TOP-motif", "GCG-rich", "TISU"), replacement = c("GGAGG_motif_GGAGG", "TOP-motif_1_TOP-motif", "GCG-rich_1_GCG-rich","TISU_1_TISU"))
  

# Total RBPs in all classes 
#class == "three_prime_utr" | "five_prime_utr" | "5UTR-CDS"|"CDS"
my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS"|class == "CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Total 182 RBP motifs are enriched
length(unique(my_data_new$RBP_motif_name)) # Total 93 unique RBPs motifs are enriched in all categories

a3 <- unique(sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3))
RBP_motif_all <- a3
length(RBP_motif_all) # Total 63 unique RBPs in 3'UTR (including CDS)
RBP_motif_all

# Filter RBPs based on class 
#class == "three_prime_utr"|"five_prime_utr"|"5UTR-CDS"

my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # total 56 
length(unique(my_data_new$RBP_motif_name)) # Total 49 RBPs HOMER motifs are enriched all 3'UTR, 5'UTR and 5'UTR-CDS

a3 <- unique(sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3))
RBP_motif_all_no_CDS <- a3
length(RBP_motif_all_no_CDS) # Total 38 unique RBPs in 3'UTR (excluding CDS)
RBP_motif_all_no_CDS

######################### ************************* ############################
                        ### 3' UTR motif results ###
######################### ************************* ############################
# 3'UTR RBP motifs in both geneset == "unt48.tgfb48.DEtranslation" | "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  #filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "three_prime_utr") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Total 37 RBPs motifs are enriched
length(unique(my_data_new$RBP_motif_name)) # Total 36 unique RBPs motifs are enriched
# HNRNPR is enriched in both "unt48.tgfb48.DEtranslation" and "reversible.translation.CX"

a3 <- unique(sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3))
utr3_RBP_motif_all <- a3
length(utr3_RBP_motif_all) # Total 27 unique RBPs in 3'UTR geneset == "unt48.tgfb48.DEtranslation" | "reversible.translation.CX"
utr3_RBP_motif_all

######################### ************************* ############################
######################### ************************* ############################
# RBP motif in 3'UTR in geneset == "unt48.tgfb48.DEtranslation"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "three_prime_utr") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
length(unique(my_data_new$RBP_motif_name)) # Total five motif

# Name of all (five) RBP motifs for class == "three_prime_utr" && geneset == "unt48.tgfb48.DEtranslation"
a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
utr3_RBP_motif_DEtranslation <- paste0(a1, "_", a3)

utr3_RBPs_all_DEtranslation <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) # 5 RBPs
utr3_RBPs_all_DEtranslation
table(utr3_RBPs_all_DEtranslation) # SRSF12 has two motifs

utr3_RBPs_uniq_DEtranslation <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # 4 unique RBPs
utr3_RBPs_uniq_DEtranslation


######################### ************************* ############################
######################### ************************* ############################
# RBP motif in 3'UTR in geneset == "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "reversible.translation.CX") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "three_prime_utr") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
length(unique(my_data_new$RBP_motif_name)) # Total five motif

# Name of all (five) RBP motifs for class == "three_prime_utr" && geneset == "unt48.tgfb48.DEtranslation"
a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
utr3_RBP_motif_reversible <- paste0(a1, "_", a3)

utr3_RBPs_all_reversible <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) 
length(utr3_RBPs_all_reversible) # 32 RBP motifs are enriched
utr3_RBPs_all_reversible
table(utr3_RBPs_all_reversible) # SRSF1=5 and  EIF4B, SF3B4, and TIA1 has two motifs
# SRSF1 = 5, EIF4B = 2, SF3B4 = 2, and TIA1 = 2

utr3_RBPs_uniq_reversible <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # 25 unique RBPs have motif enrichment
length(utr3_RBPs_uniq_reversible) # Total 4 RBPs, all are unique
utr3_RBPs_uniq_reversible


######################### ************************* ############################
                        ### 5' UTR motif results ###
######################### ************************* ############################
# 5'UTR RBP motifs in both geneset == "unt48.tgfb48.DEtranslation" | "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  #filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "five_prime_utr") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)

table(my_data_new$RBP_motif_name) # Total 17 
length(unique(my_data_new$RBP_motif_name)) # Total 16 RBPs HOMER motifs are enriched
# PABPN1 is enriched in both "unt48.tgfb48.DEtranslation" and "reversible.translation.CX"

######################### ************************* ############################
######################### ************************* ############################
# RBP motif in 5'UTR in geneset == "unt48.tgfb48.DEtranslation"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "five_prime_utr") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Total 16 RBP motifs are enriched
table(my_data_new$RBP_motif_name)
length(unique(my_data_new$RBP_motif_name)) # Total 16 RBP motifs are enriched

# Name of all (five) RBP motifs for class == "three_prime_utr" && geneset == "unt48.tgfb48.DEtranslation"
a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
utr5_RBP_motif_DEtranslation <- paste0(a1, "_", a3) 
length(utr5_RBP_motif_DEtranslation) # Total 16 RBP motifs

utr5_RBPs_all_DEtranslation <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) # 16 RBPs
length(utr5_RBPs_all_DEtranslation) # Total 16 RBPs
utr5_RBPs_all_DEtranslation
table(utr5_RBPs_all_DEtranslation) # RBM4B has two motifs

utr5_RBPs_uniq_DEtranslation <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # 15 unique RBPs
length(utr5_RBPs_uniq_DEtranslation) # Total 15 unique RBPs
utr5_RBPs_uniq_DEtranslation


######################### ************************* ############################
######################### ************************* ############################
# RBP motif in 5'UTR in geneset == "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "reversible.translation.CX") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "five_prime_utr") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Only RBP PABPN1 motif
length(unique(my_data_new$RBP_motif_name)) # Only RBP PABPN1 motif

# Name of all (five) RBP motifs for class == "three_prime_utr" && geneset == "unt48.tgfb48.DEtranslation"
a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
utr5_RBP_motif_reversible <- paste0(a1, "_", a3) # Only M148_PABPN1

utr5_RBPs_all_reversible <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) 
length(utr5_RBPs_all_reversible) # Only RBP PABPN1 motifs are enriched
utr5_RBPs_all_reversible
table(utr5_RBPs_all_reversible) # Only RBP PABPN1

utr3_RBPs_uniq_reversible <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # only RBPs PABPN1 has motif enrichment
length(utr3_RBPs_uniq_reversible)
utr3_RBPs_uniq_reversible # only PABPN1



######################### ************************* ############################
                        ## 5UTR-CDS motif results ##
######################### ************************* ############################
# 5'UTR-CDS RBP motifs in both geneset == "unt48.tgfb48.DEtranslation" | "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  #filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "5UTR-CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)

table(my_data_new$RBP_motif_name) # Total 2 RBP motif
length(unique(my_data_new$RBP_motif_name)) # Total 2 RBPs motifs are enriched

######################### ************************* ############################
######################### ************************* ############################
# RBP motif in 5'UTR-CDS in geneset == "unt48.tgfb48.DEtranslation"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "5UTR-CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Total 2 RBP motifs are enriched
table(my_data_new$RBP_motif_name)
length(unique(my_data_new$RBP_motif_name)) # Total 2 RBP motifs are enriched


a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
utr5_CDS_RBP_motif_DEtranslation <- paste0(a1, "_", a3) # Total 2 RBP motifs

utr5_CDS_RBPs_all_DEtranslation <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) # 2 RBPs
length(utr5_CDS_RBPs_all_DEtranslation) # Tota 2 RBP
utr5_CDS_RBPs_all_DEtranslation
table(utr5_CDS_RBPs_all_DEtranslation) 

utr5_CDS_RBPs_uniq_DEtranslation <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # 2 unique RBPs
length(utr5_CDS_RBPs_uniq_DEtranslation) # Both proteins are unique
utr5_CDS_RBPs_uniq_DEtranslation


######################### ************************* ############################
######################### ************************* ############################
# RBP motif in 5'UTR-CDS in geneset == "reversible.translation.CX"
# No RBP motif enriched in "reversible.translation.CX" class == "5UTR-CDS"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "reversible.translation.CX") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "5UTR-CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Only RBP PABPN1 motif
length(unique(my_data_new$RBP_motif_name)) # Only RBP PABPN1 motif

# No RBP motif enriched in "reversible.translation.CX" class == "5UTR-CDS"



######################### ************************* ############################
                          ### CDS motif results ###
######################### ************************* ############################
# CDS RBP motifs in both geneset == "unt48.tgfb48.DEtranslation" | "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  #filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Total 126 RBP motifs
table(my_data_new$RBP_motif_name) # Total 79 RBS have 126 notifs 
length(unique(my_data_new$RBP_motif_name)) # Total 126 RBPs motifs are enriched
# Total 47 RBP motifs are enriched in both "unt48.tgfb48.DEtranslation" and "reversible.translation.CX"
# Total 32 enriched in only one class

######################### ************************* ############################
######################### ************************* ############################
# RBP motif in CDS in geneset == "unt48.tgfb48.DEtranslation"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "unt48.tgfb48.DEtranslation") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Total 63 RBP motifs are enriched
table(my_data_new$RBP_motif_name)
length(unique(my_data_new$RBP_motif_name)) # Total 63 RBP motifs are enriched

# Name of all (five) RBP motifs for class == "three_prime_utr" && geneset == "unt48.tgfb48.DEtranslation"
a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
CDS_RBP_motif_DEtranslation <- paste0(a1, "_", a3) # Total 63 RBP motifs

CDS_RBPs_all_DEtranslation <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) # 63 RBPs
length(CDS_RBPs_all_DEtranslation)
CDS_RBPs_all_DEtranslation
table(CDS_RBPs_all_DEtranslation) # [SRSF1, SRSF12, RBM38] have 4; [PCBP4 and SF3B4] have 3
# A2BP1, CSDA, EIF4B, HNRNPL, RBFOX2, RBM47, and SRSF9 have two RBP motifs

CDS_RBPs_uniq_DEtranslation <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # 43 unique RBPs
length(CDS_RBPs_uniq_DEtranslation) # 43 unique RBPs
CDS_RBPs_uniq_DEtranslation


######################### ************************* ############################
######################### ************************* ############################
# RBP motif in CDS in geneset == "reversible.translation.CX"

my_data_new <- 
  my_data |> 
  #filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(geneset == "reversible.translation.CX") |>
  #filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  filter(class == "CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)
dim(my_data_new) # Only RBP PABPN1 motif
length(unique(my_data_new$RBP_motif_name)) # Only RBP PABPN1 motif

# Name of all (five) RBP motifs for class == "three_prime_utr" && geneset == "unt48.tgfb48.DEtranslation"
a3 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3)
a1 <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',1)
utr5_CDS_motif_reversible <- paste0(a1, "_", a3)

CDS_RBPs_all_reversible <- sapply(strsplit(my_data_new$RBP_motif_name, "_"), '[[',3) 
length(CDS_RBPs_all_reversible) # Total 63 RBP motifs are enriched
CDS_RBPs_all_reversible
sort(table(CDS_RBPs_all_reversible)) # SRSF1:6; CSDA and RBM38 : 4
# EIF4B    HNRNPL    LIN28A     PCBP4     RBM47      RBM5     SNRPA     SRSF9 :: 2

CDS_RBPs_uniq_reversible <- sort(unique(sapply(strsplit(unique(my_data_new$RBP_motif_name), "_"), '[[',3))) # Total 42 RBPs have 62 motif enrichment
length(CDS_RBPs_uniq_reversible)
CDS_RBPs_uniq_reversible # Total 42 RBPs



################################################################################
# Find splicing proteins and splicing repressor HNRN proteins

#my_data_new[grep("SRSF", my_data_new$RBP_motif_name),] # 10 RBP motifs (4 SRSFs proteins)
#my_data_new[grep("HNRN", my_data_new$RBP_motif_name),] # 5 RBP motifs (3 HNRNPs proteins)

my_data_new <- 
  my_data |> 
  filter(geneset == "unt48.tgfb48.DEtranslation"|geneset == "reversible.translation.CX") |>
  filter(class == "three_prime_utr" | class=="five_prime_utr" |class == "5UTR-CDS") |>
  group_by(Extra) |> 
  mutate(MotifNumber= n()) |> 
  distinct(MotifNumber, .keep_all=TRUE) |>
  filter(MotifNumber >= 2) |> # Select RBPs with p-value <= 0.01 in both up/down geneset
  select(!direc)

my_data_new |>
  filter(grepl("SRSF|HNRN", RBP_motif_name))
# Unique 16 HNRNPR are enriched in both "unt48.tgfb48.DEtranslation" and "reversible.translation.CX"

################################################################################
################################################################################

save.image("HOMER_Results_Filter.RData")

################################################################################
################################################################################