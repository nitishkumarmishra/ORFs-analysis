# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2022
# Email:  nitishimtech@gmail.com
# Date: 05 December, 2022
# Script Name: make_GRCm39_104_transcript_gtf.R
#####################################################################
# ##################  Script Description: ###########################
# Script is based on Matt Park's "/home/map2085/riboprof/src/make_eukaryotic_transcript_gtf.R"
##################### Notes: ########################################
#
#####################################################################
################## SET WORKING DIRECTORY ############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/jupyter codes"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################################
# source("make_GRCm39_104_transcript_gtf.R")

## Modified by Nitish on Dec 05, 2022
suppressMessages(suppressWarnings(source("./differential_expression_analysis.R")))
suppressMessages(suppressWarnings(source("./stjude_riboprof_common_tests.R")))
suppressMessages(suppressWarnings(source("./functions_homer.R")))

# Original Matt/Hyunsoo code is for Linux system.
# The function read_gtf in differential_expression_analysis.R is for linux zcat based, so it will not work on windows
# So I make read_gtf for windows
## Modified by Nitish on Dec 05, 2022

# GTF file "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/jupyter codes/data/mouse/Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz")
# outfname <- "Mus_musculus.GRCm39.104"
# source("make_GRCm39_104_transcript_gtf.R")


##########  load R packages
library(stringr)
library(dplyr)
library(ggplot2)

################################################################################
########## ************* Modified by Nitish on Dec 05, 2022 ********** #########

gtf  <-  "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/jupyter codes/data/mouse/Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz"
outfname  <-  "/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/jupyter codes/data/mouse/Mus_musculus.GRCm39.104.rdna_rn18s.txt"

# read gtf file
read_gtf  <-  function( file , header = FALSE ) {
  message('[read_gtf][file] ', file)
  in.gtf  <-  fread( file = file , sep = "\t" , header = FALSE ,  stringsAsFactors = FALSE , quote = "" , colClasses = c( rep("character",3) , rep("numeric",2) , rep("character",4) ) )
  
  def.colnames  <-  c("seqname" , "source" , "feature" , "start" , "end" , "score" , "strand" , "frame" , "attribute")
  colnames(in.gtf)  <-  def.colnames
  
  in.gtf  <-  as.data.frame(in.gtf)
  
  return(in.gtf)
  
} # read_gtf

## Modified by Nitish on Dec 05, 2022

make_eukaryotic_transcript_gtf  <-  function(
    gtf = NULL ,
    outfname = NULL )  {
  
  ### DESCRIPTION:
  #	Given an Ensembl gtf for a eukaryotic organism, this transcript will create a transcript gtf.

  ######################################## DEBUG
  ############## NOTE
  # 5' & 3' UTRs are nested inside of the first & last exons.
  
  ######################## load
  verb("\tload.\n")
  
  myData.gtf  <-  read_gtf( file = gtf )
  
  myData.colnames  <-  colnames(myData.gtf)
  
  myData.tranGTF  <-  data.frame()
  
  # DEBUG
  #myData.SAVE  <-  myData.gtf  
  
  
  ### feature per feature
  verb("\tfeature per feature.\n")
  
  myData.gtf$transcript_id  <-  get_gtf_attribute_field( gtf = myData.gtf , field = "transcript_id" )
  myData.gtf$exon_number  <-  get_gtf_attribute_field( gtf = myData.gtf , field = "exon_number" , num = TRUE )
  
  
  ### feat width
  myData.gtf$width  <-  myData.gtf$end - myData.gtf$start + 1
  
  
  ###### transcript length
  verb("\ttranscript length.\n")
  
  # Note that "exon" is anything other than intron, i.e. anything transcribed, even pseudogenes.  So the exons = 5' UTR & CDS & 3' UTR.
  # Thus, summing up exon lengths will give you the full transcript length.
  sub.feat  <-  myData.gtf[ myData.gtf$feature %in% c("exon") , ,drop=FALSE]
  sub.feat  <-  sub.feat  %>%  dplyr::group_by(transcript_id)  %>%  dplyr::mutate(transcript.length = sum(width) )
  sub.feat  <-  as.data.frame(sub.feat)
  sub.feat  <-  unique(sub.feat[, c("transcript_id","transcript.length")])
  
  myData.gtf  <-  merge( myData.gtf , sub.feat , by = "transcript_id" , all.x = TRUE)
  
  
  
  #### basic filter
  verb("\tbasic filter.\n")
  
  #myData.gtf  <-  myData.gtf[ !is.na(myData.gtf$transcript_id)  &  !is.na(myData.gtf$transcript.length)  &  myData.gtf$transcript.length > 0 ,]
  myData.gtf  <-  myData.gtf[ !is.na(myData.gtf$transcript_id) , ,drop=FALSE]
  
  if ( any(is.na(myData.gtf$transcript.length))  ||  any(myData.gtf$transcript.length == 0) ) {
    verb("\n\n\nERROR!  NA or 0 transcript lengths!!!!\n")
    
    verb("NA lengths:\n")
    show(head(myData.gtf[ is.na(myData.gtf$transcript.length) , ,drop=FALSE ] ))
    
    verb("\n\n\n0 lengths:\n")
    show(head(myData.gtf[ myData.gtf$transcript.length == 0 , ,drop=FALSE ] ))
    
    stop()
  } # check
  
  
  
  #### UTR cs
  verb("\tUTR widths.\n")
  
  ### UTRs may be split into several pieces with introns inbetween,
  # spliced together during splicing.
  
  for (utrx  in  c("five_prime_utr","three_prime_utr")) {
    verb("\t\t%s\n", utrx)
    
    widfield  <-  paste( utrx , "width" , sep=".")
    
    sub.utr  <-  myData.gtf[ myData.gtf$feature == utrx , ,drop=FALSE]
    
    if (nrow(sub.utr) == 0) {
      myData.gtf[[widfield]]  <-  0
    } else {
      sub.utr  <-  sub.utr[,c("transcript_id","width")]
      sub.utr  <-  sub.utr  %>%  dplyr::group_by(transcript_id)  %>%  dplyr::summarize(utr.width = sum(width))
      sub.utr  <-  as.data.frame(sub.utr)
      sub.utr[[widfield]]  <-  sub.utr$utr.width
      sub.utr$utr.width  <-  NULL
      
      myData.gtf  <-  merge( myData.gtf , sub.utr , by = "transcript_id" , all.x = TRUE)
      myData.gtf[[widfield]][ is.na(myData.gtf[[widfield]]) ]  <-  0
    } # has this UTR
  } # utrx
  
  
################################################################################
############## ********* cumulative relative position ********* ################
  
verb("\tcumulative relative position.\n")
  
sub.exon  <-  myData.gtf[ myData.gtf$feature == "exon" , ]
sub.exon  <-  sub.exon[ order(sub.exon$transcript_id , sub.exon$exon_number) , ,drop=FALSE]
  
### orient start
verb("\t\torient start.\n")
  
sub.exon$exon.orient.start  <-  ifelse( sub.exon$strand == "+" , sub.exon$start , sub.exon$end )
  
  
### cumsum
verb("\t\tcumsum.\n")
  
sub.exon  <-  sub.exon  %>%  dplyr::group_by(transcript_id)  %>%  dplyr::mutate(cs.pos = cumsum(width) - width)
sub.exon  <-  as.data.frame(sub.exon)
  
#logi.na  <-  is.na(sub.exon$five_prime_utr.width)
#sub.exon$cs.pos[!logi.na]  <-  sub.exon$cs.pos[!logi.na] + sub.exon$five_prime_utr.width[!logi.na]
  
myData.gtf  <-  merge( sub.exon[,c("transcript_id","exon_number","cs.pos","exon.orient.start")] , myData.gtf , by = c("transcript_id","exon_number") , all.y = TRUE , fill = 0 )
  

  ############### 3' UTR ##############
  verb("\t\t3.\n")
  
  logi.utr  <-  myData.gtf$feature == "three_prime_utr"
  if (any(logi.utr)) {
    myData.gtf$cs.pos[logi.utr]  <-  myData.gtf$transcript.length[logi.utr] - myData.gtf$three_prime_utr.width[logi.utr]
  } # exists
  
  ######### 5' UTR # is already = 0 (by fill) 
  
  #### add info to attributes
  verb("\tadd info to attributes.\n")
  
  p.3  <-  paste( 'three_prime_utr_length "' , myData.gtf$three_prime_utr.width , '"; ' , sep="")
  p.5  <-  paste( 'five_prime_utr_length "' , myData.gtf$five_prime_utr.width , '"; ' , sep="")
  p.L  <-  paste( 'transcript_length "' , myData.gtf$transcript.length , '"; ' , sep="")
  
  myData.gtf$attribute  <-  paste( myData.gtf$attribute , " " , p.3 , p.5 , p.L , sep = "")
  
  
  # feature = "transcript","three_prime_utr","five_prime_utr" have "NA" for exon_number and hence
  # will lack exon.cs.   This is easy to infer though
  
  
  ########## to transcript coords.
  verb("\tto transcript coords.\n")
  

################################################################################
####################### ************ 5' UTR ************  ######################
  verb("\t\t5.\n")
  
  sub.feat  <-  myData.gtf[ myData.gtf$feature == "five_prime_utr" , ,drop=FALSE]
  if (nrow(sub.feat) > 0) {
    logi.dup  <-  duplicated( sub.feat[, c("transcript_id") ])
    sub.feat  <-  sub.feat[ !logi.dup ,]
    
    sub.feat$seqname  <-  sub.feat$transcript_id
    sub.feat$strand  <-  "+"
    sub.feat$start  <-  1
    sub.feat$end  <-  sub.feat$five_prime_utr.width
    
    myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
  } # has 5 utr
  
  
  
  ### 3' UTR
  verb("\t\t3.\n")
  
  sub.feat  <-  myData.gtf[ myData.gtf$feature == "three_prime_utr" , ,drop=FALSE]
  if (nrow(sub.feat) > 0) {
    logi.dup  <-  duplicated( sub.feat[, c("transcript_id") ])
    sub.feat  <-  sub.feat[ !logi.dup ,]
    
    sub.feat$seqname  <-  sub.feat$transcript_id
    sub.feat$strand  <-  "+"
    sub.feat$start  <-  sub.feat$transcript.length - sub.feat$three_prime_utr.width + 1
    sub.feat$end  <-  sub.feat$transcript.length
    
    myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
  } # has 3 utr
  

  
  
  ###### start codon
  verb("\tstart codon.\n")
  
  sub.feat  <-   myData.gtf[ myData.gtf$feature == "start_codon" , ,drop=FALSE]
  sub.feat$seqname  <-  sub.feat$transcript_id
  sub.feat$strand  <-  "+"
  sub.feat$start  <-  sub.feat$five_prime_utr.width + 1
  sub.feat$end  <-  sub.feat$five_prime_utr.width + 3
  sub.feat  <-  sub.feat[ !duplicated(sub.feat[,c("feature","transcript_id")])  ,]
  
  myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
  
  
  ###### stop codon
  verb("\tstop codon.\n")
  
  sub.feat  <-   myData.gtf[ myData.gtf$feature == "stop_codon" , ,drop=FALSE]
  sub.feat$seqname  <-  sub.feat$transcript_id
  sub.feat$strand  <-  "+"
  sub.feat$start  <-  sub.feat$transcript.length - sub.feat$three_prime_utr.width - 2
  sub.feat$end  <-  sub.feat$transcript.length - sub.feat$three_prime_utr.width
  sub.feat  <-  sub.feat[ !duplicated(sub.feat[,c("feature","transcript_id")])  ,]
  
  myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
  
  
  ###### CDS
  verb("\t\tCDS.\n")
  
  tid.start  <-   myData.gtf$transcript_id[ myData.gtf$feature == "start_codon" ]
  tid.stop  <-   myData.gtf$transcript_id[ myData.gtf$feature == "stop_codon" ]
  tid.both  <-  intersect( tid.start , tid.stop )
  
  sub.feat  <-  myData.gtf[ myData.gtf$feature == "exon" , ,drop=FALSE]
  logi.dup  <-  duplicated(sub.feat[, c("transcript_id")])
  sub.feat  <-  sub.feat[ !logi.dup ,]
  sub.feat  <-  sub.feat[ sub.feat$transcript_id  %in%  tid.both , ,drop = FALSE ]
  
  sub.feat$seqname  <-  sub.feat$transcript_id
  sub.feat$strand  <-  "+"
  sub.feat$feature  <-  "CDS"
  sub.feat$start  <-  sub.feat$five_prime_utr.width + 1
  sub.feat$end  <-  sub.feat$transcript.length - sub.feat$three_prime_utr.width
  
  myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
  
  
  ####### exon
  verb("\t\texon.\n")
  # Note:  the definition of exon is anything that is transcribed and retained after splicing.
  # Thus, the 5' and 3' UTR are considered exons.
  # Exons are NOT restricted to just the coding sequence.
  #
  # Thus, the following gets all the non protein coding genes as well (things without 5'/3' UTR and CDS annotations)
  
  sub.feat  <-   myData.gtf[ myData.gtf$feature == "exon" , ,drop=FALSE]
  sub.feat$seqname  <-  sub.feat$transcript_id
  sub.feat$strand  <-  "+"
  sub.feat$start  <-  sub.feat$cs.pos + 1
  sub.feat$end  <-  sub.feat$start + sub.feat$width - 1
  sub.feat  <-  sub.feat[ !duplicated(sub.feat[,c("feature","transcript_id","exon_number"),drop=FALSE])  , ,drop=FALSE]
  
  myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
  
  
  ####################### sanity checks
  verb("\t\tsanity check.\n")
  
  tdf  <-  myData.tranGTF  
  tdf$transcript_id <- get_gtf_attribute_field( gtf = tdf , field = "transcript_id" , num = FALSE )
  tdf$gene_name <- get_gtf_attribute_field( gtf = tdf , field = "gene_name" , num = FALSE )
  tdf$gene_id <- get_gtf_attribute_field( gtf = tdf , field = "gene_id" , num = FALSE )
  tdf$transcript_length  <-  get_gtf_attribute_field( gtf = tdf , field = "transcript_length" , num = TRUE )
  tdf$five_prime_utr_length  <-  get_gtf_attribute_field( gtf = tdf , field = "five_prime_utr_length" , num = TRUE )
  tdf$three_prime_utr_length  <-  get_gtf_attribute_field( gtf = tdf , field = "three_prime_utr_length" , num = TRUE )
  tdf$inferred.CDS.length  <-  tdf$transcript_length - (tdf$five_prime_utr_length + tdf$three_prime_utr_length)
  tdf$CDS_length  <-  tdf$end - tdf$start + 1
  tdf$gene_biotype <- get_gtf_attribute_field( gtf = tdf , field = "gene_biotype" , num = FALSE )
  tdf$transcript_biotype <- get_gtf_attribute_field( gtf = tdf , field = "transcript_biotype" , num = FALSE )
  
  ##### UTR length matches
  verb("\t\t\tUTR lengths agree.\n")
  
  for (featx  in  c("five_prime_utr","three_prime_utr")) {
    subdf  <-  tdf[ tdf$feature == featx , ,drop=FALSE]
    logi.eq  <-  subdf$CDS_length == subdf[[paste(featx,"length",sep="_")]]
    if (!all(logi.eq)) {
      verb("\n\n\n\nERROR!  not all [%s] lengths agree!\n\n\n\n" , featx)
      stop()
    }
  } # featx
  
  
  ##### inferred CDS length is real
  verb("\t\t\tinferred CDS length is real.\n")
  
  subdf  <-  tdf[ tdf$feature == "CDS" ,]
  logi.infer  <-  subdf$CDS_length == subdf$inferred.CDS.length
  if (!all(logi.infer)) {
    verb("\n\n\n\nERROR!  inferred CDS length != annotated length!!!!\n\n\n")
    stop()
  } # infer
  
  
  #### plots
  verb("\t\tplots.\n")
  
  fname  <-  sprintf("%s.pdf", outfname )
  pdf( file = fname )
  
  feats  <-  c("transcript_length" , "five_prime_utr_length" , "three_prime_utr_length" , "CDS_length" )
  for (featx  in  feats) {
    gp  <-  ggplot( data = tdf , aes_string( featx ) )   +
      geom_histogram( bins = 200 )   +
      theme_bw()   +
      ggtitle(sprintf("Histogram of %s", featx))
    print(gp)
  } # featx
  
  dummydev  <-  dev.off()
  
  
  ########### sort
  verb("\tsort.\n")
  
  myData.tranGTF  <-  myData.tranGTF[ order(myData.tranGTF$seqname , myData.tranGTF$start , myData.tranGTF$end) , ]
  return(tdf)
  
  
  ##### save
  verb("\tsave.\n")
  
  fname  <-  sprintf("%s", outfname)
  write.table( myData.tranGTF , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
  
} # make_eukaryotic_transcript_gtf


Mus_musculus.GRCm39.gtf <- make_eukaryotic_transcript_gtf(gtf = gtf, outfname = "Mus_musculus.GRCm39.104.txt")
#rm(list=setdiff(ls(), "Mus_musculus.GRCm39.gtf"))
save(Mus_musculus.GRCm39.gtf, file = "Mus_musculus.GRCm39.gtf.RData")

readr::write_tsv(x = Mus_musculus.GRCm39.gtf, file = "Mus_musculus.GRCm39.104.rdna_rn18s.transcriptome.gtf", col_names = FALSE, quote = "none", escape = "none")
# Transfer file on HPC and use dos2unix and make .gz file

################################################################################
###### ********** This part will make BED file for UTR & CDS ********** ########

# Convert GTF file in BED for sequence retrieval for HOMER motif analysis.
# This BED file will be used for making fasta file of UTRs, CDS.
# We have to use start-1 for start position, as BED start with 0
# https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/


Mus.musculus.GRCm39.BED <- Mus_musculus.GRCm39.gtf %>%
  dplyr::mutate(name=seqname) %>%
  dplyr::select(seqname, start, end, name, feature, strand) %>%
  dplyr::mutate(start=start-1) %>%
  readr::write_tsv(file = "Mus_musculus.GRCm39.BED", col_names = FALSE)

####################### This is the end of UTR BED files #######################


################################################################################
###### ********** This part will make BED file for 5'UTR-CDS ********** ########

# In this code I already added function "five_prime_utr_to_5UTR_CDS", which will make two BED file along with "Mus_musculus.GRCm39.BED". 
# Details are in below lines -----
# This part have issue, if 5'UTR is smaller than 20 then start will be negative. Matt Park use only 20 bp of CDS if 5'UTR is smaller than 20.
# Here I make two BED if 5'UTR is smaller tahn 20)
# 1. Mus_musculus.GRCm39.Matt.5UTR-CDS.bed :: Similar to Matt, if 5'UTR is smaller than 20 BP, ignore UTR use used 20 BP from CDS
# 2. Mus_musculus.GRCm39.Nitish.5UTR-CDS.bed :: Modified, if 5'UTR is smaller than 20 BP still use all 5'UTR sequences with 20 BP CDS in BED.

five_prime_utr_to_5UTR_CDS = function(data, All_UTR_Bases=TRUE) {
  data <- data %>%
    dplyr::filter(feature=="five_prime_utr") %>%
    dplyr::mutate(feature = gsub("five_prime_utr","five_prime_utr_CDS", feature))
  
  if(All_UTR_Bases==TRUE)
  {
    data %>%
      dplyr::select(c( "seqname", "start", "end", "name", "feature", "strand" )) %>%
      dplyr::rename(Chr=seqname, Start=start, End=end, Gene=name, Feature=feature, Strand=strand) %>%
      dplyr::mutate(Start=End-20) %>%
      dplyr::mutate(End=End+20) %>%
      dplyr::mutate(Start = if_else(Start < 0, 0, Start)) %>%
      readr::write_tsv(file = "Mus_musculus.GRCm39.Nitish.5UTR-CDS.bed", col_names = FALSE)
  }
  else
  {
    data %>%
      dplyr::select(c( "seqname", "start", "end", "name", "feature", "strand" )) %>%
      dplyr::rename(Chr=seqname, Start=start, End=end, Gene=name, Feature=feature, Strand=strand) %>%
      dplyr::mutate(Start=End-20) %>%
      dplyr::mutate(End=End+20) %>%
      dplyr::mutate(Start = if_else(Start < 0, End-20, Start)) %>%
      readr::write_tsv(file = "Mus_musculus.GRCm39.Matt.5UTR-CDS.bed", col_names = FALSE)
  }
}

five_prime_utr_to_5UTR_CDS.bed <- five_prime_utr_to_5UTR_CDS(data = Mus.musculus.GRCm39.BED, All_UTR_Bases = TRUE)
five_prime_utr_to_5UTR_CDS.1.bed <- five_prime_utr_to_5UTR_CDS(data = Mus.musculus.GRCm39.BED, All_UTR_Bases = FALSE)

################################################################################
################################################################################
