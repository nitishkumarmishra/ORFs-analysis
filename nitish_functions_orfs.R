# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2023
# Email:  nitishimtech@gmail.com
# Date: 03 March, 2023
# Script Name: code_name.R
#####################################################################
# ##################  Script Description: ###########################
#
##################### Notes: ########################################
#
#####################################################################
################## SET WORKING DIRECTORY ############################

# R functions for the ORF read count analysis

library(tidyverse)
library(Biostrings)
library(plyr)
library(dplyr)
library(stringr)
library(limma)
library(cqn)
library(edgeR)
library(methods)
library(utils)
library(stats)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(Matrix)  # rank
library(gplots)
library(cowplot)
library(scales) # muted
library(GGally) # scatmat
library(data.table)
library(AnnotationDbi)

###############################################################################
###############################################################################

suppressMessages(suppressWarnings(source("/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/data/mgi/jupyter_common.R")))

################################################################################
################################################################################

# verb
verb <- function(...) cat(sprintf(...), sep='', file=stdout())

################################################################################
################################################################################


# get transcriptome file

get_transcriptome_feature_fasta_filename  <-  function( species , feat ) {
  
  #items <- strsplit(feat, "_")[[1]]
  items <- feat
  for (i in length(items):1) {
    
    file  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/m39/Mus_musculus.GRCm39.104.rdna_rn18s.transcriptome.gtf.gz"
    basef  <-  gsub( "\\.gtf\\.gz$", "" , file)
    #featf  <-  paste( basef , feat , "fa.gz" , sep=".")
    featf  <-  paste( basef , paste(items[1:i], collapse='_') , "fa.gz" , sep=".")
    
    if (file.exists(featf)) break
  }
  
  return(featf)
  
} # get_transcriptome_feature_fasta_filename

################################################################################
################################################################################

read_gtf  <-  function( file , header = FALSE ) {
  
  in.gtf  <-  fread( cmd = sprintf("zcat  -f  %s", file) , sep = "\t" , header = FALSE ,  stringsAsFactors = FALSE , quote = "" , colClasses = c( rep("character",3) , rep("numeric",2) , rep("character",4) ) )
  
  def.colnames  <-  c("seqname" , "source" , "feature" , "start" , "end" , "score" , "strand" , "frame" , "attribute")
  colnames(in.gtf)  <-  def.colnames
  
  in.gtf  <-  as.data.frame(in.gtf)
  
  return(in.gtf)
  
} # read_gtf



################################################################################
################################################################################

load_transcriptome_gtf  <-  function( species , verbose = FALSE ) {
  
  ### transcriptome
  if (species == "m39") {
    fname  <-  sprintf("/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/m39/Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz")
  } else {
    verb("\n\n\nERROR!  unrecognized species=[%s]!!!\n", species)
    stop()
  }
  if (verbose) { verb("\t\t\tloading transcriptome gtf: [%s]\n", fname ) }
  
  indf  <-  read_gtf( file = fname )
  
  return(indf)
  
} # load_transcriptome_gtf



################################################################################
################################################################################

get_gtf_attribute_field  <-  function( gtf , field , num = FALSE ) {
  # always returns as a string
  # num = TRUE -> will convert to numeric (via as.numeric)
  
  library(stringr)
  
  fmt  <-  paste( field , ' "[^"]+' , sep="" )
  exf  <-  str_extract( gtf$attribute , fmt )
  
  fmt  <-  paste( field , ' "' , sep="")
  exf  <-  str_replace( exf , fmt , '' )
  
  if (num) { exf  <-  as.numeric(exf) }
  return(exf)
  
} # get_gtf_attribute_field


################################################################################
################################################################################


prepare_transcriptome_gtf  <-  function( gtf ) {
  
  if (!("transcript_id" %in% colnames(gtf))) {
    gtf$transcript_id  <-  get_gtf_attribute_field( gtf = gtf , field = "transcript_id" )
  }
  if (!("gene_name" %in% colnames(gtf))) {
    gtf$gene_name  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_name" )
  }
  if (!("gene_id" %in% colnames(gtf))) {
    gtf$gene_id  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_id" )
  }
  
  return(gtf)
  
} # prepare_transcriptome_gtf

################################################################################
################################################################################



extract_and_append_attribute_field_to_gtf  <-  function( gtf , field , num = FALSE , overwrite = FALSE ) {
  
  if (overwrite  ||  !(field %in% colnames(gtf))) {
    gtf[[field]]  <-  get_gtf_attribute_field( gtf = gtf , field = field , num = num )
  } # write
  
  return(gtf)
  
} # extract_and_append_attribute_field_to_gtf



################################################################################
################################################################################

# matt_get_info  <-  function( feat , species = NULL , strip = TRUE ) {
#   
#   ######## prepare
#   general  <-  is.null(species)
#   
#   if (!general) {
#     if (species == "mouse") { species  <-  "mm10" }
#   } # mouse
#   
#   ######## make command
#   shcmd  <-  sprintf("/home/map2085/riboprof/src/matt_get_info.sh  -f %s " , feat )
#   if (general) {
#     shcmd  <-  sprintf("%s  -g ", shcmd )
#   }
#   if (!is.null(species)) {
#     shcmd  <-  sprintf("%s  -s %s", shcmd , species)
#   }
#   if (!strip) {
#     shcmd  <-  sprintf("%s  -z " , shcmd)
#   }
#   
#   
#   ####### execute
#   val  <-  system( command = shcmd , intern = TRUE )
#   
#   return(val)
#   
# } # matt_get_info



################################################################################
################################################################################

get_project_info  <-  function( project = NULL , individual = NULL , runid = NULL , strategy = NULL , verbose = FALSE ) {
  
  
  
  dbdir  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/data2/database"
  poss.names  <-  list.files( path = dbdir , pattern = "*.txt" )
  poss.names  <-  poss.names[!(poss.names %in% c("sample.db.txt","sample.db.notes.txt","public.db.txt","bad.1000G.txt",
                                                 "general.database.txt","organism.databases.txt","bad.hipsci.txt" , "inhouse.txt","ecoliOperonExpr.txt" ,
                                                 "protoVariant.database.txt" , "OLD.PRJEB9586.txt"))]
  
  
  db  <-  paste( dbdir , poss.names , sep = "/" )
  
  if (verbose) { show(db) }
  
  indf  <-  data.frame()
  for (dbx  in  db) {
    if (verbose) { verb("\t\tdbx = [%s]\n", dbx ) }
    
    invardb  <-  read.table( file = dbx , header = TRUE , sep = "\t" ,  row.names = NULL , stringsAsFactors = FALSE , quote = "" , fill = TRUE , comment.char = "" )
    
    invardb$project  <-  as.character(invardb$project)
    invardb$individual  <-  as.character(invardb$individual)
    invardb$runid  <-  as.character(invardb$runid)
    
    ### find
    sdf  <-  data.frame()
    if (!is.null(runid)) {	
      sdf  <-  invardb[ invardb$runid %in% runid ,,drop=FALSE]
    }
    if (!is.null(individual)) {
      sdf  <-  invardb[ invardb$individual %in% individual ,,drop=FALSE]
    }
    if (!is.null(project)) {
      sdf  <-  invardb[ invardb$project %in% project ,,drop=FALSE]
    } # project
    
    if (verbose) { verb("\t\t\t%d rows\n", nrow(sdf) ) }
    
    if (nrow(sdf) > 0) {
      indf  <-  rbind.fill( indf , sdf )
    } # sdf
  } # dbx
  
  
  if (!is.null(strategy)) {
    stopifnot("strategy"  %in%  colnames(indf))
    indf  <-  indf[ indf$strategy %in% strategy  &  !is.na(indf$strategy) ,, drop=FALSE]
  } # strategy
  
  
  
  ### check
  if (nrow(indf) == 0) {
    verb("\n\n\nERROR!  unable to find entries!!!\n" )
    show(project)
    show(runid)
    stop()
  } # check
  
  return(indf)
  
} # get_project_info



################################################################################
################################################################################

process_project_input  <-  function( project = NULL , runid = NULL , individual = NULL , strategy = NULL , verbose = FALSE ) {
  
  if ( (is.character(project)  ||  is.character(runid))  &&  !is.data.frame(project)  && !is.data.table(project)  ) {
    project  <-  get_project_info( project = project , individual = individual , runid = runid , strategy = strategy , verbose = verbose )
  }
  
  stopifnot(!any(is.na(project$condition)))
  stopifnot(!any(is.na(project$runid)))
  stopifnot(!any(is.na(project$individual)))
  stopifnot(!any(is.na(project$strategy)))
  stopifnot(!any(is.na(project$layout)))
  stopifnot(!any(is.na(project$species)))
  
  if ("order" %in% colnames(project)) {
    stopifnot(!any(is.na(project$order)))
    
    if (is.data.table(project)) {
      setkey(project , order , runid)
    } else {
      project  <-  project[ order(project$order , project$runid) , ,drop=FALSE]
    } 
    
    ### each condition must have a distinct order number, with a single order number per condition.
    if ("condition" %in% colnames(project)) {
      uco  <-  unique(project[,c("condition","order"),drop=FALSE])
      #stopifnot(!any(duplicated(uco$condition)))
      stopifnot(!any(duplicated(uco$order)))
    } # condition
  } # order
  
  rownames(project)  <-  project$runid
  
  return(project)
  
} # process_project_input



################################################################################
################################################################################

load_accumulated_RPF_count_per_ORF  <-  function( project , classAccum , verbose = FALSE ) {
  
  project  <-  process_project_input( project = project )
  
  rdf  <-  data.frame()
  for (rx  in  1:nrow(project)) {
    projx  <-  project$project[rx]
    indivx  <-  project$individual[rx]
    runx  <-  project$runid[rx]
    
    if (verbose) { verb("\t\t%s\n", runx) }
    
    bdir  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM"
    bdir  <-  sprintf("%s/%s", bdir, runx)
    
    if (classAccum) {
      fname  <-  sprintf("%s.class.bed" , bdir)
      indf  <-  read.table( file = fname , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "" , colClasses = c("character" , rep("numeric",2) , "character" , "numeric"  )  )
      colnames(indf)  <-  c("transcript_id","start","end", "classification" , "count")
      
    } else {
      fname  <-  sprintf("%s.bed" , bdir)
      indf  <-  read.table( file = fname , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "" , colClasses = c("character" , rep("numeric",3))  )
      colnames(indf)  <-  c("transcript_id","start","end", "classification" , "count")
    } # classAccum
    
    indf$runid  <-  runx
    
    rdf  <-  rbind( indf , rdf )
  } # rx
  
  return(rdf)
  
} # load_accumulated_RPF_count_per_ORF


################################################################################
################################################################################



#' Find all instance of a pattern in a set of sequences
#'
#' This function searched a set of sequences for a fixed string pattern and reports
#' start and end positions (both 1-based) for all matches of pattern in all sequences.
#'
#' @param  seq		DNAStringSet of sequences from Biostrings
#' @param  pattern	a character string. Wildcards not in regex format, but in IUPAC format.
#'			Will be passed to "vmatchPattern" from Biostrings
#' @param  name		name of varialbe containg sequence names.  e.g. "chr", "transcript_id", etc.
#'
#' @return	data.frame with columns "ix","start","end",name,"pattern"
#'
#' @export
find_all_pattern_indices_in_string_set  <-  function( seq , pattern , name = "transcript_id" ) {
  
  ### get hits
  mh1  <-  vmatchPattern(pattern = pattern , subject = seq )
  
  ### long format
  hdf  <-  melt(mh1) # melt works on this class, miraculously.
  colnames(hdf)  <-  c("ix","NA","start","end","width")
  hdf  <-  hdf[, c("ix","start","end") ,drop=FALSE]
  
  # name correctly
  tnames  <-  data.frame( ix = 1:length(names(mh1)) , transcript_id = names(mh1) )
  colnames(tnames)[2]  <-  name
  
  hdf  <-  merge( hdf , tnames , by = "ix" , all.x = TRUE )
  hdf$pattern  <-  pattern
  
  return(hdf)
  
} # find_all_pattern_indices_in_string_set




################################################################################
################################################################################


#' Test for differences in RPF ratios between ORF types for specific gene sets across conditions.
#'
#' Given a gene set of interest and RPF data for various conditions, for each ORF type classification C
#' (e.g. uORF, truncated ORF, etc.), this function will test if the ratio of total RPF counts for that  
#' ORF classification over the total RPF count for the annotated, canonical ORFs changes between conditions,
#' where the total RPF counts are restricted to the gene set of interest.
#'
#' @param  orfs		data.frame output of function "find_all_ORFs_in_transcriptome". Loaded from file if missing.
#' @param  rpf		data.frame output of function "accumulate_rpf_counts_per_putative_ORFs.sh".  Loaded from file if missing.
#' @param  species	character string giving species name.  used to load orfs and rpf if they are missing.
#' @param  groups	data frame with columns "genes" and "group"
#' @param  project	project info
#'
#' @return	data.frame of hypergeometric test results.
#'
#' @export
test_ORF_classification_read_count_between_conditions  <-  function( gtf = NULL , rpf = NULL , species = NULL , groups , project , classes = NULL , outfbase = NULL , verbose = FALSE ) { 
  
  # DEBUG
  #fname  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\ codes/out/limma-voom.mrna/161021.rdna_rn18s/NMuMG/blancgrp_161021_Riboseq.unt48--vs--tgfbCX5461.diff-all.txt.gz"
  #groups  <-  read.table( file = fname , header = TRUE , row.names = 1 , sep="\t" , stringsAsFactors = FALSE , quote = "" )
  #groups$gene_name  <-  rownames(groups)
  #groups$group  <-  ifelse( groups$log2FC > 0 , "up" , "down" )
  
  
  project  <-  process_project_input( project = project )
  
  groups  <-  unique(groups)
  
  ### load
  if (verbose) { verb("\tload.\n") }
  
  
  ### rpf
  if (is.null(rpf)) {
    rpf  <-  load_accumulated_RPF_count_per_ORF( project = project , classAccum = TRUE , verbose = verbose )
  } # rpf
  stopifnot(is.data.frame(rpf))
  
  if (!("condition" %in% colnames(rpf))) {
    if (verbose) { verb("\t\tcond.\n") }
    rpf  <-  merge( rpf , unique(project[,c("runid","condition"),drop=FALSE]) , by = "runid" )
  }
  
  ### gtf
  if (is.null(gtf)) {
    if (verbose) { verb("\t\tload gtf.\n") }
    gtf  <-  load_transcriptome_gtf( species = species )
    gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
  } # proteinCodingTranscripts
  
  ### prepare gtf
  if (verbose) { verb("\t\tprepare gtf.\n") }
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_name" ) # Here I am usig "gene_name" rather than "gene" in GTF file
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_id" )
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_biotype" )
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_biotype" )
  
  
  
  ### test each group
  if (verbose) { verb("\ttest each group.\n") }
  
  conpairs  <-  combn( unique(project$condition) , 2)
  
  uniq.classes  <-  unique(rpf$classification)
  uniq.classes  <-  uniq.classes[!is.na(uniq.classes)]
  if (!is.null(classes)) {
    if (verbose) { verb("\trestricting to classes.\n") }
    uniq.classes  <-  intersect(uniq.classes , classes)
  } # 
  show(uniq.classes)
  
  pctids  <-  unique(gtf$transcript_id[gtf$gene_biotype == "protein_coding"  &  gtf$transcript_biotype == "protein_coding"])
  
  resdf  <-  data.frame()
  for (groupx  in  unique(groups$group)) {
    if (verbose) { verb("\t\t%s\n", groupx) }
    
    subgenes  <-  unique(groups$gene_name[groups$group == groupx])
    subtids  <-  unique(gtf$transcript_id[gtf$gene_name  %in%  subgenes])
    
    if (length(subgenes) < 2) { next }
    
    ### summarize by class, protein coding
    pcdf  <-  rpf[ rpf$transcript_id %in% pctids  &  rpf$transcript_id %in% subtids ,,drop=FALSE]
    pcdf  <-  pcdf  %>%  group_by(classification,runid,condition)  %>%  summarise(total.count = sum(count))
    pcdf  <-  as.data.frame(pcdf)
    
    anndf  <-  pcdf[ pcdf$classification == "annotated ORF"  &  !is.na(pcdf$classification) ,,drop=FALSE]
    anndf$classification  <-  NULL
    colnames(anndf)[colnames(anndf) == "total.count"]  <-  "annotated.ORF.count"
    pcdf  <-  merge( pcdf , anndf , by = c("runid","condition") , all = TRUE )
    pcdf$log.odds  <-  log2(pcdf$total.count) - log2(pcdf$annotated.ORF.count)
    
    
    for (classx  in  uniq.classes) {
      if (classx == "annotated ORF") { next }
      if (verbose) { verb("\t\t\t%s\n", classx) }
      
      for (cx  in  1:ncol(conpairs)) {
        cond1  <-  conpairs[1,cx]
        cond2  <-  conpairs[2,cx]
        if (verbose) { verb("\t\t\t\t%s  %s\n", cond1, cond2) }
        
        subdf  <-  pcdf[ pcdf$classification == classx  &  !is.na(pcdf$classification)  &  pcdf$condition %in% c(cond1,cond2) ,,drop=FALSE]
        if (nrow(subdf) < 2  ||  length(unique(subdf$condition)) < 2) { next }
        
        subdf$condition  <-  factor( subdf$condition , levels = c(cond1,cond2) )
        
        if (any(is.na(subdf$log.odds))) { next }
        if (any(!is.finite(subdf$log.odds))) { next }
        
        # test	
        tres  <-  t.test( formula = log.odds ~ condition , data = subdf )
        pval  <-  tres$p.value
        
        # summarize
        sdf  <-  subdf  %>%  group_by(condition)  %>% summarise(mean.log.odds = mean(log.odds) , sd.log.odds = sd(log.odds) )
        sdf  <-  as.data.frame(sdf)
        rownames(sdf)  <-  sdf$condition
        
        rdf  <-  data.frame( classification = classx , group = groupx , condition1 = cond1 , condition2 = cond2 ,
                             num.genes = length(unique(subgenes)) , p.value = pval ,
                             mean.log.odds.1 = sdf[cond1,"mean.log.odds"] , sd.log.odds.1 = sdf[cond1,"sd.log.odds"] ,
                             mean.log.odds.2 = sdf[cond2,"mean.log.odds"] , sd.log.odds.2 = sdf[cond2,"sd.log.odds"]  )
        resdf  <-  rbind( resdf , rdf )
      } # cx
    } # classx
  } # groupx
  
  resdf  <-  resdf[order(resdf$classification , resdf$p.value) , ,drop=FALSE]
  resdf$log2FC  <-  (resdf$mean.log.odds.2 - resdf$mean.log.odds.1)
  
  if (!is.null(outfbase)) {
    fname  <-  sprintf("%s.test.txt" , outfbase )
    write.table( resdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
    system(sprintf("gzip  -f  %s" , fname ))
  } # outfbase
  
  return(resdf)	
  
} # test_ORF_classification_read_count_between_conditions



################################################################################
################################################################################


#' Test for differences in RPF ratios between ORF types for specific gene sets across conditions.
#'
#' Given a gene set of interest and RPF data for various conditions, for each ORF type classification C
#' (e.g. uORF, truncated ORF, etc.), this function will test if the ratio of total RPF counts for that  
#' ORF classification over the total RPF count for the annotated, canonical ORFs changes between conditions,
#' where the total RPF counts are restricted to the gene set of interest.
#'
#' @param  rpf		data.frame output of function "accumulate_rpf_counts_per_putative_ORFs.sh".  Loaded from file if missing.
#' @param  species	character string giving species name.  used to load orfs and rpf if they are missing.
#' @param  project	project info
#'
#' @return	data.frame of hypergeometric test results.
#'
#' @export

test_ORF_classification_read_count_between_conditions_per_gene  <-  function( gtf = NULL , rpf = NULL , species = NULL , project , outfbase = NULL , verbose = FALSE ) { 
  
  project  <-  process_project_input( project = project )
  
  ### load
  if (verbose) { verb("\tload.\n") }
  
  ### rpf
  if (is.null(rpf)) {
    rpf  <-  load_accumulated_RPF_count_per_ORF( project = project , classAccum = TRUE , verbose = verbose )
  } # rpf
  stopifnot(is.data.frame(rpf))
  
  if (!("condition" %in% colnames(rpf))) {
    if (verbose) { verb("\t\tcond.\n") }
    rpf  <-  merge( rpf , unique(project[,c("runid","condition"),drop=FALSE]) , by = "runid" )
  }
  
  ### gtf
  if (is.null(gtf)) {
    if (verbose) { verb("\t\tload gtf.\n") }
    gtf  <-  load_transcriptome_gtf( species = species )
    gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
  } # proteinCodingTranscripts
  
  ### prepare gtf
  if (verbose) { verb("\t\tprepare gtf.\n") }
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_name" ) # Here I am usig "gene_name" rather than "gene" in GTF file
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_id" )
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_biotype" )
  gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_biotype" )
  colnames(gtf) <- gsub("gene_name", "gene",colnames(gtf))
  
  #### merge gene name
  if (verbose) { verb("\tmerge gene name.\n") }
  
  rpf  <-  merge( rpf , unique(gtf[,c("gene","transcript_id"),drop=FALSE]) , all.x = TRUE )
  
  conpairs  <-  combn( unique(project$condition) , 2)
  
  uniq.classes  <-  unique(rpf$classification)
  uniq.classes  <-  uniq.classes[!is.na(uniq.classes)]
  show(uniq.classes)
  
  pctids  <-  unique(gtf$transcript_id[gtf$gene_biotype == "protein_coding"  &  gtf$transcript_biotype == "protein_coding"])
  
  resdf  <-  data.frame()
  
  ### summarize by gene, class, protein coding
  pcdf  <-  rpf[ rpf$transcript_id %in% pctids  ,,drop=FALSE]
  pcdf  <-  pcdf  %>%  group_by(gene,classification,runid,condition)  %>%  summarise(total.count = sum(count))
  pcdf  <-  as.data.frame(pcdf)
  
  anndf  <-  pcdf[ pcdf$classification == "annotated ORF"  &  !is.na(pcdf$classification) ,,drop=FALSE]
  anndf$classification  <-  NULL
  colnames(anndf)[colnames(anndf) == "total.count"]  <-  "annotated.ORF.count"
  pcdf  <-  merge( pcdf , anndf , by = c("gene","runid","condition") , all = TRUE )
  pcdf$log.odds  <-  log2(pcdf$total.count) - log2(pcdf$annotated.ORF.count)
  logi.zn  <-  pcdf$total.count == 0  |  is.na(pcdf$total.count)
  logi.zd  <-  pcdf$annotated.ORF.count == 0  |  is.na(pcdf$annotated.ORF.count)
  pcdf$log.odds[logi.zn & !logi.zd]  <-  log2(pcdf$total.count[logi.zn & !logi.zd] + 1e-6) - log2(pcdf$annotated.ORF.count[logi.zn & !logi.zd])
  pcdf$log.odds[!logi.zn & logi.zd]  <-  log2(pcdf$total.count[!logi.zn & logi.zd]) - log2(pcdf$annotated.ORF.count[!logi.zn & logi.zd] + 1e-6 )
  
  
  for (classx  in  uniq.classes) {
    if (classx == "annotated ORF") { next }
    if (verbose) { verb("\t\t\t%s\n", classx) }
    
    for (cx  in  1:ncol(conpairs)) {
      cond1  <-  conpairs[1,cx]
      cond2  <-  conpairs[2,cx]
      if (verbose) { verb("\t\t\t\t%s  %s\n", cond1, cond2) }
      
      subdf  <-  pcdf[ pcdf$classification == classx  &  !is.na(pcdf$classification)  &  pcdf$condition %in% c(cond1,cond2) ,,drop=FALSE]
      if (nrow(subdf) < 2  ||  length(unique(subdf$condition)) < 2) { next }
      
      
      ### matrix
      if (verbose) { verb("\t\t\t\tmatrix.\n") }
      
      gmat  <-  reshape2::acast( data = subdf , formula = gene ~ runid , value.var = "log.odds" )
      gmat  <-  gmat[complete.cases(gmat),,drop=FALSE]
      if (nrow(gmat) == 0) { next }
      
      ridc1  <-  project$runid[project$condition == cond1]
      logi.c1  <-  colnames(gmat) %in% ridc1
      tres  <-  apply( gmat , 1 , FUN = function(r) t.test( x = r[logi.c1] , y = r[!logi.c1] , alternative = "two.sided" )$p.value )
      tres  <-  data.frame( gene = rownames(gmat) , p.value = tres )
      
      tcdf  <-  reshape2::dcast( data = subdf , formula = gene ~ runid , value.var = "total.count" , fill = 0 )
      aodf  <-  reshape2::dcast( data = subdf , formula = gene ~ runid , value.var = "annotated.ORF.count" , fill = 0 )
      lodf  <-  reshape2::dcast( data = subdf , formula = gene ~ runid , value.var = "log.odds" , fill = 0 )
      colnames(lodf)  <-  paste( colnames(lodf) , "log.odds" , sep="." )
      
      
      ### summarize
      if (verbose) { verb("\t\t\t\tsummarize.\n") }
      
      mdf  <-  merge( tcdf , aodf , by = "gene" , suffixes = c(sprintf(".%s",classx) , ".CDS" ) , all = TRUE )
      mdf  <-  merge( mdf , lodf , by.x = "gene" , by.y = "gene.log.odds" , all = TRUE )
      mdf  <-  merge( mdf , tres , by = "gene" , all = TRUE )
      
      
      ### mean, sd
      if (verbose) { verb("\t\t\t\tmean, sd.\n") }
      
      sdf  <-  subdf  %>%  group_by(gene,classification,condition)  %>%  summarise(mean.log.odds = mean(log.odds) , sd.log.odds = sd(log.odds) )
      sdf  <-  as.data.frame(sdf)
      modf  <-  reshape2::dcast( data = sdf , formula = gene ~ condition , value.var = "mean.log.odds")
      modf$log2FC  <-  modf[[cond2]] - modf[[cond1]]
      sodf  <-  reshape2::dcast( data = sdf , formula = gene ~ condition , value.var = "sd.log.odds")
      
      smdf  <-  merge( modf , sodf , by = "gene" , suffixes = c(".mean",".sd") , all = TRUE )
      mdf  <-  merge( mdf , smdf , by = "gene" , all = TRUE )
      
      
      if (!is.null(outfbase)) {
        ### write
        if (verbose) { verb("\t\t\t\twrite.\n") }
        
        fname  <-  sprintf("%s.%s.%s--%s.summary.txt", outfbase , classx , cond1 , cond2 )
        write.table( mdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
      } # outfbase
      
      rdf  <-  mdf[,c("gene","log2FC","p.value"),drop=FALSE]
      rdf$classification  <-  classx
      rdf$condition1  <-  cond1
      rdf$condition2  <-  cond2
      resdf  <-  rbind( resdf , rdf )
    } # cx
  } # classx
  
  return(resdf)	
  
} # test_ORF_classification_read_count_between_conditions_per_gene



################################################################################
################################################################################


#'  Test enrichment of ORF types in gene sets
#'
#'  This function will perform hypergeometric tests to determine if
#'  genes with certain types of ORFs are overrepresented in any of the gene sets (groups).
#'
#' @param  orfs		data.frame output of function "find_all_ORFs_in_transcriptome"
#' @param  groups	data frame with columns "genes" and "group"
#'
#' @return	data frame with enrichment and depletion  p values from Hypergeometric test
#'		for each group, as well as summary statistics.
#'
#' @export
#' 
hypergeometric_test_ORF_type_enrichment  <-  function( orfs , groups , universe, gtf = NULL , species = NULL , verbose = FALSE, project="161021") {
  
  universe  <-  unique(universe)
  stopifnot(all(groups$gene_name  %in%  universe))
  
  ### load
  if (verbose) { verb("\tload.\n") }
  if (is.null(gtf)) {
    gtf  <-  load_transcriptome_gtf( species = species )
  } # gtf
  
  
  #### process
  if (verbose) { verb("\tprocess.\n") }
  gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
  gtf$gene_biotype  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_biotype" )
  gtf$transcript_biotype  <-  get_gtf_attribute_field( gtf = gtf , field = "transcript_biotype" )
  
  gtf_selected <- gtf |>
    dplyr::select(transcript_id, gene_name) |>
    na.omit() |>
    unique()
  if (is.null(orfs)) {
  orfs <- load_accumulated_RPF_count_per_ORF(project=project, classAccum=TRUE)}
  
  orfs <- merge(orfs , unique(gtf[,c("transcript_id","gene_name", "transcript_biotype")]) , by = "transcript_id" , all.x = TRUE )
  orfs <- orfs[!is.na(orfs$classification),]
  
  
  resdf  <-  data.frame()
  for (groupx  in  unique(groups$group)) {
    if (verbose) { verb("\t\t%s\n", groupx) }
    
    subdf  <-  groups[ groups$group == groupx , ,drop=FALSE]
    subdf  <-  unique(subdf)
    
    for (classx  in  unique(orfs$classification)) {
      if (verbose) { verb("\t\t\t%s\n", classx) }
      
      sodf  <-  orfs[ orfs$classification == classx  &  orfs$transcript_biotype == "protein_coding" , ,drop=FALSE]
      #sodf  <-  merge( sodf , gtf_selected[,c("transcript_id","gene_name")] , by = "gene_name" , all.x = TRUE ) # Here I am usig "gene_name" rather than "gene" in GTF file
      sodf  <-  sodf[ sodf$gene_name %in% universe , ,drop=FALSE]
      
      ### test
      if (verbose) { verb("\t\t\t\ttest\n") }
      
      x.samp.succ  <-  length(intersect(subdf$gene_name , sodf$gene_name))
      m.pop.succ  <-  length(unique(sodf$gene_name))
      n.pop.fail  <-  length(setdiff(universe , sodf$gene_name))
      k.samp.size  <-  nrow(subdf)
      
      p.low  <-  phyper( q = x.samp.succ , m = m.pop.succ , n = n.pop.fail , k = k.samp.size , lower.tail = TRUE )
      p.upp  <-  phyper( q = x.samp.succ , m = m.pop.succ , n = n.pop.fail , k = k.samp.size , lower.tail = FALSE )
      
      rdf  <-  data.frame( classification = classx , group = groupx ,
                           sample.success = x.samp.succ , sample.size = k.samp.size , pop.success = m.pop.succ , pop.size = m.pop.succ + n.pop.fail ,
                           perc.success.sample = x.samp.succ / k.samp.size * 100 ,
                           perc.success.pop = m.pop.succ / (m.pop.succ + n.pop.fail) * 100 ,
                           p.overrepresented = p.upp ,
                           p.underrepresented = p.low )
      resdf  <-  rbind( resdf , rdf )
    } # classx
  } # groupx
  
  return(resdf)
  
} # hypergeometric_test_ORF_type_enrichment



################################################################################
################################################################################
############## ************    Make file for analysis *********** ##############

load("orf.RData") # This is the file generated on Window. Code to generate this file is in below part of script in commented lines.

# ************    Unt Vs TGFb *********** #
classAccum <- load_accumulated_RPF_count_per_ORF(project="161021", classAccum=TRUE)
gtf <- classAccum
fname  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\ codes/out/limma-voom.mrna/161021.rdna_rn18s/NMuMG/blancgrp_161021_Riboseq.unt48--vs--tgfbCX5461.diff-all.txt.gz"
groups  <-  read.table( file = fname , header = TRUE , row.names = 1 , sep="\t" , stringsAsFactors = FALSE , quote = "" )
groups$gene_name  <-  rownames(groups)
groups$group  <-  ifelse( groups$log2FC > 0 , "up" , "down" )


#colnames(gtf) <- gsub("gene_name", "gene",colnames(gtf))
fname  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\ codes/out/limma-voom.mrna/161021.rdna_rn18s/NMuMG/blancgrp_161021_Riboseq.unt48--vs--tgfb48.all.txt.gz"
universe  <-  read.table( file = fname , header = TRUE , row.names = 1 , sep="\t" , stringsAsFactors = FALSE , quote = "" )
universe <- rownames(universe)



################## ************     Run analysis *********** ###################
result.ORFs_unt_vs_tfgb <- test_ORF_classification_read_count_between_conditions(project="161021", rpf=classAccum, groups=groups, species="m39")
result.ORFs_hypergeometric_test_ORF_unt_vs_tfgb <- hypergeometric_test_ORF_type_enrichment(orfs = orfs, universe = universe, gtf = gtf, verbose = FALSE, groups = groups)



# ************    TGFb Vs CX5461 *********** #

fname  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter codes/out/limma-voom.mrna/161021.rdna_rn18s/NMuMG/blancgrp_161021_Riboseq.tgfb48--vs--tgfbCX5461.diff-all.txt.gz"
groups  <-  read.table( file = fname , header = TRUE , row.names = 1 , sep="\t" , stringsAsFactors = FALSE , quote = "" )
groups$gene_name  <-  rownames(groups)
groups$group  <-  ifelse( groups$log2FC > 0 , "up" , "down" )
fname  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter codes/out/limma-voom.mrna/161021.rdna_rn18s/NMuMG/blancgrp_161021_Riboseq.tgfb48--vs--tgfbCX5461.all.txt.gz"
universe  <-  read.table( file = fname , header = TRUE , row.names = 1 , sep="\t" , stringsAsFactors = FALSE , quote = "" )
universe <- rownames(universe)


################## ************     Run analysis *********** ###################
result.ORFs_tgfb_vs_cx5461 <- test_ORF_classification_read_count_between_conditions(project="161021", rpf=classAccum, groups=groups, species="m39")
result.ORFs_hypergeometric_test_ORF_tgfb_vs_cx5461 <- hypergeometric_test_ORF_type_enrichment(orfs = orfs, universe = universe, gtf = gtf, verbose = FALSE, groups = groups)


################################################################################
################################################################################

#test_ORF_classification_read_count_between_conditions_per_gene(project="161021", rpf=classAccum, species="m39")
# Rather than runninh test_ORF_classification_read_count_between_conditions_per_gene function I run code manually
# Because function make so many output files for summary and dump all lines with NaN/NA

project = "161021" ; verbose=FALSE;  species = "m39"

project  <-  process_project_input( project = project )

### load
if (verbose) { verb("\tload.\n") }

### rpf
if (is.null(rpf)) {
  rpf  <-  load_accumulated_RPF_count_per_ORF( project = project , classAccum = TRUE , verbose = verbose )
} # rpf
stopifnot(is.data.frame(rpf))

if (!("condition" %in% colnames(rpf))) {
  if (verbose) { verb("\t\tcond.\n") }
  rpf  <-  merge( rpf , unique(project[,c("runid","condition"),drop=FALSE]) , by = "runid" )
}

### gtf
if (is.null(gtf)) {
  if (verbose) { verb("\t\tload gtf.\n") }
  gtf  <-  load_transcriptome_gtf( species = species )
  gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
} # proteinCodingTranscripts

### prepare gtf
if (verbose) { verb("\t\tprepare gtf.\n") }
gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_name" ) # Here I am usig "gene_name" rather than "gene" in GTF file
gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_id" )
gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_biotype" )
gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_biotype" )
colnames(gtf) <- gsub("gene_name", "gene",colnames(gtf))

#### merge gene name
if (verbose) { verb("\tmerge gene name.\n") }

rpf  <-  merge( rpf , unique(gtf[,c("gene","transcript_id"),drop=FALSE]) , all.x = TRUE )

conpairs  <-  combn( unique(project$condition) , 2)

uniq.classes  <-  unique(rpf$classification)
uniq.classes  <-  uniq.classes[!is.na(uniq.classes)]
show(uniq.classes)

pctids  <-  unique(gtf$transcript_id[gtf$gene_biotype == "protein_coding"  &  gtf$transcript_biotype == "protein_coding"])

resdf  <-  data.frame()

### summarize by gene, class, protein coding
pcdf  <-  rpf[ rpf$transcript_id %in% pctids  ,,drop=FALSE]
pcdf  <-  pcdf  %>%  group_by(gene,classification,runid,condition)  %>%  summarise(total.count = sum(count))
pcdf  <-  as.data.frame(pcdf)

anndf  <-  pcdf[ pcdf$classification == "annotated ORF"  &  !is.na(pcdf$classification) ,,drop=FALSE]
anndf$classification  <-  NULL
colnames(anndf)[colnames(anndf) == "total.count"]  <-  "annotated.ORF.count"
pcdf  <-  merge( pcdf , anndf , by = c("gene","runid","condition") , all = TRUE )
pcdf$log.odds  <-  log2(pcdf$total.count) - log2(pcdf$annotated.ORF.count)
logi.zn  <-  pcdf$total.count == 0  |  is.na(pcdf$total.count)
logi.zd  <-  pcdf$annotated.ORF.count == 0  |  is.na(pcdf$annotated.ORF.count)
pcdf$log.odds[logi.zn & !logi.zd]  <-  log2(pcdf$total.count[logi.zn & !logi.zd] + 1e-6) - log2(pcdf$annotated.ORF.count[logi.zn & !logi.zd])
pcdf$log.odds[!logi.zn & logi.zd]  <-  log2(pcdf$total.count[!logi.zn & logi.zd]) - log2(pcdf$annotated.ORF.count[!logi.zn & logi.zd] + 1e-6 )


for (classx  in  uniq.classes) {
  if (classx == "annotated ORF") { next }
  if (verbose) { verb("\t\t\t%s\n", classx) }
  
  for (cx  in  1:ncol(conpairs)) {
    cond1  <-  conpairs[1,cx]
    cond2  <-  conpairs[2,cx]
    if (verbose) { verb("\t\t\t\t%s  %s\n", cond1, cond2) }
    
    subdf  <-  pcdf[ pcdf$classification == classx  &  !is.na(pcdf$classification)  &  pcdf$condition %in% c(cond1,cond2) ,,drop=FALSE]
    if (nrow(subdf) < 2  ||  length(unique(subdf$condition)) < 2) { next }
    
    
    ### matrix
    if (verbose) { verb("\t\t\t\tmatrix.\n") }
    
    gmat  <-  reshape2::acast( data = subdf , formula = gene ~ runid , value.var = "log.odds" )
    gmat  <-  gmat[complete.cases(gmat),,drop=FALSE]
    if (nrow(gmat) == 0) { next }
    
    ridc1  <-  project$runid[project$condition == cond1]
    logi.c1  <-  colnames(gmat) %in% ridc1
    tres  <-  apply( gmat , 1 , FUN = function(r) t.test( x = r[logi.c1] , y = r[!logi.c1] , alternative = "two.sided" )$p.value )
    tres  <-  data.frame( gene = rownames(gmat) , p.value = tres )
    
    tcdf  <-  reshape2::dcast( data = subdf , formula = gene ~ runid , value.var = "total.count" , fill = 0 )
    aodf  <-  reshape2::dcast( data = subdf , formula = gene ~ runid , value.var = "annotated.ORF.count" , fill = 0 )
    lodf  <-  reshape2::dcast( data = subdf , formula = gene ~ runid , value.var = "log.odds" , fill = 0 )
    colnames(lodf)  <-  paste( colnames(lodf) , "log.odds" , sep="." )
    
    
    ### summarize
    if (verbose) { verb("\t\t\t\tsummarize.\n") }
    
    mdf  <-  merge( tcdf , aodf , by = "gene" , suffixes = c(sprintf(".%s",classx) , ".CDS" ) , all = TRUE )
    mdf  <-  merge( mdf , lodf , by.x = "gene" , by.y = "gene.log.odds" , all = TRUE )
    mdf  <-  merge( mdf , tres , by = "gene" , all = TRUE )
    
    
    ### mean, sd
    if (verbose) { verb("\t\t\t\tmean, sd.\n") }
    
    sdf  <-  subdf  %>%  group_by(gene,classification,condition)  %>%  summarise(mean.log.odds = mean(log.odds) , sd.log.odds = sd(log.odds) )
    sdf  <-  as.data.frame(sdf)
    modf  <-  reshape2::dcast( data = sdf , formula = gene ~ condition , value.var = "mean.log.odds")
    modf$log2FC  <-  modf[[cond2]] - modf[[cond1]]
    sodf  <-  reshape2::dcast( data = sdf , formula = gene ~ condition , value.var = "sd.log.odds")
    
    smdf  <-  merge( modf , sodf , by = "gene" , suffixes = c(".mean",".sd") , all = TRUE )
    mdf  <-  merge( mdf , smdf , by = "gene" , all = TRUE )
    
    
    #if (!is.null(outfbase)) {
    ### write
    #  if (verbose) { verb("\t\t\t\twrite.\n") }
    
    #  fname  <-  sprintf("%s.%s.%s--%s.summary.txt", outfbase , classx , cond1 , cond2 )
    #  write.table( mdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
    # } # outfbase
    
    rdf  <-  mdf[,c("gene","log2FC","p.value"),drop=FALSE]
    rdf$classification  <-  classx
    rdf$condition1  <-  cond1
    rdf$condition2  <-  cond2
    resdf  <-  rbind( resdf , rdf )
  } # cx
} # classx
resdf <- na.omit(resdf)
resdf_pval <- resdf[resdf$p.value <= 0.05,]
write.csv(resdf_pval, file = "ORF_classification_read_count_between_conditions_per_gene.csv", row.names = FALSE)


################################################################################
################################################################################



################################################################################
################################################################################
#This part is done on Windiws system
# Just load workspace of gtf, orfs :: load("orf.RData") 

# suppressMessages(suppressWarnings(source("functions_rsem.R")))
# 
# read.transcriptome.feature.fasta  <-  function( species = NULL , file = NULL , feat )  {
#   
#   library(Biostrings)
#   
#   ### get base transcriptome file name
#   if (is.null(file)) {
#     file  <-  sprintf("%s/%s.gz", matt_get_info( feat = "transcriptomeGTFdir" , species = species ) , matt_get_info( species = species , feat = "transcriptomeGTF"))
#   }
#   
#   basef  <-  gsub( "\\.gtf\\.gz$", "" , file)
#   
#   #	if (feat == "transcript") {
#   #		featf  <-  paste( basef , "fa.gz" , sep=".")
#   #	} else {
#   featf  <-  paste( basef , feat , "fa.gz" , sep=".")
#   #	}
#   
#   ### read fasta
#   infa  <-  readDNAStringSet( featf )
#   
#   return(infa)
# } # read.transcriptome.feature.fasta
# 
# 
# 
# 
# ################################################################################
# ################################################################################
# ################################################################################
# 
# 
# 
# Mus_musculus.GRCm39.gtf <- NULL
# 
# rsem.transcriptome.fasta <- readDNAStringSet("Mus_musculus.GRCm39.104.rdna_rn18s.gtf.rsem.transcriptome.transcripts.fa")
# names(rsem.transcriptome.fasta) <- 
#   names(rsem.transcriptome.fasta) %>%
#   strsplit(., split="::",fixed=TRUE) %>%
#   sapply(., '[', 1) 
# 
# tran <- rsem.transcriptome.fasta
# gtf <- Mus_musculus.GRCm39.gtf
# #outfbase <- "orfs.mm39.Mus_musculus.GRCm39.104.gtf"
# 
# ################################################################################
# ################################################################################
# 
# find_all_pattern_indices_in_string_set  <-  function( seq , pattern , name = "transcript_id" ) {
#   
#   ### get hits
#   mh1  <-  vmatchPattern(pattern = pattern , subject = seq )
#   
#   ### long format
#   hdf  <-  reshape2::melt(mh1) # melt works on this class, miraculously.
#   colnames(hdf)  <-  c("ix","NA","start","end","width")
#   hdf  <-  hdf[, c("ix","start","end") ,drop=FALSE]
#   
#   # name correctly
#   tnames  <-  data.frame( ix = 1:length(names(mh1)) , transcript_id = names(mh1) )
#   colnames(tnames)[2]  <-  name
#   
#   hdf  <-  merge( hdf , tnames , by = "ix" , all.x = TRUE )
#   hdf$pattern  <-  pattern
#   
#   return(hdf)
#   
# } # find_all_pattern_indices_in_string_set
# 
# 
# ################################################################################
# ################################################################################
# # This part I did in Windows
# 
# 
# if(is.null(Mus_musculus.GRCm39.gtf)) {suppressMessages(suppressWarnings(source("make_GRCm39_104_transcript_gtf.R")))}
# 
# if (is.null(tran)) { tran  <-  read.transcriptome.feature.fasta( species = species , feat = "transcripts" )  } 
# stopifnot(class(tran) == "DNAStringSet")  # tran
# 
# ### find all start/stop codons in sequence
# 
# # start
# startdf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "ATG" )
# 
# # stop
# tagdf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "TAG" )
# taadf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "TAA" )
# tgadf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "TGA" )
# stopdf  <-  rbind( tagdf , taadf , tgadf )
# 
# foundf  <-  rbind( startdf , stopdf )
# foundf$id  <-  paste( foundf$transcript_id , foundf$start , foundf$end )
# foundf$codon.type  <-  ifelse( foundf$pattern == "ATG" , "start" , "stop" )
# foundf$mod3  <-  foundf$start %% 3
# 
# ### get all possible start/stop pairs.
# 
# 
# # get start rank and stop rank per entry.  This is very efficient for ensuring that a putative ORF does not contain another STOP inside of it.
# # Note that this MUST be first sorted by mod3 to ensure that the start/stop rank pairs are within the same coding frame.
# foundf  <-  foundf[ order(foundf$mod3 , foundf$transcript_id , foundf$start , foundf$end ) , ,drop=FALSE]
# foundf$rank.stop  <-  cumsum( foundf$codon.type == "stop" )
# 
# # This allows start and stop codons to be matched such that each start codon is ONLY matched with the closest downstream stop codon within the same coding frame (sorted by mod3, above).
# foundf$rank.stop[foundf$codon.type == "start"]  <-  foundf$rank.stop[foundf$codon.type == "start"] + 1 
# 
# # separate into start and stop
# f1df  <-  foundf[ foundf$codon.type == "start" , ,drop=FALSE]
# f2df  <-  foundf[ foundf$codon.type == "stop" , ,drop=FALSE]
# 
# mdf  <-  merge( f1df , f2df , by = c("transcript_id","mod3","rank.stop") , suffixes = c(".start",".stop") )
# mdf$ORF.length  <-  (mdf$start.stop - mdf$start.start) / 3  # in number of codons, includes start but not stop codon
# 
# stopifnot(all(mdf$start.start < mdf$start.stop))
# 
# # must be non-trivial
# #mdf  <-  mdf[ mdf$ORF.length > 1 , ,drop=FALSE]
# ### get known start and stop codons 
# 
# 
# startc  <-  gtf[ gtf$feature == "start_codon" , c("seqname","start","end","feature","frame") , drop=FALSE ]
# stopc  <-  gtf[ gtf$feature == "stop_codon" , c("seqname","start","end","feature","frame") , drop=FALSE ]
# 
# startc$id  <-  paste( startc$seqname , startc$start , startc$end )
# stopc$id  <-  paste( stopc$seqname , stopc$start , stopc$end )
# 
# stopifnot(length(unique(startc$id)) == nrow(startc))
# stopifnot(length(unique(stopc$id)) == nrow(stopc))
# 
# ### get transcript biotype
# 
# 
# tbiot  <-  get_gtf_attribute_field( gtf = gtf , field = "transcript_biotype" )
# gbiot  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_biotype" )
# tdf  <-  unique(data.frame( transcript_id = gtf$transcript_id , transcript_biotype = tbiot , gene_biotype = gbiot ))
# mdf  <-  merge( mdf , tdf , by = "transcript_id" , all.x = TRUE )
# 
# ### identify annotated ORFs
# 
# mdf$isCanonicalORF  <-  (mdf$id.start %in% startc$id  &  mdf$id.stop %in% stopc$id)
# 
# #### classify by relationship to canonical
# 
# colnames(startc)  <-  paste( colnames(startc) , "canonicalStart" , sep  = "." )
# colnames(stopc)  <-  paste( colnames(stopc) , "canonicalStop" , sep  = "." )
# 
# mdf  <-  merge( mdf , startc , by.x = "transcript_id" , by.y = "seqname.canonicalStart" , all.x = TRUE )
# mdf  <-  merge( mdf , stopc , by.x = "transcript_id" , by.y = "seqname.canonicalStop" , all.x = TRUE )
# 
# # classify
# mdf$mod3.canonical  <-  mdf$start.canonicalStart %% 3
# logi.sameFrameAsCanonical  <-  mdf$mod3  ==  mdf$mod3.canonical
# 
# ### Classification groups as defined by Ingolia "Ribosome profiling: new views of translation, from single codons to genome scale" 2014, Figure 3.
# 
# 
# logi.uORF  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  <  mdf$start.canonicalStart)
# logi.overlapping.uORF  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  >  mdf$start.canonicalStart)  &  (mdf$start.stop < mdf$start.canonicalStop)  &  !logi.sameFrameAsCanonical
# logi.alt.ORF  <-  (mdf$start.start  >  mdf$start.canonicalStart)  &  (mdf$start.stop  <  mdf$start.canonicalStop)  &  !logi.sameFrameAsCanonical
# logi.truncation  <-  (mdf$start.start  >  mdf$start.canonicalStart)  &  (mdf$start.stop  ==  mdf$start.canonicalStop)  &  logi.sameFrameAsCanonical
# logi.extension  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  ==  mdf$start.canonicalStop)  &  logi.sameFrameAsCanonical
# 
# # i've generalized a few more
# #prot.cod.biotypes  <-  c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_M_gene","IG_V_gene","IG_Z_gene","nonsense_mediated_decay","nontranslating_CDS","non_stop_decay","polymorphic_pseudogene","protein_coding","TR_C_gene","TR_D_gene","TR_gene","TR_J_gene","TR_V_gene")
# 
# logi.overlapping.uORF.all  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  >  mdf$start.canonicalStart)  &  (mdf$start.stop < mdf$start.canonicalStop)
# logi.dORF  <-  (mdf$start.start  >  mdf$start.canonicalStop)
# logi.overlapping.dORF  <-  (mdf$start.start > mdf$start.canonicalStart)  &  (mdf$start.start < mdf$start.canonicalStop)  &  (mdf$start.stop > mdf$start.canonicalStop) 
# logi.nc  <-  !(mdf$transcript_biotype %in% c("protein_coding"))
# logi.internal  <-  (mdf$start.start  >  mdf$start.canonicalStart)  &  (mdf$start.stop  <  mdf$start.canonicalStop)  &  logi.sameFrameAsCanonical
# 
# # Modification by Nitish on December 20, 2022. Add Rendleman et al. classification
# # Add start-stop uORFs/dORFs based on https://www.biorxiv.org/content/10.1101/2021.07.26.453809v2.full 
# # I will use only uORFs/dORFs for start-stop ORF analysis
# logi.ss.uORF  <-  (logi.uORF == TRUE  &  mdf$ORF.length == 1)
# logi.ss.dORF  <-  (logi.dORF == TRUE  &  mdf$ORF.length == 1)
# 
# ## apply
# ## 
# 
# 
# mdf$classification  <-  NA
# mdf$classification[logi.uORF]  <-  "uORF"
# mdf$classification[logi.overlapping.uORF]  <-  "overlapping uORF"
# mdf$classification[logi.alt.ORF]  <-  "alternative RF"
# mdf$classification[logi.truncation]  <-  "truncation"
# mdf$classification[logi.extension]  <-  "extension"
# mdf$classification[logi.dORF]  <-  "dORF"
# mdf$classification[logi.overlapping.dORF]  <-  "overlapping dORF"
# mdf$classification[logi.overlapping.uORF.all]  <-  "overlapping uORF any frame"
# mdf$classification[logi.nc]  <-  "noncoding ORF"
# mdf$classification[logi.internal]  <-  "internal ORF"
# mdf$classification[mdf$isCanonicalORF]  <-  "annotated ORF"
# # modification by Nitish
# mdf$classification[logi.ss.uORF] <- "start stop uORF"
# mdf$classification[logi.ss.dORF] <- "start stop dORF"
# 
# subg  <-  unique(gtf[,c("gene_name","transcript_id","gene_id"),drop=FALSE])
# mdf  <-  merge( mdf , subg , by = "transcript_id" , all.x = TRUE )
# 
# if(!is.null(outfbase)) {
#   fname  <-  sprintf("%s.orfs.txt" , outfbase )
#   write.table( mdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
#   system(sprintf("gzip  -f  %s" , fname ))
#   # bed
#   orbed  <-  data.frame( chr = mdf$transcript_id , start = mdf$start.start - 1 , end = mdf$end.stop , classification = mdf$classification ) # , paste( mdf$id.start , mdf$id.stop ) , mdf$mod3 , strand = "+" )
#   #orbed  <-  cbind( orbed , mdf[,c("pattern.start","pattern.stop","transcript_biotype","gene_biotype","isCanonicalORF","start.canonicalStart","end.canonicalStop","mod3.canonical","classification"),drop=FALSE] )
#   orbed <- orbed[order(orbed[,1], orbed[,2], orbed[,3]),] # substiture of sort at line 248 in accumulate_rpf_counts_per_putative_ORFs.sh
#   # sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o $TMPWF.orfs.bed  $TMPWF.orfs.bed
#   fname  <-  sprintf("%s.orfs.bed" , outfbase )
#   write.table( orbed , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
#   system(sprintf("gzip  -k -f  %s" , fname ))
#   
# } # outfbase
# 
# ################################################################################
# ################################################################################
# 
# datalist=list(c(gtf=gtf, orfs=mdf, bed=Mus.musculus.GRCm39.BED))
# saveRDS(datalist, "orfs.RDS")
# save(list = c("gtf", "mdf", "Mus.musculus.GRCm39.BED"), file = "orf.RData")
# gtf=gtf; orfs=mdf; GTFbed=Mus.musculus.GRCm39.BED
# save(list = c("gtf", "orfs", "GTFbed"), file = "orf.RData")
# load("orf.RData")



