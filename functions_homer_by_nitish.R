################################################################################
################################################################################
# R functions for the HOMER analysis
# I am trying to make single script for all functions which required for the mofif analysis

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



################################################################################
################################################################################

# verb
verb <- function(...) cat(sprintf(...), sep='', file=stdout())

################################################################################
################################################################################

################################################################################
################################################################################
# R functions for the HOMER analysis
# I am trying to make single script for all functions which required for the mofif analysis

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

de_analysis_type_to_directory_and_tag  <-  function( type , project ) { 
  
  if (type == "mrna") {
    
    #bdir  <-  "out/limma-voom.mrna"
    bdir  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\ codes/out/limma-voom.mrna"
    btag  <-  "limma-voom.mrna"
  } else if (type == "rrna-var") {
    
    #bdir  <-  "out/limma-voom.rrna-var"
    bdir  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/jupyter\ codes/out/limma-voom.mrna"
    btag  <-  "limma-voom.rrna-var"
  } else {
    verb("\n\n\nERROR!  unrecognized type=[%s]!!!\n", type)
    stop()
  } # type
  
  res  <-  list(dir = bdir , tag = btag )
  return(res)
  
} # de_analysis_type_to_directory_and_tag

################################################################################
################################################################################

# load_project_DE_test_tables
# input:
#   file_names: vector
#   project:
#   individual:
# output:
#   data.frame:
#     gene_name	condition1	condition2	expression1	expression2	log2FC	FDR	project	individual
# usage:
# myData.deTest <- load_project_DE_test_tables(c('/home/jupyter/hkim/stjude-nmumg_riboprof/data/2016-11-01--14.45.30.new.mrna.de/emt.161021/emt.161021.limma-voom.mrna.limma.gene.unt48--vs--tgfb48.all.txt','/home/jupyter/hkim/stjude-nmumg_riboprof/data/2016-11-01--14.45.30.new.mrna.de/emt.161021/emt.161021.limma-voom.mrna.limma.gene.tgfb48--vs--tgfbCX5461100nm.all.txt'), project='161021', individual='NMuMG')
#
load_project_DE_test_tables <- function(file_names, project, individual) {
  
  df_de <- NULL
  for (filename in file_names) {
    df1 <- read.table(file=filename,
                      header=TRUE, sep="\t", row.names=1,
                      quote="", comment.char="#", stringsAsFactors=F)
    # unt48	tgfb48	log2FC	FDR
    cols <- colnames(df1)
    condition1 <- cols[1]
    condition2 <- cols[2]
    # convert to simpler condition names. 
    # Modified by Nitish
    ## Modified by Nitish:: my RNAseq result data have unt, tgfb, and tgfbCX. So replace condition1=unt48, and tgfb48, tgfbCX5461 for function "build_gene_set". Riboseq DE data doesn't have any issue.
    ## build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "tgfb48.tgfbCX5461.DEtranscription", verbose = T)
    condition1 <- ifelse(condition1=="tgfb", "tgfb48", ifelse(condition1=="unt", "unt48", ifelse(condition1=="tgfb", "tgfb48", condition1)))
    condition2 <- ifelse(condition2=="tgfb", "tgfb48", ifelse(condition2=="tgfbCX", "tgfbCX5461", ifelse(condition2=="tgfbCX5461100nm", "tgfbCX5461", condition2)))
    
    #condition2 <- gsub('tgfbCX|tgfbCX5461100nm','tgfbCX5461',condition2)
    colnames(df1) <- c('expression1','expression2','log2FC','FDR','p.value')
    df1$project <- project
    df1$individual <- individual
    df0 <- data.frame(gene_name=rownames(df1), condition1=condition1, condition2=condition2)
    df1 <- cbind(df0,df1)
    rownames(df1) <- NULL
    df_de <- rbind(df_de, df1)
  }
  
  return(df_de)
  
}


################################################################################
################################################################################

# build_gene_set
# input:
#   gene_set: 
# output:
#   gene_set$sets: list of list
#     gene_set$sets[['unt48.tgfb48']]$up
#     gene_set$sets[['unt48.tgfb48']]$down
#     gene_set$sets[['unt48.tgfb48']]$notSig
#   gene_set$factors: list of factor: up, down, notSig
#   gene_set$tables: list of data.frame
#     gene_name	condition1	condition2	expression1	expression2	log2FC	FDR	project	individual
#     
# usage:
# gene_set <- build_gene_set(de = myData.deTest, project = "161021", conditions = c("unt48", "tgfb48"), gname="unt48.tgfb48", verbose=T)
# gene_set <- build_gene_set(de = myData.deTest, project = "161021", conditions = c("unt48", "tgfb48", "tgfbCX5461"), gene_set=gene_set, gname="reversible.CX", verbose=T)
build_gene_set  <-  function( de , project , conditions , reversible = length(conditions) > 2 , include.notSig = !reversible , FDR.thresh = 0.05 , element = "gene_name" , gene_set=NULL, gname='default', verbose = FALSE ) {
  
  if (is.null(gene_set)) {
    gene_set <- list()
    gene_set$sets <- list()
    gene_set$factors <- list()
    gene_set$tables <- list()
  } 
  
  if (is.data.frame(project)) { project  <-  unique(project$project) }
  
  gset  <-  list()
  
  if (!reversible) {
    if (verbose) { verb("\tpairwise.\n") }
    
    stopifnot(length(conditions) == 2)
    
    ### DE
    gset$up  <-  de[[element]][ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC > 0  &  de$project %in% project ]
    gset$down  <-  de[[element]][ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC < 0  &  de$project %in% project ]
    
    ### factor and notsig
    if (include.notSig) {
      gset$notSig  <-  de[[element]][ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR > FDR.thresh  &  de$project %in% project ]
      gfac  <-  factor(names(gset) , levels = c("notSig","down","up"))
    } else {
      gfac  <-  factor(names(gset) , levels = c("down","up"))
    } # n.s.
    
    ### table
    gtab  <-  de[ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]   &  de$project %in% project ,,drop=FALSE]
    
  } else {
    if (verbose) { verb("\treversible.\n") }
    
    stopifnot(length(conditions) == 3)
    stopifnot(!include.notSig)
    
    ### rev
    up1  <-  de[[element]][(de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC > 0  &  de$project %in% project)]
    down1  <-  de[[element]][(de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC < 0  &  de$project %in% project)]
    
    up2  <-  de[[element]][(de$condition1 == conditions[2]  &  de$condition2 == conditions[3]  &  de$FDR <= FDR.thresh  &  de$log2FC > 0  &  de$project %in% project)]
    down2  <-  de[[element]][(de$condition1 == conditions[2]  &  de$condition2 == conditions[3]  &  de$FDR <= FDR.thresh  &  de$log2FC < 0  &  de$project %in% project)]
    
    gset$upDown  <-  intersect( up1 , down2 )
    gset$downUp  <-  intersect( down1 , up2 )
    
    ### factor
    gfac  <-  factor(names(gset) , levels = c("downUp","upDown"))
    
    ### table
    gtab  <-  de[ ( (de$condition1 == conditions[1]  &  de$condition2 == conditions[2])  | (de$condition1 == conditions[2]  &  de$condition2 == conditions[3]) )  & de$project %in% project ,,drop=FALSE]
  } # rev
  
  gene_set$sets[[gname]] <- gset
  gene_set$factors[[gname]] <- gfac
  gene_set$tables[[gname]] <- gtab
  return(gene_set)
  
} # build_gene_set



################################################################################
################################################################################

# build_gene_set_after_comparison
# input:
#   type: {'agree','complement','opposite'}
# output:
#   gene_set$sets
#   gene_set$factors
#   gene_set$tables
# usage:
# gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rna, name1="unt48.tgfb48", name2="rna.unt48.tgfb48", type="complement", gname="unt48.tgfb48.ribo.only")
build_gene_set_after_comparison <- function(gene_set, gene_set2, name1, name2, type, gname) {
  
  
  res  <-  compare_gene_set( set1 = gene_set$sets[[name1]] , set2 = gene_set2$sets[[name2]] , fac = gene_set$factors[[name1]] , table1 = gene_set$tables[[name1]] , table2 = gene_set2$tables[[name2]] , type = type )
  
  gene_set$sets[[gname]] <- res$set
  gene_set$factors[[gname]] <- res$factor
  gene_set$tables[[gname]] <- res$table
  
  return(gene_set)
  
}


################################################################################
################################################################################

# compare_gene_set
# input:
# output:
# usage:
# gname  <-  "translationONLY.unt48.tgfb48"
# name1  <-  "riboProf.unt48.tgfb48"
# name2  <-  "polyA.unt48.tgfb48"
# res  <-  compare_gene_set( set1 = gene_set$sets[[name1]] , set2 = gene_set$sets[[name2]] , factor1 = gene_set$factors[[name1]] , table1 = gene_set$tables[[name1]] , table2 = gene_set$tables[[name2]] , type = "complement" )
#
compare_gene_set  <-  function( set1 , set2 , fac , table1 , table2 , type ) {
  
  
  common.names  <-  intersect( names(set1) , names(set2) )
  common.names  <-  setdiff( common.names , "notSig" )
  stopifnot(length(common.names) > 1)
  
  gset  <-  list()
  if (type == "agree") {
    for (namex  in  common.names) {
      gset[[namex]]  <-  intersect( set1[[namex]] , set2[[namex]] )
    } # namex
    
  } else if (type == "opposite") {
    stopifnot(length(common.names) == 2)
    for (namex  in  common.names) {
      name2  <-  setdiff( common.names , namex )
      gset[[namex]]  <-  intersect( set1[[namex]] , set2[[name2]] )
    } # namex
  } else if (type == "complement") {
    for (namex  in  common.names) {
      gset[[namex]]  <-  setdiff( set1[[namex]] , set2[[namex]] )
    } # namex
  } else {
    verb("\n\n\nERROR!  unrecognized type=[%s] for comparing gene sets!!!\n", type)
    stop()
  } # type
  
  ### factor
  gfac  <-  as.character(fac)
  gfac  <-  setdiff( gfac , "notSig" )
  gfac  <-  factor( gfac , levels = setdiff(levels(fac) , "notSig" ) )
  stopifnot(!any(is.na(gfac)))
  
  ### table
  gtab  <-  rbind( table1 , table2 )
  
  res  <-  list( set = gset , factor = gfac , table = gtab )
  
  return(res)
  
} # compare_gene_set



################################################################################
################################################################################

get_gene_biotypes_and_protein_coding  <-  function( gtf ) {
  # returns a list with fields "biotype" and "proteinCoding"
  #       "biotype"  =  character vector of unique biotypes from gtf
  #       "proteinCoding" = character vector of transcript ids that are protein coding
  
  biot  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_biotype" )
  myData.geneBiotype  <-  data.frame( gene_id = gtf$gene_id , biotype = biot )
  myData.geneBiotype  <-  unique(myData.geneBiotype)
  
  myData.proteinCoding  <-  unique(myData.geneBiotype$gene_id[myData.geneBiotype$biotype == "protein_coding"])
  
  return(list(biotype = myData.geneBiotype , proteinCoding = myData.proteinCoding))
  
} # get_gene_biotypes_and_protein_coding



################################################################################
################################################################################

get_canonical_transcript_per_gene_name  <-  function( gtf ) {
  ## DESCRIPTION
  #       This function will convert "gene names/symbol" to a canonical transript id.
  # The problem is that, outrageously, gene name to gene id (ensembl) is one-to-many for some gene names.
  # Thus, while gene id to canonical gene transcript is one-to-one, to get from gene name to canonical gene
  # transcript is nontrivial.
  #
  # INPUT
  #       gtf             = transcriptome gtf (data frame)
  #       canon           = a data frame containing canonical transcript_id per gene_id.
  #                        if NULL, then will be inferred from the gtf.
  #
  # OUTPUT
  #       dataframe with columns:
  #                gene_id , transcript_id , gene , CDS.length , five_prime_utr.length , three_prime_utr.length
  
  
  gtf  <-  prepare_transcriptome_gtf( gtf )
  # 'seqname' 'source' 'feature' 'start' 'end' 'score' 'strand' 'frame' 'attribute' 'transcript_id' 'gene_name' 'gene_id'
  canon  <-  get_canonical_transcripts(gtf)
  utab  <-  unique(gtf[,c("gene_name","gene_id","transcript_id") ,drop=FALSE])
  
  gdf  <-  merge( utab , canon , by = c("gene_id","transcript_id") )
  
  gdf  <-  gdf[ order(gdf$CDS.length,gdf$five_prime_utr.length,gdf$three_prime_utr.length) , ,drop=FALSE]
  gdf  <-  gdf[ !duplicated(gdf$gene_name) , ,drop=FALSE]
  
  return(gdf)
  
} # get_canonical_transcript_per_gene_name



################################################################################
################################################################################
get_canonical_transcripts  <-  function( gtf ) {
  # returns a data frame with 3 columns:  "gene_id", "transcript_id", "cds.length".
  # each gene appears only once.
  
  # order by longest CDS, 5'UTR, 3'UTR (in that order)
  gtf$feat.length  <-  gtf$end - gtf$start + 1
  gtf  <-  gtf[ gtf$feature %in% c("five_prime_utr","three_prime_utr","CDS") , ,drop=FALSE]
  gtf  <-  gtf[ order(gtf$feature, -gtf$feat.length) , ,drop=FALSE]
  
  # canonincal set
  uniqdf  <-  gtf[ !duplicated(gtf$gene_id) , c("gene_id","transcript_id") ,drop=FALSE]
  
  # feat length table
  tfeat  <-  reshape2::dcast( gtf , transcript_id ~ feature , value.var = "feat.length" , fill = 0 )
  loc.o  <-  which(colnames(tfeat)  %in%  c("five_prime_utr","three_prime_utr","CDS"))
  colnames(tfeat)[loc.o]  <-  paste( colnames(tfeat)[loc.o] , "length" , sep=".")
  
  # put feature lengths onto canonical df
  uniqdf  <-  merge( uniqdf , tfeat , by = "transcript_id" , all.x = TRUE )
  
  return(uniqdf)
  
} # get_canonical_transcripts




################################################################################
################################################################################

prepare_transcriptome_gtf  <-  function( gtf ) {
  
  if (!("transcript_id" %in% colnames(gtf))) {
    gtf$transcript_id  <-  gtf$seqname
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

get_all_transciptome_classification_info  <-  function( gtf ) {
  
  utab  <-  unique(gtf[,c("gene_name","gene_id","transcript_id") ,drop=FALSE])
  biot  <-  get_gene_biotypes_and_protein_coding( gtf = gtf )
  geneBiotype  <-  biot$biotype
  # consider only protein coding genes...
  proteinCodingGenes  <-  biot$proteinCoding
  
  canonTrans  <-  get_canonical_transcripts( gtf = gtf )
  protCod.canon  <-  canonTrans[ canonTrans$gene_id %in% proteinCodingGenes , ,drop=FALSE]
  
  canonTransName  <-  get_canonical_transcript_per_gene_name( gtf = gtf )
  
  retval  <-  list( utab = utab,
                    geneBiotype = geneBiotype ,
                    proteinCodingGenes = proteinCodingGenes ,
                    canonicalTranscripts = canonTrans ,
                    canonicalTranscriptsGeneName = canonTransName ,
                    proteinCodingCanonTrans = protCod.canon )
  
  return(retval)
  
} # get_all_transciptome_classification_info

################################################################################
################################################################################

# build_transcriptome_info
# input:
#   species:
# output:
#   list$myData.transcriptome.gtf: data.frame of columns ("seqname", "source", "feature", "start", "end", "frame", "attribute", "transcript_id", "gene_name", "gene_id")
#      seqname: transcript
#   list$myData.transcriptome.info
#      utab: data.frame of columns ("gene_name", "gene_id", "transcript_id")
#      geneBiotype: data.frame of columns ("gene_id", "biotype")
#      proteinCodingGenes: vector of gene_id
#      canonicalTranscripts: data.frame of columns ("transcript_id", "gene_id", "CDS.length", "five_prime_utr.length", "three_prime_utr.length")
#      canonicalTranscriptsGeneName: data.frame of columns ("gene_id", "transcript_id", "gene_name", "CDS.length", "five_prime_utr.length", "three_primer_utr.length")
#      proteinCodingCanonTrans: data.frame of columns ("transcript_id", "gene_id", "CDS.length", "five_prime_utr.length", "three_primer_utr.length")
# comment:
# consider using load(sprintf("%s/rdata/jupyter_common_with_gtf.rdata", dir_jupyter))
# see make_jupyter_common.ipynb
# usage:
# list_out <- build_transcriptome_info(species="m39")
build_transcriptome_info <- function(species = species_notebook) {
  
  myData.species  <-  species
  myData.transcriptome.gtf  <-  load_transcriptome_gtf( species = myData.species )
  myData.transcriptome.gtf  <-  prepare_transcriptome_gtf( gtf = myData.transcriptome.gtf )
  
  # see functions_transcriptome.R
  #gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_length" , num = TRUE )
  #gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "five_prime_utr_length" , num = TRUE )
  #gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "three_prime_utr_length" , num = TRUE )
  
  myData.transcriptome.info  <-  get_all_transciptome_classification_info( gtf = myData.transcriptome.gtf )
  #myData.transcriptome.feats.fa <- load_all_transcriptome_gtf_feature_fastas(species = myData.species, gtf = myData.transcriptome.gtf)
  #myData.genome.gtf  <-  load_genome_gtf( species = myData.species)
  
  list(myData.transcriptome.gtf=myData.transcriptome.gtf, myData.transcriptome.info=myData.transcriptome.info)
  
} # build_transcriptome_info


################################################################################
################################################################################

convert_ensembl_gene_id_to_gene_name_and_aggregate_matrix  <-  function( mat , species = NULL , fun = sum , gtf = NULL , verbose = FALSE ) {
  
  if (is.null(gtf)) {
    if (verbose) { verb("\t\tloading gtf.\n") }
    gtf  <-  load_genome_gtf( species = species )
    
    if (verbose) { verb("\t\tpreparing gtf.\n") }
    gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
  } # gtf
  
  
  #### table
  if (verbose) { verb("\t\tconvert table.\n") }
  
  namedf  <-  unique(gtf[,c("gene_name","gene_id"),drop=FALSE])
  
  if (verbose) { show(head(namedf)) }
  
  
  ### melt
  if (verbose) { verb("\t\tmelt.\n") }
  
  mdf  <-  reshape2::melt(mat)
  colnames(mdf)  <-  c("gene_id" , "col" , "value" )
  
  if (verbose) { show(head(mdf)) }
  
  mdf  <-  merge( mdf , namedf , by = "gene_id" , all.x = TRUE )
  
  if (verbose) { verb("\t\tmerge.d\n") ; show(head(mdf)) }
  
  # begin of modification by H. Kim
  #mdf$gene_name[is.na(mdf$gene_name)]  <-  mdf$gene_id
  f_no_conversion <- is.na(mdf$gene_name)
  if (any(f_no_conversion)) {
    # ENS1,ENS2 --> sym1,sym2
    mdf[f_no_conversion, "gene_name"] <- sapply(as.character(mdf[f_no_conversion, "gene_id"]), function(x) {
      items <- strsplit(x, ',')[[1]]
      idx <- fmatch(items, gtf$gene_id)
      if (all(!is.na(idx))) {
        paste(gtf[idx,"gene_name"], collapse=',')
      } else {
        NA
      }
    })
    f_no_conversion <- is.na(mdf$gene_name)
    mdf[f_no_conversion, "gene_name"]  <-  as.character(mdf[f_no_conversion, "gene_id"])
  }
  # end of modification
  mdf  <-  mdf  %>%  group_by(gene_name,col)  %>%  summarise( total.val = fun(value) )
  mdf  <-  as.data.frame(mdf)
  
  
  if (verbose) { verb("\t\tcast.\n") ; show(head(mdf)) }
  
  amat  <-  acast( data = mdf , formula = gene_name ~ col , value.var = "total.val" )
  
  return(amat)
  
} # convert_ensembl_gene_id_to_gene_name_and_aggregate_matrix


################################################################################
################################################################################

#' Restrict orthologous gene table to given species
#'
#' This is an internal function.  
#'
#' @param	ortho		data frame of ortholgoous gene table, e.g. as returned by function load_orthologous_gene_maps_MouseHuman
#' @param	from, to	character vectors of species names
#'
#' @return			data frame, restricted to species
#'
#' @export
restrict_orthologous_genes_by_species  <-  function( ortho , from = NULL , to = NULL ) {
  
  # legacy mouse nomenclature
  if (from == "m39") { from  <-  "mouse" }
  if (to == "m39") { to  <-  "mouse" }
  
  ortho$species  <-  str_extract( ortho$Common.Organism.Name , "^[^,]+" )
  ortho  <-  subset( ortho , species %in% c(from , to) )
  
  return(ortho)
  
} # restrict_orthologous_genes_by_species


################################################################################
################################################################################

#' Loads the orthologous gene table provided by MGI Jax
#'
#' From: http://www.informatics.jax.org/homology.shtml
#'	http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
#'
#' @param from, to	if specified, the table will be restricted to these species
#'
#' @return		data frame
#'
#' @export
load_orthologous_gene_maps_MouseHuman  <-  function( from = NULL , to = NULL) {
  
  orth.dir  <-  "/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/data/mgi"
  orth.name  <-  "HOM_MouseHumanSequence.rpt.gz"
  
  fname  <-  sprintf("%s/%s", orth.dir , orth.name )
  
  indf  <-  read.table( file = fname , header = TRUE , sep = "\t" , row.names = NULL , stringsAsFactors = FALSE , quote = "" , fill = TRUE )
  
  if (!is.null(from)  ||  !is.null(to))  {
    indf  <-  restrict_orthologous_genes_by_species( ortho = indf , from = from , to = to )
  } # from, to
  
  return(indf)
  
} # load_orthologous_gene_maps_MouseHuman

################################################################################
################################################################################

#' Get orthologous gene table
#'
#' @param  from		[char] species name 
#' @param  to		[char] species name
#' @param  conv		[data frame] orthologous gene table from MGI (e.g. use function load_orthologous_gene_maps_MouseHuman)
#' @param  genes	[character vector] gene name ("symbols") corresponding to the "from" species.
#'
#' @return		[data frame] with columns given by "from","to".  Note that mappings may not be unique, so a gene may appear multiple times on either list.
get_orthologous_gene_table  <-  function( from , to , conv , genes ) {
  
  # legacy mouse nomenclature
  if (from == "m39") { from <- "mouse" }
  if (to == "m39") { to <- "mouse" }
  
  # prepare
  conv$species  <-  str_extract( conv$Common.Organism.Name , "^[^,]*" )
  conv  <-  subset( conv , species %in% c(from , to) )
  
  # homol ids, for genes in "from"
  fdf  <-  subset( conv , species == from  &  Symbol %in% genes )
  hgi  <-  fdf$HomoloGene.ID
  
  # get same homol ids, for "to"
  tdf  <-  subset( conv , species == to  &  HomoloGene.ID %in% hgi )
  
  return(tdf)
  
} # get_orthologous_gene_table


################################################################################
################################################################################

#' Map orthologous gene names between two species
#'
#' @param  from		[char] species name 
#' @param  to		[char] species name
#' @param  conv		[data frame] orthologous gene table from MGI (e.g. use function load_orthologous_gene_maps_MouseHuman)
#' @param  genes	[character vector] gene name ("symbols") corresponding to the "from" species.
#'
#' @return		character vector of gene names homologous to those given in "genes", but not necessary in same order.
get_orthologous_gene_names  <-  function( from , to , conv , genes ) {
  
  tdf  <-  get_orthologous_gene_table( from = from , to = to , conv = conv , genes = genes )
  tog  <-  tdf$Symbol
  tog  <-  unique(tog)
  
  return(tog)
  
} # get_orthologous_gene_names

################################################################################
################################################################################

# add_ortholog_info
# add columns of ortholog information (e.g., HomoloGene.ID, mouse.sym, mouse.eid, human.sym, human.eid)
add_ortholog_info <- function(df, from="mouse", to="human") {
  
  df_from <- orthodf[orthodf$specie==from,]
  idx <- match(rownames(df), df_from[,"Symbol"])
  df$HomoloGene.ID <- df_from[idx,"HomoloGene.ID"]
  df[, sprintf("%s.sym",from)] <- df_frome[idx,"Symbol"]
  df[, sprintf("%s.eid",from)] <- df_from[idx,"EntrezGene.ID"]
  
  df_to <- orthodf[orthodf$specie==to,]
  idx <- match(df$HomoloGene.ID , df_to$HomoloGene.ID)
  df[, sprintf("%s.sym",to)] <- df_to[idx,"Symbol"]
  df[, sprintf("%s.eid",to)] <- df_to[idx,"EntrezGene.ID"]
  
  return(df)
  
}
################################################################################
################################################################################

write_fasta_with_boundary <- function( fname_fasta, list_boundary, featx, outpref, alphabet) {
  
  f_search_codon_pattern_in_cds <- FALSE
  if (grepl("CDS", fname_fasta) && alphabet=="AA") {
    f_search_codon_pattern_in_cds <- TRUE
    alphabet <- "RNA"
    if (any(names(list_boundary) == featx)) {
      # update list_boundary
      s <- list_boundary[[featx]][1]
      e <- list_boundary[[featx]][2]
      if (s < 0) s <- s*3        else s <- (s-1)*3+1
      if (e < 0) e <- (e+1)*3-1  else e <- e*3
      list_boundary[[featx]] <- c(s,e)
    }
  }
  
  switch(alphabet,
         "DNA"={ ss <- readDNAStringSet(fname_fasta) },
         "RNA"={
           #ss <- readRNAStringSet(fname_fasta)
           ss <- readDNAStringSet(fname_fasta)
         },
         "AA"={ ss <- readAAStringSet(fname_fasta) },
         {}
  )
  
  if (f_search_codon_pattern_in_cds) {
    mycode <- GENETIC_CODE
    #Selenocysteine  Sec/U   UGA     https://en.wikipedia.org/wiki/Selenocysteine
    mycode[["TGA"]] <- "U"
    #Pyrrolysine     Pyl/O   UAG     https://en.wikipedia.org/wiki/Pyrrolysine
    mycode[["TAG"]] <- "O"
    
    # trim untranslated tail sequences.
    n_nucleotides_of_tail <- width(ss) %% 3
    ss <- subseq(ss, 1, -1-n_nucleotides_of_tail) 
    # translation
    suppressWarnings( aa  <-  Biostrings::translate( ss, genetic.code=mycode, no.init.codon=FALSE, if.fuzzy.codon="error" ) )
    
    # remove the last stop codon for interproscan.sh.
    # see make_transcriptome_feature_protein_fasta.R
    f <- as.character(subseq(aa, -1,-1)) %in% c("*","U","O")
    ss[f] <- subseq(ss[f], 1, -4)
    aa[f] <- subseq(aa[f], 1, -2)
    
    # remove polymorphic pseudogenes that contain stop codons inside aa sequences.
    # see make_transcriptome_feature_protein_fasta.R
    # ENSMUST00000126664.7
    # ENSMUST00000134773.3
    # ENSMUST00000135748.1
    f <- grepl("\\*", as.character(aa))
    ss <- ss[!f]
  }
  
  if (any(names(list_boundary) == featx)) {
    n <- width(ss)
    s <- list_boundary[[featx]][1]
    e <- list_boundary[[featx]][2]
    if (s < 0) s <- pmax(-n, s) else s <- pmin(s, n)
    if (e < 0) e <- pmax(-n, e) else e <- pmin(e, n)
    
    f.na <- (sign(s*e) < 0) & (abs(s)+abs(e) > n)
    ss[f.na] <- subseq(ss[f.na], 1, width=0)
    ss[!f.na] <- subseq(ss[!f.na], s[!f.na], e[!f.na])
  }
  
  # write a temporary fasta file
  fname_out <- sprintf("%s.fa", outpref)
  writeXStringSet(ss, filepath=fname_out, compress=FALSE, format="fasta")
  
  fname_out
  
}

################################################################################
################################################################################

# INTERNAL
get_acceptable_iupac_letters  <-  function() {
  let  <-  c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N")
  return(let)
} # get_acceptable_iupac_letters

################################################################################
################################################################################


find_motif_in_fasta  <-  function( motif , name , fasta , outpref , mismatch = 0 , discard.target = TRUE , verbose = TRUE ) {
  
  ### make Homer motif file
  if (verbose) verb("\tmake Homer motif file.\n")
  
  # begin of addition by H. Kim
  if (grepl('\\.', motif)) {
    mo.fname <- motif
    inmotif  <-  load_homer_motif( motif = mo.fname )
  } else {
    # end of addition
    # check
    iupac  <-  get_acceptable_iupac_letters()
    stopifnot(all(sapply( str_split(motif , "" ) , FUN=function(x)  x %in% iupac )))
    
    mo.fname  <-  sprintf("%s.motif" , outpref )
    
    cmd  <-  sprintf("$HOMER/seq2profile.pl  %s  %d  %s > %s" , motif , mismatch , name , mo.fname )
    system(cmd)
    
    # check
    inmotif  <-  load_homer_motif( motif = mo.fname )
    stopifnot(nrow(inmotif$matrix) == nchar(motif) )
    # begin of addition by H. Kim
  }
  # end of addition
  
  
  local.fa.fname  <-  sprintf("%s.target.fa" , outpref )
  
  if (class(fasta) == "DNAStringSet") {
    writeXStringSet( fasta , local.fa.fname )
  } else if (is.character(fasta)  &&  length(fasta) == 1) {
    system(sprintf("zcat  -f  %s  >  %s" , fasta , local.fa.fname))
  } else {
    verb("\n\n\nunrecognized fasta!!!!\n")
    show(class(fasta))
    stop()
  } # fasta
  
  
  ### find motifs
  if (verbose) verb("\tfind motifs in fasta.\n")
  
  ffound  <-  sprintf("%s.found" , outpref)
  #find.dir  <-  sprintf("%s.find" , outpref )
  #cmd  <-  sprintf("$HOMER/findMotifs.pl  %s  fasta  %s  -fasta %s  -find %s  >  %s.found.bed" , fasta , find.dir , fasta , mo.fname , outpref )
  cmd  <-  sprintf("$HOMER/homer2  find  -i %s  -m %s  -o %s  -offset 1 " , local.fa.fname , mo.fname , ffound )
  # with offset = 1, then if the motif starts at the first position in the FASTA and is on the forward strand, then in the results file it will have position = 1
  
  system(cmd)
  
  
  if (discard.target) {
    file.remove(local.fa.fname)
  } # discard.target
  
} # find_motif_in_fasta


################################################################################
################################################################################
# Internal
# This loads the file created by "seq2profile.pl" function
load_homer_motif  <-  function( motif ) {
  
  header  <-  readLines( con = motif , n = 1 )
  
  ### matrix
  inmat  <-  read.table( file = motif , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , skip = 1 )
  inmat  <-  as.matrix(inmat)
  
  inmotif  <-  list( header = header , matrix = inmat )
  
  return(inmotif)
  
} # load_homer_motif

################################################################################
################################################################################

# This loads the file resulting from "findMotifs.pl" with the "-find" option.
#' Load motif hits from searching for a motif within target sequence.
#'
#' After searching for motifs in a target fasta, e.g. with "find_motif_in_fasta",
#' this function will load the results.
#' Note that the "sequence" given by a HOMER hit always corresponds to the forward strand, even if the motif was found on the reverse strand.
#'
#' @param  hits		filename giving the results from HOMER motif search, e.g. from function "find_motif_in_fasta"
#' @param  bed		if TRUE, will convert the results to BED format (no loss of information).
#'
#' @return  if bed = TRUE, then will return a dataframe in BED format.
#'		if bed = FALSE, will return a dataframe in HOMER results format.
load_homer_motif_annotation  <-  function( hits , bed = FALSE ) {
  
  # begin of modification by H. Kim
  #indf  <-  read.table( file = hits , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , fill = TRUE )
  indf  <-  tryCatch( read.table( file = hits , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , fill = TRUE ), error=function(e) NULL )
  if (is.null(indf)) return(indf)
  # end of modification
  
  colnames(indf)  <-  c("chr","offset","sequence","motif","strand","score")
  
  if (bed) {
    indf  <-  homer_motif_annotation_to_bed( annot = indf )
  } # bed
  
  return(indf)
  
} # load_homer_motif_annotation



################################################################################
################################################################################

# Internal
# This loads the file created by "seq2profile.pl" function
load_homer_motif  <-  function( motif ) {
  
  header  <-  readLines( con = motif , n = 1 )
  
  ### matrix
  inmat  <-  read.table( file = motif , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , skip = 1 )
  inmat  <-  as.matrix(inmat)
  
  inmotif  <-  list( header = header , matrix = inmat )
  
  return(inmotif)
  
} # load_homer_motif


################################################################################
################################################################################
# Internal
homer_motif_annotation_to_bed  <-  function( annot ) { 
  
  annot$start  <-  ifelse( annot$strand == "+" , annot$offset - 1 , annot$offset - nchar(annot$sequence) )
  annot$end  <-  annot$start + nchar(annot$sequence)
  
  abed  <-  data.frame( chr = annot$chr , start = annot$start , end = annot$end , name = annot$motif , score = annot$score , strand = annot$strand , sequence = annot$sequence )
  return(abed)
  
} # homer_motif_annotation_to_bed



################################################################################
################################################################################

#' Test enrichment and depletion via hypergeometric distribution
#'
#' @param  universe		character vector of element names for the entire universe
#' @param  classes		Data frame. Master table containing all gene <-> class annotations.
#'					Must have columns "element","group", where element is e.g. gene name,
#'					and where group is e.g. the regulating transcript factor, etc.
#'					The universe may have the same element in multiple different classifications 
#'					(e.g. the same gene regulated by multiple different transcriptiion factors.)
#'					It may also have group = NA. These will not be tested as a class.
#' @param  groups		data frame: must have columns "element" and "group". element is e.g. gene name, group is e.g. "up", "down", etc.
#'					The same element may appear in multiple different groups.
#' @return dataframe of p-value info
#'
#' @export
test_hypergeometric_enrichment  <-  function( universe , classes , groups , verbose = FALSE ) {
  
  ##### check
  if (verbose) { verb("\tchecking.\n") }
  
  stopifnot(is.data.frame(groups))
  stopifnot(is.data.frame(classes))
  stopifnot(all(complete.cases(groups)))
  stopifnot(all(complete.cases(classes)))
  stopifnot(all(c("element","group") %in% colnames(classes)))
  stopifnot(all(c("element","group") %in% colnames(groups)))
  stopifnot(length(setdiff(groups$element , universe)) == 0)
  stopifnot(length(setdiff(classes$element , universe)) == 0)
  
  
  ##### prepare
  if (verbose) { verb("\tpreparing.\n") }
  
  # unique
  classes  <-  unique(classes[,c("group","element"),drop=FALSE])
  groups  <-  unique(groups[,c("group","element"),drop=FALSE])
  
  # prepare
  pop.size  <-  length(unique(universe))
  
  classes$group  <-  factor( as.character(classes$group) , levels = unique(as.character(classes$group)) )
  
  ctab  <-  table(classes$group)
  classdf  <-  data.frame( group = names(ctab) , pop.success = as.numeric(ctab) )
  
  
  ##### check again
  if (verbose) { verb("\tchecking again.\n") }
  
  stopifnot(!any(is.na(groups)))
  stopifnot(!any(is.na(classes)))
  stopifnot(!any(is.na(classdf)))
  
  
  resdf  <-  data.frame()
  for (groupx  in  unique(groups$group)) {
    if (verbose) { verb("\t\t%s\n" , groupx ) }
    
    genes  <-  groups$element[groups$group == groupx]
    
    subclass  <-  classes[ classes$element %in% genes  ,,drop=FALSE]
    stab  <-  table(subclass$group)
    
    scdf  <-  data.frame( group = names(stab) , sample.success = as.numeric(stab) )
    overdf  <-  merge( classdf , scdf , by = "group" , all = TRUE)
    stopifnot(!any(is.na(overdf)))
    
    q  <-  overdf$sample.success
    m  <-  overdf$pop.success
    n  <-  pop.size - overdf$pop.success
    k  <-  length(genes)
    
    # Recall that:
    #	lower.tail = TRUE  -->  P[X <= x]
    #	lower.tail = FALSE -->  P[X  > x]
    
    utail  <-  phyper( q = q , m = m , n = n , k = k , lower.tail = FALSE , log.p = FALSE )
    dens.x  <-  dhyper( x = q , m = m , n = n , k = k , log = FALSE )
    p.over  <-  utail + dens.x
    p.under  <-  1.0 - utail
    
    rdf  <-  data.frame( class = overdf$group , group = groupx ,
                         success.sample = q , sample.size = k , success.population = m , population.size = pop.size ,
                         frac.success.sample = q / k , frac.success.population = m / pop.cd ize ,
                         p.overrepresented = p.over , p.underrepresented = p.under )
    resdf  <-  rbind( resdf , rdf )
  } # groupx
  
  resdf  <-  resdf[order(resdf$group , resdf$p.overrepresented , resdf$success.population , resdf$sample.size) ,,drop=FALSE]
  
  return(resdf)
  
} # test_hypergeometric_enrichment



################################################################################
################################################################################
# build_gene_set_with_rnaseq_riboseq
# input:
#   pattern_rnaseq_gene: (e.g. n-R5s)
#   pattern_riboseq_gene: (e.g. n-R5s)
# output:
#   gene_set: list
#   gene_set$sets: list of list
#     gene_set$sets[['unt48.tgfb48']]$up
#     gene_set$sets[['unt48.tgfb48']]$down
#     gene_set$sets[['unt48.tgfb48']]$notSig
#   gene_set$factors: list of factor: up, down, notSig
#   gene_set$tables: list of data.frame
#     gene_name condition1      condition2      expression1     expression2     log2FC  FDR     project individual
# usage:
# gene_set <- build_gene_set_with_rnaseq_riboseq()
# gene_set <- build_gene_set_with_rnaseq_riboseq(pattern_rnaseq_gene="n-R5s", pattern_riboseq_gene="n-R5s")
#build_gene_set_with_rnaseq_riboseq <- function(project_rnaseq="170224", project_riboseq="161021", rundata_appendix=".rdna_rn18s", individual="NMuMG", level="htseq_gene", f_include_tgfb_vs_cx5461=FALSE, pattern_rnaseq_gene=NULL, pattern_riboseq_gene=NULL) {
build_gene_set_with_rnaseq_riboseq <- function(project_rnaseq="211613", project_riboseq="161021", rundata_appendix=".rdna_rn18s", individual="NMuMG", level="htseq_gene", f_include_tgfb_vs_cx5461=FALSE, pattern_rnaseq_gene=NULL, pattern_riboseq_gene=NULL) {
  
  
  list_out <- de_analysis_type_to_directory_and_tag(type = "mrna", project_rnaseq)
  bdir <- list_out$dir
  btag <- list_out$tag
  
  dir_de <- sprintf("%s/%s%s/%s", bdir, project_rnaseq, rundata_appendix, individual)
  fname1 <- sprintf("%s/blancgrp_%s_RNAseq_total_stranded.unt--vs--tgfb.diff-all.txt.gz", dir_de, project_rnaseq)
  fname2 <- sprintf("%s/blancgrp_%s_RNAseq_total_stranded.tgfb--vs--tgfbCX.diff-all.txt.gz", dir_de, project_rnaseq)
  rnaseq.deTest <- load_project_DE_test_tables(c(fname1, fname2), project = project_rnaseq, individual = individual)
  
  
  if (!is.null(pattern_rnaseq_gene)) {
    f <- grepl(pattern_rnaseq_gene, rnaseq.deTest$gene_name)
    rnaseq.deTest <- rnaseq.deTest[f,]
  }
  
  dir_de <- sprintf("%s/%s%s/%s", bdir, project_riboseq, rundata_appendix, individual)
  fname1 <- sprintf("%s/blancgrp_%s_Riboseq.unt48--vs--tgfb48.all.txt.gz", dir_de, project_riboseq)
  fname2 <- sprintf("%s/blancgrp_%s_Riboseq.tgfb48--vs--tgfbCX5461.all.txt.gz", dir_de, project_riboseq)
  riboseq.deTest <- load_project_DE_test_tables(c(fname1, fname2), project = project_riboseq, individual = individual)
  
  if (!is.null(pattern_riboseq_gene)) {
    f <- grepl(pattern_riboseq_gene, riboseq.deTest$gene_name)
    riboseq.deTest <- riboseq.deTest[f,]
  }
  
  # gene_set_rnaseq
  gene_set_rnaseq <- NULL
  gene_set_rnaseq <- build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("unt48", "tgfb48"), gene_set = gene_set_rnaseq, gname = "unt48.tgfb48.DEtranscription", verbose = T)
  if (f_include_tgfb_vs_cx5461) {
    gene_set_rnaseq <- build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "tgfb48.tgfbCX5461.DEtranscription", verbose = T)
  }
  gene_set_rnaseq <- build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("unt48", "tgfb48", "tgfbCX5461"), gene_set = gene_set_rnaseq, gname = "reversible.transcription.CX", verbose = T)
  
  
  # gene_set
  gene_set <- NULL
  gene_set <- build_gene_set(de = riboseq.deTest, project = project_riboseq, conditions = c("unt48", "tgfb48"), gene_set = gene_set, gname = "unt48.tgfb48.DEtranslation", verbose = T)
  if (f_include_tgfb_vs_cx5461) {
    gene_set <- build_gene_set(de = riboseq.deTest, project = project_riboseq, conditions = c("tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "tgfb48.tgfbCX5461.DEtranslation", verbose = T)
  }
  gene_set <- build_gene_set(de = riboseq.deTest, project = project_riboseq, conditions = c("unt48", "tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "reversible.translation.CX", verbose = T)
  
  # translationONLY
  gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rnaseq, name1="unt48.tgfb48.DEtranslation", name2="unt48.tgfb48.DEtranscription", type="complement", gname="unt48.tgfb48.DEtranslationONLY")
  if (f_include_tgfb_vs_cx5461) {
    gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rnaseq, name1="tgfb48.tgfbCX5461.DEtranslation", name2="tgfb48.tgfbCX5461.DEtranscription", type="complement", gname="tgfb48.tgfbCX5461.DEtranslationONLY")
  }
  gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rnaseq, name1="reversible.translation.CX", name2="reversible.transcription.CX", type="complement", gname="reversible.translationONLY.CX")
  
  gene_set
  
} # build_gene_set_with_rnaseq_riboseq




################################################################################
################################################################################



perform_tests <- function(input, feats.to.test, input_type, gene_set, f_protein_coding_only=TRUE, df_gene_info=NULL, max_gset=Inf, dir_out, fname_prefix='', list_boundary=list("protein_nterminus"=c(1,200), "protein_middle"=c(201,-201), "protein_cterminus"=c(-200,-1), "CDS_nterminus"=c(1,200), "CDS_middle"=c(201,-201), "CDS_cterminus"=c(-200,-1)), verbose=TRUE) {
  
  list_cudf <- list()
  alltestdf  <-  data.frame(); n_feat <- 0
  for (featx  in  feats.to.test) {
    if (verbose) verb("%s\n", featx)
    
    outLocalDir <- sprintf("%s/%s" , dir_out, featx )
    dir.create(outLocalDir, recursive = TRUE , showWarnings = FALSE)
    outLocalPrefix  <-  sprintf("%s/%s.%s", outLocalDir, fname_prefix, featx)
    
    ### sumdf, pcdf
    n_feat <- n_feat+1
    switch(input_type,
           
           # Ybx1
           "rna-binding_protein_mrna_motif"={
             # /home/hkim5/riboprof/out/transcriptome-gtf/m38/Mus_musculus.GRCm38.97.rdna_rn18s.transcriptome.five_prime_utr.fa.gz
             # >ENSMUST00000000001
             tran.fname  <-  get_transcriptome_feature_fasta_filename( species = myData.species , feat = featx )
             fname_fasta <- write_fasta_with_boundary( fname_fasta=tran.fname, list_boundary, featx, outLocalPrefix, alphabet="RNA")
             if (verbose) verb("\tfind_motifs_in_fasta: %s\n", tran.fname)
             Sys.setenv(HOMER="/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/HOMER/bin")
             find_motif_in_fasta( motif=input[['motif']], name=input[['motif_name']], fasta = fname_fasta , outpref = outLocalPrefix  , mismatch = 0 , verbose = verbose )
             if (verbose) verb("\tload_homer_motif_annotation\n")
             cudf  <-  load_homer_motif_annotation( hits = sprintf("%s.found" , outLocalPrefix ) , bed = TRUE )
             if (is.null(cudf)) {
               system(sprintf("rm -f %s", fname_fasta))
               next
             }
             cudf  <-  subset( cudf , strand == "+" )
             suppressMessages(suppressWarnings( sumdf  <-  cudf  %>%  group_by(chr,name)  %>%  summarise(num.hits = n()) ))
             sumdf  <-  as.data.frame(sumdf)
             if (f_protein_coding_only) {
               pcdf <- merge(myData.transcriptome.info$proteinCodingCanonTrans, sumdf, by.x = "transcript_id", by.y = "chr", all.x = TRUE)
             } else {
               pcdf <- merge(myData.transcriptome.info$canonicalTranscripts, sumdf, by.x = "transcript_id", by.y = "chr", all.x = TRUE)
             }
             pcdf$num.hits[is.na(pcdf$num.hits)] <- 0
             # hasCUbox, noCUbox
             #pcdf$feature  <-  ifelse(pcdf$num.hits > 0 , input[['motif_name']] , paste0("no_",input[['motif_name']]))
             pcdf$feature  <-  ifelse( pcdf$num.hits > 0 , featx , paste0("no_",featx))
             pcdf <- merge(pcdf, myData.transcriptome.info$canonicalTranscriptsGeneName[, c("transcript_id", "gene_name"), drop = FALSE], by = "transcript_id", all.x = TRUE)
             # pcdf: data.frame
             # gene_name	gene_id	transcript_id	CDS.length	five_prime_utr.length	three_prime_utr.length	name	num.hits	feature
             stopifnot(nrow(pcdf) == length(unique(pcdf$transcript_id)) )
           },
           
           # Larp1
           "table_mrna_5utr_cds_3utr"={
             if (verbose) verb("\tload %s.\n", input[n_feat])
             df <- read.table(input[n_feat], header=T, row.names=1, sep='\t', quote='')
             strcol <- make.names(tail(strsplit(featx, '_')[[1]],1))
             idx <- which(df[,strcol]==1)
             sym_human <- df[idx,'Gene.Name']
             sym_human <- sym_human[!is.na(sym_human)]
             sym_mouse <- get_orthologous_gene_names( from = "human" , to = "mouse" , conv = orthodf , genes = sym_human )
             cudf <- data.frame(gene_name=sym_mouse, name=featx, num.hits=1)
             sumdf <- cudf
             pcdf  <-  merge( myData.transcriptome.info$canonicalTranscriptsGeneName , sumdf , by = "gene_name", all.x = TRUE)
             pcdf$num.hits[is.na(pcdf$num.hits)]  <-  0
             # mTOR_active_Larp1_5'UTR, no_mTOR_active_Larp1_5'UTR
             pcdf$feature  <-  ifelse( pcdf$num.hits > 0 , featx , paste0("no_",featx))
             # pcdf: data.frame
             # gene_name	gene_id	transcript_id	CDS.length	five_prime_utr.length	three_prime_utr.length	name	num.hits	feature
             stopifnot(nrow(pcdf) == length(unique(pcdf$transcript_id)) )
           },
           
           "table_gene_1st_column"={
             if (verbose) verb("\tload %s.\n", input[n_feat]) 
             sumdf <- read.table(input[n_feat], header=T, row.names=1, sep='\t', quote='')
             sym_human <- sumdf[,1]
             sym_mouse <- get_orthologous_gene_names( from = "human" , to = "mouse" , conv = orthodf , genes = sym_human )
             cudf <- data.frame(gene_name=sym_mouse, name=featx, num.hits=1)
             sumdf <- cudf
             pcdf  <-  merge( myData.transcriptome.info$canonicalTranscriptsGeneName , sumdf , by = "gene_name", all.x = TRUE)
             pcdf$num.hits[is.na(pcdf$num.hits)]  <-  0
             # mTOR_inactive_TR, no_mTOR_inactive_TR
             pcdf$feature  <-  ifelse( pcdf$num.hits > 0 , featx , paste0("no_",featx))
             # pcdf: data.frame
             # gene_name	gene_id	transcript_id	CDS.length	five_prime_utr.length	three_prime_utr.length	name	num.hits	feature
             stopifnot(nrow(pcdf) == length(unique(pcdf$transcript_id)) )
           },
           {}
    )
    
    # store pcdf
    list_cudf[[featx]] <- cudf
    
    # gene_set
    switch(class(gene_set),
           "enrichResult"={
             df_enricher <- as.data.frame(gene_set)
             df_enricher <- df_enricher[order(df_enricher$p.adjust),]
             gene_set_names <- gsub(' ','_',df_enricher$Description)
           },
           { gene_set_names <- names(gene_set$sets) }
    ) 
    
    alldf  <-  data.frame()
    
    n_gset <- 0
    for (gsetx  in  gene_set_names) {
      n_gset <- n_gset + 1
      if (n_gset > max_gset) break
      if (verbose) verb("\t\t%s\n", gsetx)
      
      outLocalDir <-  sprintf("%s/%s/%s", dir_out , featx, gsetx)
      dir.create(outLocalDir ,  recursive = TRUE , showWarnings = FALSE)
      outLocalPrefix  <-  sprintf("%s/%s.%s.%s" , outLocalDir, fname_prefix, featx, gsetx )
      
      # gtab: define universe (protein coding genes detected in a platform)
      # mdf: define genes selected, colnames=c('gene_name','direc')
      switch(class(gene_set),
             "enrichResult"={
               if (is.null(df_gene_info)) {
                 gtab <- data.frame(gene_name=pcdf$gene_name, biotype='protein_coding')
               } else {
                 gtab <- df_gene_info
               }
               syms <- strsplit(df_enricher[n_gset,'geneID'],'/')[[1]]
               mdf <- data.frame(gene_name=syms, direc=df_enricher[n_gset,'Description'])
             },
             {
               gtab  <-  gene_set$tables[[gsetx]]
               mdf  <-  reshape2::melt( gene_set$sets[[gsetx]] )
             }
      )
      gtab  <-  merge( gtab , pcdf , by = "gene_name" , all.y = TRUE )
      # begin addition by H. Kim
      f <- !is.na(gtab$gene_name)
      gtab <- gtab[f, , drop = F]
      # end addition
      colnames(mdf)  <-  c("gene_name","direc")
      mdf  <-  mdf[ !is.na(mdf$direc)  &  mdf$direc != "notSig" , ,drop=FALSE]
      
      # protein coding, with 3' UTR
      mdf  <-  subset( mdf , gene_name  %in%  gtab$gene_name )
      groups  <-  mdf
      colnames(groups)  <-  c("element","group")
      
      # begin modification by H. Kim
      #tdf  <-  test_hypergeometric_enrichment( universe = gtab , groups = groups , prop = "has.CUbox" , element = "gene_name" , strict = TRUE )
      classes <- gtab[, c("gene_name", "feature")]
      colnames(classes) <- c("element", "group")            
      tdf <- test_hypergeometric_enrichment(universe = gtab$gene_name, classes = classes, groups = groups)            
      # end modification
      tdf$geneset  <-  gsetx
      tdf$feature  <-  featx
      alldf  <-  rbind( alldf , tdf )
      
      annotdf  <-  merge( gtab ,  mdf , by = "gene_name" , all = TRUE )
      fname  <-  sprintf("%s.annotation.txt" , outLocalPrefix)
      write.table( annotdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
      
    } # gsetx
    
    
    alldf  <-  subset( alldf , class == featx )
    
    # begin modification by H. Kim
    #alldf  <-  alldf[ order(alldf$p.over.represented) , ,drop=FALSE]
    alldf  <-  alldf[ order(alldf$p.overrepresented) , ,drop=FALSE]
    # end modification 
    alltestdf  <-  rbind( alltestdf , alldf )
    
    # save
    outLocalDir <-  sprintf("%s/%s", dir_out, featx )
    dir.create(outLocalDir ,  recursive = TRUE , showWarnings = FALSE)
    outLocalPrefix <- sprintf("%s/%s.%s", outLocalDir, fname_prefix, featx)
    
    fname  <-  sprintf("%s.tests.txt", outLocalPrefix)
    write.table( alldf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
    
    if (exists("fname_fasta")) {
      # delete temporaty files
      system(sprintf("rm -f %s", fname_fasta))
    }
    
  } # featx
  
  ### save
  outLocalDir <-  dir_out
  dir.create(outLocalDir ,  recursive = TRUE , showWarnings = FALSE)
  outLocalPrefix  <-  sprintf("%s/%s", outLocalDir, fname_prefix)
  
  fname  <-  sprintf("%s.tests.txt", outLocalPrefix)
  write.table( alltestdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
  # begin of addition by H. Kim
  if (nrow(alltestdf) == 0 || (!"geneset" %in% colnames(alltestdf))) {
    return(list(tdf=NULL, tmat=NULL, list_cudf=list_cudf))
  }
  # end of addition
  
  ### plot
  tdf  <-  alltestdf
  tdf$minp  <-  pmin(tdf$p.underrepresented , tdf$p.overrepresented)
  tdf$log10p  <-  log10(tdf$minp)
  # begin modification by H. Kim
  #tdf  <-  tdf[ tdf$minp <= 1e-2  &&  tdf$class == "hasCUbox" ,, drop=FALSE]
  #tdf  <-  tdf[ tdf$minp <= 1e-2 ,, drop=FALSE]
  # include all lines instead of selecting enriched lines for drawing
  #tdf$direc  <-  ifelse( tdf$p.underrepresented < 1e-2 , "underrepresented" , "overrepresented" )
  #tdf$dlog10p  <-  ifelse( tdf$p.overrepresented < 1e-2 , abs(tdf$log10p) , tdf$log10p )
  tdf$direc  <-  ifelse( tdf$p.overrepresented > tdf$p.underrepresented, "underrepresented" , "overrepresented" )
  tdf$dlog10p  <-  ifelse( tdf$p.overrepresented < tdf$p.underrepresented, abs(tdf$log10p), tdf$log10p )
  # when success.population is too small, set nonSig 
  f <- tdf$success.population < 10
  tdf[f, "dlog10p"] <- 0
  # end modification 
  
  tmat  <-  acast( data = tdf , formula = geneset + group ~ feature , value.var = "dlog10p" , fill = 0 )
  # begin of addition by H. Kim
  rownames(tmat) <- gsub('_', ' ', rownames(tmat))
  f <- feats.to.test %in% colnames(tmat)
  feats.to.test <- feats.to.test[f]
  # end of addition
  tmat  <-  tmat[, feats.to.test,drop=FALSE]
  tmat  <-  tmat[ !str_detect(rownames(tmat) , "no_") ,, drop=FALSE]
  f <- tdf$group == gsub('_',' ',tdf$geneset)
  if (all(f)) {
    rownames(tmat) <- sapply(rownames(tmat), function(x) {
      vec <- str_split(x,'_')[[1]]
      tail(vec,1) } )
  }
  
  thresh  <-  15
  tmat[tmat >  thresh]  <-  thresh
  tmat[tmat < -thresh]  <-  -thresh
  
  # save
  fname  <-  sprintf("%s.tests.plot.tdf.txt", outLocalPrefix)
  write.table( tdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
  fname  <-  sprintf("%s.tests.plot.tmat.txt", outLocalPrefix)
  write.table( tmat , file = fname , quote = FALSE , sep = "\t" , row.names = TRUE , col.names = TRUE )
  
  return(list(tdf=tdf, tmat=tmat, list_cudf=list_cudf))
  
  
} # performt_tests

################################################################################
################################################################################

# process_orf_resdf
# input:
#   resdf: generated by test_ORF_classification_read_count_between_conditions()
#   gene_set: 
#   good.comps: vector
#   classx: {"uORF", ... }
# ouput:
#   list$tmat
# usage:
# resdf <- test_ORF_classification_read_count_between_conditions( gtf = NULL , rpf = NULL , species = "m38" , groups = groups, project = "161021" , classes = NULL , outfbase = NULL , verbose = FALSE )
# list_tests <- process_orf_resdf(resdf, gene_set, good.comps=c("unt48 tgfb48" , "tgfb48 tgfbCX5461100nm"), classx="uORF")
process_orf_resdf <- function(resdf, gene_set, good.comps, classx) {
  
  gsdf <- reshape2::melt(gene_set$sets)
  colnames(gsdf) <- c("gene_name", "direc", "geneset")
  gsdf$group <- paste(gsdf$geneset, gsdf$direc, sep = " ")
  gsdf <- gsdf[gsdf$direc != "notSig", , drop = FALSE]
  groups <- gsdf
  
  rdf  <-  merge( resdf , unique(gsdf[,c("geneset","group","direc"),drop=FALSE]) , by = "group" , all.x = TRUE )
  rdf$comparison  <-  paste( rdf$condition1 , rdf$condition2 , sep = " " )
  
  rdf  <-  rdf[rdf$comparison  %in% good.comps ,,drop=FALSE]
  rdf$dlog10p  <-  ifelse( rdf$log2FC > 0 , -log10(rdf$p.value) , log10(rdf$p.value) )
  rdf$dlog10p[rdf$p.value >= 0.05]  <-  0
  rdf  <-  rdf[order(rdf$geneset,rdf$direc),,drop=FALSE]
  
  
  subdf <- rdf[rdf$classification == classx & !is.na(rdf$classification), , drop = FALSE]
  amat <- acast(data = subdf, formula = group ~ comparison, value.var = "dlog10p", fill = 0)
  if (nrow(amat) == 0 || ncol(amat) == 0) {
    return(NULL)
  }
  amat <- amat[intersect(unique(rdf$group), rownames(amat)), good.comps, drop = FALSE]
  
  thresh <- 10
  tmat <- amat
  tmat[tmat > thresh] <- thresh
  tmat[tmat < -thresh] <- -thresh
  
  heat.max.break <- max(abs(tmat)) + 1
  heat.max.break <- max(heat.max.break, 4)
  heat.by.break <- 2 * heat.max.break/100
  
  list(tmat=tmat)
  
} # process_orf_resdf



