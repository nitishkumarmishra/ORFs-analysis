# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2022
# Email:  nitishimtech@gmail.com
# Date: 15 December, 2022
# Script Name: code_name.R
#####################################################################
# ##################  Script Description: ###########################
# This is code for ORF analysis, ORF_analysis_in_progress.R
# Modified by Nitish, still under process
##################### Notes: ########################################
#source("ORF_analysis_in_progress.R")
#####################################################################
################## SET WORKING DIRECTORY ############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/jupyter codes"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################################
#source("functions_transcriptome.R ") # This is giving some error, I have to check probably reshape::melt in place of melt
#Or it may be calling another function, which I need to put in directory

suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(dplyr)))
Mus_musculus.GRCm39.gtf <- NULL

if(is.null(Mus_musculus.GRCm39.gtf)) {suppressMessages(suppressWarnings(source("make_GRCm39_104_transcript_gtf.R")))}

suppressMessages(suppressWarnings(source("functions_rsem.R")))

# five_prime_utr.fasta <- readDNAStringSet("five_prime_utr.bed.fasta")
# names(five_prime_utr.fasta) <- 
#   names(five_prime_utr.fasta) %>%
#   strsplit(., split="::",fixed=TRUE) %>%
#   sapply(., '[', 1) 
# 
# tran <- five_prime_utr.fasta
# outfbase <- "five_prime_utr"

rsem.transcriptome.fasta <- readDNAStringSet("Mus_musculus.GRCm39.104.rdna_rn18s.gtf.rsem.transcriptome.transcripts.fa")
names(rsem.transcriptome.fasta) <- 
  names(rsem.transcriptome.fasta) %>%
  strsplit(., split="::",fixed=TRUE) %>%
  sapply(., '[', 1) 

tran <- rsem.transcriptome.fasta
gtf <- Mus_musculus.GRCm39.gtf
outfbase <- "orfs.mm39.Mus_musculus.GRCm39.104.gtf"



find_all_pattern_indices_in_string_set  <-  function( seq , pattern , name = "transcript_id" ) {
  
  ### get hits
  mh1  <-  vmatchPattern(pattern = pattern , subject = seq )
  
  ### long format
  hdf  <-  reshape2::melt(mh1) # melt works on this class, miraculously.
  colnames(hdf)  <-  c("ix","NA","start","end","width")
  hdf  <-  hdf[, c("ix","start","end") ,drop=FALSE]
  
  # name correctly
  tnames  <-  data.frame( ix = 1:length(names(mh1)) , transcript_id = names(mh1) )
  colnames(tnames)[2]  <-  name
  
  hdf  <-  merge( hdf , tnames , by = "ix" , all.x = TRUE )
  hdf$pattern  <-  pattern
  
  return(hdf)
  
} # find_all_pattern_indices_in_string_set


  if (is.null(tran)) { tran  <-  read.transcriptome.feature.fasta( species = species , feat = "transcripts" )  } 
  stopifnot(class(tran) == "DNAStringSet")  # tran
  
  ### find all start/stop codons in sequence
  
  # start
  startdf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "ATG" )
  
  # stop
  tagdf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "TAG" )
  taadf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "TAA" )
  tgadf  <-  find_all_pattern_indices_in_string_set( seq = tran , pattern = "TGA" )
  stopdf  <-  rbind( tagdf , taadf , tgadf )
  
  foundf  <-  rbind( startdf , stopdf )
  foundf$id  <-  paste( foundf$transcript_id , foundf$start , foundf$end )
  foundf$codon.type  <-  ifelse( foundf$pattern == "ATG" , "start" , "stop" )
  foundf$mod3  <-  foundf$start %% 3
  
  ### get all possible start/stop pairs.
  
  
  # get start rank and stop rank per entry.  This is very efficient for ensuring that a putative ORF does not contain another STOP inside of it.
  # Note that this MUST be first sorted by mod3 to ensure that the start/stop rank pairs are within the same coding frame.
  foundf  <-  foundf[ order(foundf$mod3 , foundf$transcript_id , foundf$start , foundf$end ) , ,drop=FALSE]
  foundf$rank.stop  <-  cumsum( foundf$codon.type == "stop" )
  
  # This allows start and stop codons to be matched such that each start codon is ONLY matched with the closest downstream stop codon within the same coding frame (sorted by mod3, above).
  foundf$rank.stop[foundf$codon.type == "start"]  <-  foundf$rank.stop[foundf$codon.type == "start"] + 1 
  
  # separate into start and stop
  f1df  <-  foundf[ foundf$codon.type == "start" , ,drop=FALSE]
  f2df  <-  foundf[ foundf$codon.type == "stop" , ,drop=FALSE]
  
  mdf  <-  merge( f1df , f2df , by = c("transcript_id","mod3","rank.stop") , suffixes = c(".start",".stop") )
  mdf$ORF.length  <-  (mdf$start.stop - mdf$start.start) / 3  # in number of codons, includes start but not stop codon
  
  stopifnot(all(mdf$start.start < mdf$start.stop))
  
  # must be non-trivial
  #mdf  <-  mdf[ mdf$ORF.length > 1 , ,drop=FALSE]
  ### get known start and stop codons 
  
  
  startc  <-  gtf[ gtf$feature == "start_codon" , c("seqname","start","end","feature","frame") , drop=FALSE ]
  stopc  <-  gtf[ gtf$feature == "stop_codon" , c("seqname","start","end","feature","frame") , drop=FALSE ]
  
  startc$id  <-  paste( startc$seqname , startc$start , startc$end )
  stopc$id  <-  paste( stopc$seqname , stopc$start , stopc$end )
  
  stopifnot(length(unique(startc$id)) == nrow(startc))
  stopifnot(length(unique(stopc$id)) == nrow(stopc))
  
  ### get transcript biotype
 
  
  tbiot  <-  get_gtf_attribute_field( gtf = gtf , field = "transcript_biotype" )
  gbiot  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_biotype" )
  tdf  <-  unique(data.frame( transcript_id = gtf$transcript_id , transcript_biotype = tbiot , gene_biotype = gbiot ))
  mdf  <-  merge( mdf , tdf , by = "transcript_id" , all.x = TRUE )
  
  ### identify annotated ORFs
  
  mdf$isCanonicalORF  <-  (mdf$id.start %in% startc$id  &  mdf$id.stop %in% stopc$id)

  #### classify by relationship to canonical
  
  colnames(startc)  <-  paste( colnames(startc) , "canonicalStart" , sep  = "." )
  colnames(stopc)  <-  paste( colnames(stopc) , "canonicalStop" , sep  = "." )
  
  mdf  <-  merge( mdf , startc , by.x = "transcript_id" , by.y = "seqname.canonicalStart" , all.x = TRUE )
  mdf  <-  merge( mdf , stopc , by.x = "transcript_id" , by.y = "seqname.canonicalStop" , all.x = TRUE )
  
  # classify
  mdf$mod3.canonical  <-  mdf$start.canonicalStart %% 3
  logi.sameFrameAsCanonical  <-  mdf$mod3  ==  mdf$mod3.canonical
  
  ### Classification groups as defined by Ingolia "Ribosome profiling: new views of translation, from single codons to genome scale" 2014, Figure 3.
  
  
  logi.uORF  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  <  mdf$start.canonicalStart)
  logi.overlapping.uORF  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  >  mdf$start.canonicalStart)  &  (mdf$start.stop < mdf$start.canonicalStop)  &  !logi.sameFrameAsCanonical
  logi.alt.ORF  <-  (mdf$start.start  >  mdf$start.canonicalStart)  &  (mdf$start.stop  <  mdf$start.canonicalStop)  &  !logi.sameFrameAsCanonical
  logi.truncation  <-  (mdf$start.start  >  mdf$start.canonicalStart)  &  (mdf$start.stop  ==  mdf$start.canonicalStop)  &  logi.sameFrameAsCanonical
  logi.extension  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  ==  mdf$start.canonicalStop)  &  logi.sameFrameAsCanonical
  
  # i've generalized a few more
  #prot.cod.biotypes  <-  c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_M_gene","IG_V_gene","IG_Z_gene","nonsense_mediated_decay","nontranslating_CDS","non_stop_decay","polymorphic_pseudogene","protein_coding","TR_C_gene","TR_D_gene","TR_gene","TR_J_gene","TR_V_gene")
  
  logi.overlapping.uORF.all  <-  (mdf$start.start  <  mdf$start.canonicalStart)  &  (mdf$start.stop  >  mdf$start.canonicalStart)  &  (mdf$start.stop < mdf$start.canonicalStop)
  logi.dORF  <-  (mdf$start.start  >  mdf$start.canonicalStop)
  logi.overlapping.dORF  <-  (mdf$start.start > mdf$start.canonicalStart)  &  (mdf$start.start < mdf$start.canonicalStop)  &  (mdf$start.stop > mdf$start.canonicalStop) 
  logi.nc  <-  !(mdf$transcript_biotype %in% c("protein_coding"))
  logi.internal  <-  (mdf$start.start  >  mdf$start.canonicalStart)  &  (mdf$start.stop  <  mdf$start.canonicalStop)  &  logi.sameFrameAsCanonical
  
  # Modification by Nitish on December 20, 2022. Add Rendleman et al. classification
  # Add start-stop uORFs/dORFs based on https://www.biorxiv.org/content/10.1101/2021.07.26.453809v2.full 
  # I will use only uORFs/dORFs for start-stop ORF analysis
  logi.ss.uORF  <-  (logi.uORF == TRUE  &  mdf$ORF.length == 1)
  logi.ss.dORF  <-  (logi.dORF == TRUE  &  mdf$ORF.length == 1)
  
  ## apply
  ## 
  
  
  mdf$classification  <-  NA
  mdf$classification[logi.uORF]  <-  "uORF"
  mdf$classification[logi.overlapping.uORF]  <-  "overlapping uORF"
  mdf$classification[logi.alt.ORF]  <-  "alternative RF"
  mdf$classification[logi.truncation]  <-  "truncation"
  mdf$classification[logi.extension]  <-  "extension"
  mdf$classification[logi.dORF]  <-  "dORF"
  mdf$classification[logi.overlapping.dORF]  <-  "overlapping dORF"
  mdf$classification[logi.overlapping.uORF.all]  <-  "overlapping uORF any frame"
  mdf$classification[logi.nc]  <-  "noncoding ORF"
  mdf$classification[logi.internal]  <-  "internal ORF"
  mdf$classification[mdf$isCanonicalORF]  <-  "annotated ORF"
  # modification by Nitish
  mdf$classification[logi.ss.uORF] <- "start stop uORF"
  mdf$classification[logi.ss.dORF] <- "start stop dORF"
  
  subg  <-  unique(gtf[,c("gene_name","transcript_id","gene_id"),drop=FALSE])
  mdf  <-  merge( mdf , subg , by = "transcript_id" , all.x = TRUE )
  
  if(!is.null(outfbase)) {
    fname  <-  sprintf("%s.orfs.txt" , outfbase )
    write.table( mdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
    system(sprintf("gzip  -f  %s" , fname ))
    # bed
    orbed  <-  data.frame( chr = mdf$transcript_id , start = mdf$start.start - 1 , end = mdf$end.stop , classification = mdf$classification ) # , paste( mdf$id.start , mdf$id.stop ) , mdf$mod3 , strand = "+" )
    #orbed  <-  cbind( orbed , mdf[,c("pattern.start","pattern.stop","transcript_biotype","gene_biotype","isCanonicalORF","start.canonicalStart","end.canonicalStop","mod3.canonical","classification"),drop=FALSE] )
    orbed <- orbed[order(orbed[,1], orbed[,2], orbed[,3]),] # substiture of sort at line 248 in accumulate_rpf_counts_per_putative_ORFs.sh
    # sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o $TMPWF.orfs.bed  $TMPWF.orfs.bed
    fname  <-  sprintf("%s.orfs.bed" , outfbase )
    write.table( orbed , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
    system(sprintf("gzip  -k -f  %s" , fname ))
    
  } # outfbase
# find_all_ORFs_in_transcriptome

  # Print  start-stop uORFs/dORFs head
  orbed %>%
    filter(str_detect(classification, "start stop uORF")) %>%
    head()
  
  orbed %>%
    filter(str_detect(classification, "start stop dORF")) %>%
    head()
  
########################################################################################
########################################################################################


  # test_ORF_classification_read_count_between_conditions_per_gene  <-  function( gtf = NULL , rpf = NULL , species = NULL , project , outfbase = NULL , verbose = FALSE ) { 
  #   
  #   project  <-  process_project_input( project = project )
  #   
  #   ### load
  #   if (verbose) { verb("\tload.\n") }
  #   
  #   ### rpf
  #   if (is.null(rpf)) {
  #     rpf  <-  load_accumulated_RPF_count_per_ORF( project = project , classAccum = TRUE , verbose = verbose )
  #   } # rpf
  #   stopifnot(is.data.frame(rpf))
  #   
  #   if (!("condition" %in% colnames(rpf))) {
  #     if (verbose) { verb("\t\tcond.\n") }
  #     rpf  <-  merge( rpf , unique(project[,c("runid","condition"),drop=FALSE]) , by = "runid" )
  #   }
  #   
  #   ### gtf
  #   if (is.null(gtf)) {
  #     if (verbose) { verb("\t\tload gtf.\n") }
  #     gtf  <-  load_transcriptome_gtf( species = species )
  #     gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
  #   } # proteinCodingTranscripts
  #   
  #   ### prepare gtf
  #   if (verbose) { verb("\t\tprepare gtf.\n") }
  #   gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene" )
  #   gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_id" )
  #   gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_biotype" )
  #   gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_biotype" )
  #   
  #   
  #   #### merge gene name
  #   if (verbose) { verb("\tmerge gene name.\n") }
  #   
  #   rpf  <-  merge( rpf , unique(gtf[,c("gene","transcript_id"),drop=FALSE]) , all.x = TRUE )
  #   
  #   conpairs  <-  combn( unique(project$condition) , 2)
  #   
  #   uniq.classes  <-  unique(rpf$classification)
  #   uniq.classes  <-  uniq.classes[!is.na(uniq.classes)]
  #   show(uniq.classes)
  #   
  #   pctids  <-  unique(gtf$transcript_id[gtf$gene_biotype == "protein_coding"  &  gtf$transcript_biotype == "protein_coding"])
  #   
  #   resdf  <-  data.frame()
  #   
  #   ### summarize by gene, class, protein coding
  #   pcdf  <-  rpf[ rpf$transcript_id %in% pctids  ,,drop=FALSE]
  #   pcdf  <-  pcdf  %>%  group_by(gene,classification,runid,condition)  %>%  summarise(total.count = sum(count))
  #   pcdf  <-  as.data.frame(pcdf)
  #   
  #   anndf  <-  pcdf[ pcdf$classification == "annotated ORF"  &  !is.na(pcdf$classification) ,,drop=FALSE]
  #   anndf$classification  <-  NULL
  #   colnames(anndf)[colnames(anndf) == "total.count"]  <-  "annotated.ORF.count"
  #   pcdf  <-  merge( pcdf , anndf , by = c("gene","runid","condition") , all = TRUE )
  #   pcdf$log.odds  <-  log2(pcdf$total.count) - log2(pcdf$annotated.ORF.count)
  #   logi.zn  <-  pcdf$total.count == 0  |  is.na(pcdf$total.count)
  #   logi.zd  <-  pcdf$annotated.ORF.count == 0  |  is.na(pcdf$annotated.ORF.count)
  #   pcdf$log.odds[logi.zn & !logi.zd]  <-  log2(pcdf$total.count[logi.zn & !logi.zd] + 1e-6) - log2(pcdf$annotated.ORF.count[logi.zn & !logi.zd])
  #   pcdf$log.odds[!logi.zn & logi.zd]  <-  log2(pcdf$total.count[!logi.zn & logi.zd]) - log2(pcdf$annotated.ORF.count[!logi.zn & logi.zd] + 1e-6 )
  #   
  #   
  #   for (classx  in  uniq.classes) {
  #     if (classx == "annotated ORF") { next }
  #     if (verbose) { verb("\t\t\t%s\n", classx) }
  #     
  #     for (cx  in  1:ncol(conpairs)) {
  #       cond1  <-  conpairs[1,cx]
  #       cond2  <-  conpairs[2,cx]
  #       if (verbose) { verb("\t\t\t\t%s  %s\n", cond1, cond2) }
  #       
  #       subdf  <-  pcdf[ pcdf$classification == classx  &  !is.na(pcdf$classification)  &  pcdf$condition %in% c(cond1,cond2) ,,drop=FALSE]
  #       if (nrow(subdf) < 2  ||  length(unique(subdf$condition)) < 2) { next }
  #       
  #       
  #       ### matrix
  #       if (verbose) { verb("\t\t\t\tmatrix.\n") }
  #       
  #       gmat  <-  acast( data = subdf , formula = gene ~ runid , value.var = "log.odds" )
  #       gmat  <-  gmat[complete.cases(gmat),,drop=FALSE]
  #       if (nrow(gmat) == 0) { next }
  #       
  #       ridc1  <-  project$runid[project$condition == cond1]
  #       logi.c1  <-  colnames(gmat) %in% ridc1
  #       tres  <-  apply( gmat , 1 , FUN = function(r) t.test( x = r[logi.c1] , y = r[!logi.c1] , alternative = "two.sided" )$p.value )
  #       tres  <-  data.frame( gene = rownames(gmat) , p.value = tres )
  #       
  #       tcdf  <-  dcast( data = subdf , formula = gene ~ runid , value.var = "total.count" , fill = 0 )
  #       aodf  <-  dcast( data = subdf , formula = gene ~ runid , value.var = "annotated.ORF.count" , fill = 0 )
  #       lodf  <-  dcast( data = subdf , formula = gene ~ runid , value.var = "log.odds" , fill = 0 )
  #       colnames(lodf)  <-  paste( colnames(lodf) , "log.odds" , sep="." )
  #       
  #       
  #       ### summarize
  #       if (verbose) { verb("\t\t\t\tsummarize.\n") }
  #       
  #       mdf  <-  merge( tcdf , aodf , by = "gene" , suffixes = c(sprintf(".%s",classx) , ".CDS" ) , all = TRUE )
  #       mdf  <-  merge( mdf , lodf , by.x = "gene" , by.y = "gene.log.odds" , all = TRUE )
  #       mdf  <-  merge( mdf , tres , by = "gene" , all = TRUE )
  #       
  #       
  #       ### mean, sd
  #       if (verbose) { verb("\t\t\t\tmean, sd.\n") }
  #       
  #       sdf  <-  subdf  %>%  group_by(gene,classification,condition)  %>%  summarise(mean.log.odds = mean(log.odds) , sd.log.odds = sd(log.odds) )
  #       sdf  <-  as.data.frame(sdf)
  #       modf  <-  dcast( data = sdf , formula = gene ~ condition , value.var = "mean.log.odds")
  #       modf$log2FC  <-  modf[[cond2]] - modf[[cond1]]
  #       sodf  <-  dcast( data = sdf , formula = gene ~ condition , value.var = "sd.log.odds")
  #       
  #       smdf  <-  merge( modf , sodf , by = "gene" , suffixes = c(".mean",".sd") , all = TRUE )
  #       mdf  <-  merge( mdf , smdf , by = "gene" , all = TRUE )
  #       
  #       
  #       if (!is.null(outfbase)) {
  #         ### write
  #         if (verbose) { verb("\t\t\t\twrite.\n") }
  #         
  #         fname  <-  sprintf("%s.%s.%s--%s.summary.txt", outfbase , classx , cond1 , cond2 )
  #         write.table( mdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
  #       } # outfbase
  #       
  #       rdf  <-  mdf[,c("gene","log2FC","p.value"),drop=FALSE]
  #       rdf$classification  <-  classx
  #       rdf$condition1  <-  cond1
  #       rdf$condition2  <-  cond2
  #       resdf  <-  rbind( resdf , rdf )
  #     } # cx
  #   } # classx
  #   
  #   return(resdf)	
  #   
  # } # test_ORF_classification_read_count_between_conditions_per_gene
  # 
