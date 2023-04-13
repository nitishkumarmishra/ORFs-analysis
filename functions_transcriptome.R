
# # # 			source("/home/map2085/riboprof/src/functions_transcriptome.R")









##########  packages

library(plyr)
library(dplyr)
library(methods)
library(Biostrings)
#library(Rsamtools)
#library(GenomicRanges)
library(kpmt)
library(data.table)

#source("/home/map2085/riboprof/src/matt.essentials.R")
source("/home/map2085/riboprof/src/functions_ucsc.R")
source("/home/map2085/riboprof/src/functions_project.R")
source("/home/map2085/riboprof/src/functions_enrichment_test.R")


######################################### verb

### use this function for stderr output
verb <- function(...) cat(sprintf(...), sep='', file=stderr())

















######################################################
######################################################
######################################################
######################################################
#######################  load  #######################
######################################################
######################################################
######################################################
######################################################




read_gtf  <-  function( file , header = FALSE ) {
#	in.gtf  <-  read.table( file = file , sep = "\t" , header = FALSE ,  row.names = NULL , stringsAsFactors = FALSE , quote = "" ,
#				colClasses = c( rep("character",3) , rep("numeric",2) , rep("character",4) ) )

	in.gtf  <-  fread( input = sprintf("zcat  -f  %s", file) , sep = "\t" , header = FALSE ,  stringsAsFactors = FALSE , quote = "" ,
				colClasses = c( rep("character",3) , rep("numeric",2) , rep("character",4) ) )

	def.colnames  <-  c("seqname" , "source" , "feature" , "start" , "end" , "score" , "strand" , "frame" , "attribute")
	colnames(in.gtf)  <-  def.colnames

	in.gtf  <-  as.data.frame(in.gtf)

	return(in.gtf)

} # read_gtf






load_transcriptome_gtf  <-  function( species , verbose = FALSE ) {

	### transcriptome
	tran.gtf.dir  <-  matt_get_info( species = species , feat = "transcriptomeGTFDir" )
	tran.gtf.gtf  <-   matt_get_info( species = species , feat = "transcriptomeGTF" )

	fname  <-  sprintf("%s/%s.gz", tran.gtf.dir , tran.gtf.gtf )
	if (verbose) { verb("\t\t\tloading transcriptome gtf: [%s]\n", fname ) }

	indf  <-  read_gtf( file = fname )

	return(indf)

} # load_transcriptome_gtf







load_genome_gtf  <-  function( species ) {

	genome.dir  <-  matt_get_info( species = species , feat = "genomedir" )
	gtfFname  <-  matt_get_info( species = species , feat = "gtf" )

	fname  <-  sprintf("%s/%s.gz" , genome.dir , gtfFname)
	indf  <-  read_gtf( file = fname )

	return(indf)

} # load_genome_gtf










load_all_transcriptome_gtf_feature_fastas  <-  function( species , gtf ) { 

	uniq.feats  <-  unique(gtf$feature)
	uniq.feats  <-  c( uniq.feats , "transcripts" )
	inl  <-  sapply( uniq.feats , FUN = function(featx)  read.transcriptome.feature.fasta( species = species , feat = featx ) , simplify = FALSE )

	return(inl)

} # load_all_transcriptome_gtf_feature_fastas







read.transcriptome.feature.fasta  <-  function( species = NULL , file = NULL , feat )  {

	library(Biostrings)

	### get base transcriptome file name
	if (is.null(file)) {
		file  <-  sprintf("%s/%s.gz", matt_get_info( feat = "transcriptomeGTFdir" , species = species ) , matt_get_info( species = species , feat = "transcriptomeGTF"))
	}

	basef  <-  gsub( "\\.gtf\\.gz$", "" , file)

#	if (feat == "transcript") {
#		featf  <-  paste( basef , "fa.gz" , sep=".")
#	} else {
		featf  <-  paste( basef , feat , "fa.gz" , sep=".")
#	}

	### read fasta
	infa  <-  readDNAStringSet( featf )

	return(infa)
} # read.transcriptome.feature.fasta





get_transcriptome_feature_fasta_filename  <-  function( species , feat ) {

	file  <-  sprintf("%s/%s.gz", matt_get_info( feat = "transcriptomeGTFdir" , species = species ) , matt_get_info( species = species , feat = "transcriptomeGTF"))
	basef  <-  gsub( "\\.gtf\\.gz$", "" , file)
	featf  <-  paste( basef , feat , "fa.gz" , sep=".")

	return(featf)

} # get_transcriptome_feature_fasta_filename






















######################################################
######################################################
######################################################
######################################################
#####################  process  ######################
######################################################
######################################################
######################################################
######################################################






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



extract_and_append_attribute_field_to_gtf  <-  function( gtf , field , num = FALSE , overwrite = FALSE ) {

	if (overwrite  ||  !(field %in% colnames(gtf))) {
		gtf[[field]]  <-  get_gtf_attribute_field( gtf = gtf , field = field , num = num )
	} # write

	return(gtf)

} # extract_and_append_attribute_field_to_gtf







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


prepare_genome_gtf  <-  function( gtf ) {

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

} # prepare_genome_gtf
















append_CDS_length_to_transcript_gtf  <-  function( gtf ) {

	subdf  <-  gtf[gtf$feature == "CDS" ,,drop=FALSE]
	stopifnot(nrow(subdf) == length(unique(subdf$transcript_id)))

	subdf$CDS_length  <-  subdf$end - subdf$start + 1
	subdf  <-  subdf[,c("transcript_id","CDS_length"),drop=FALSE]
	subdf  <-  unique(subdf)
	gtf  <-  merge( gtf , subdf , by = "transcript_id" , all.x = TRUE )
	gtf$CDS_length[is.na(gtf$CDS_length)]  <-  0

	return(gtf)

} # append_CDS_length_to_transcript_gtf


















######################################################
######################################################
######################################################
######################################################
###################  canonical  ######################
######################################################
######################################################
######################################################
######################################################




#' Verify relationships between gene_name, gene_id, and transcript_id for a gtf
#'
#' Surprisingly, multiple gene_ids can pertain to the same gene_name. In light of this sad reality,
#' this function checks if naming relationships between gene_name, gene_id, and transcript_id are
#' one-to-many or many-to-many or one-to-one for various pairs of these IDs.
#'
#' @param  gtf		a data frame of a gtf
#' 
#' @return	nothing.  prints some results to console.
#'
#' @export
validate_gtf_names  <-  function( gtf ) {

	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_name" )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_id" )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_id" )

	gtf  <-  gtf[,c("gene_name","gene_id","transcript_id"),drop=FALSE]
	gtf  <-  unique(gtf)

	### number of gene ids per gene name
	sdf  <-  unique(gtf[,c("gene_name","gene_id"),drop=FALSE])
	cvec  <-  sapply( split( sdf$gene_id , sdf$gene_name ) , length )
	cvec  <-  sort(cvec , decreasing = TRUE)
	if (sum(cvec != 1) > 0) {
		verb("\n\n\nWARNING!!! multiple gene ids per gene name!!!\n")
		show(table(cvec))
		show(head(cvec[cvec > 1]))
	} else {
		verb("One gene id per gene name.\n")
	} 

	#### number of gene names per gene id
	sdf  <-  unique(gtf[,c("gene_name","gene_id"),drop=FALSE])
	cvec  <-  sapply( split( sdf$gene_name , sdf$gene_id ) , length )
	cvec  <-  sort(cvec , decreasing = TRUE)
	if (sum(cvec != 1) > 0) {
		verb("\n\n\nWARNING!!! multiple gene names per gene id!!!\n")
		show(table(cvec))
		show(head(cvec[cvec > 1]))
	} else {
		verb("One gene name per gene id.\n")
	}

	### number of gene ids per transcript id
	sdf  <-  unique(gtf[,c("transcript_id","gene_id"),drop=FALSE])
	cvec  <-  sapply( split( sdf$gene_id , sdf$transcript_id ) , length )
	cvec  <-  sort(cvec , decreasing = TRUE)
	if (sum(cvec != 1) > 0) {
		verb("\n\n\nERROR!!!! multiple gene ids per transcript id!!!\n")
		show(table(cvec))
		stop()
	} else {
		verb("One gene id per transcript id.\n")
	} 

} # validate_gtf_names



















#' Get a table of gene name, gene id, transcript id information from a gtf
#'
#' Given a transcriptome gtf, this function will produce a table of gene names, gene ids, and transrcipt ids
#' ranked according to the "canonical transcript" metric.
#'
#' @param  gtf		data frame of transcriptome gtf. If NULL, will load using "species".
#' @param  species	Only used to load transcriptome gtf from file if gtf is NULL.
#'
#' @return	data frame of gene name, gene id, transcript id, biotypes, length, etc.
#'
#' @export
get_gene_transcript_table  <-  function( gtf = NULL , species = NULL , verbose = FALSE ) {

	if (is.null(gtf)) {
		stopifnot(!is.null(species))
		gtf  <-  load_transcriptome_gtf( species = species , verbose = verbose )
	} # 

	### extract info
	if (verbose) { verb("\t\t\textracting info.\n") }
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_name" )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_id" )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_id" )

	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "gene_biotype" )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_biotype" )

	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_length" , num = TRUE )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "five_prime_utr_length" , num = TRUE )
	gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "three_prime_utr_length" , num = TRUE )

	# CDS length
	gtf  <-  append_CDS_length_to_transcript_gtf( gtf = gtf )

	### ranking. 
	gtf  <-  gtf[ order(-gtf$CDS_length , -gtf$five_prime_utr_length , -gtf$three_prime_utr_length , -gtf$transcript_length) ,,drop=FALSE]
	# Now, the canonical versions appear highest in the list.

	good.cols  <-  c("gene_name","gene_id","transcript_id",
			"gene_biotype","transcript_biotype",
			"transcript_length","five_prime_utr_length","three_prime_utr_length","CDS_length")
	gdf  <-  gtf[, good.cols ,drop=FALSE]
	gdf  <-  unique(gdf)

	return(gdf)

} # get_gene_transcript_table









#' Get tables of canonical transcript ids per gene
#'
#' Transcript-to-gene is a many-to-one relationship. gene_id-to-gene_name is many-to-one, too.
#' But for many analyses, we need a single "canonical transcript" per gene so we can talk about CDS length, etc.
#' This function will identify single "canonical transcript" per gene according to the metric:
#' For protein coding genes: the protein coding transcript with the longest CDS, longest 5' UTR, longest 3' UTR (criteria in that order);
#' For non-protein coding genes: the non-protein-coding transcript with longest length.
#' Since gene_name-to-gene_id is one-to-many, there are separate tables depending on whether the single "canonical transcript"
#' is chosen at the level of gene_name or at the level of gene_id.
#'
#' @param  gtf		data frame of a transrcriptome gtf.  if NULL, will be loaded using "species".
#' @param  species	used to load the transcriptome gtf from file if gtf = NULL.
#' @param  table	data frame of ranked transcript summary table from function "get_gene_transcript_table".
#'			If NULL, will be created via "gtf" or "species".
#' 
#' @return	list of canonical transcript id tables.
#'
#' @export
get_canonical_transcript_tables  <-  function( gtf = NULL , species = NULL , table = NULL , verbose = FALSE ) {

	if (is.null(table)) {
		if (is.null(gtf)) {
			gtf  <-  load_transcriptome_gtf( species = species , verbose = verbose )
		} # gtf

		table  <-  get_gene_transcript_table( gtf = gtf , verbose = verbose )
	} # table

	### protein coding
	pct  <-  table[ table$gene_biotype == "protein_coding"  &  table$transcript_biotype == "protein_coding" ,,drop=FALSE] 
	npct  <-  table[ table$gene_biotype != "protein_coding"  &  table$transcript_biotype != "protein_coding"  &  !(table$gene_name %in% pct$gene_name) ,,drop=FALSE]

	# protein coding canonical
	pct  <-  pct[ order(-pct$CDS_length , -pct$five_prime_utr_length , -pct$three_prime_utr_length , -pct$transcript_length) ,,drop=FALSE]
	pctname  <-  pct[ !duplicated(pct$gene_name) ,,drop=FALSE]
	pctgid  <-  pct[ !duplicated(pct$gene_id) ,,drop=FALSE]

	# non-protein coding canonical
	npct  <-  npct[ order(-npct$transcript_length) ,,drop=FALSE]
	npctname  <-  npct[ !duplicated(npct$gene_name) ,,drop=FALSE]
	npctgid  <-  npct[ !duplicated(npct$gene_id) ,,drop=FALSE]

	stopifnot(intersect(pctname$gene_name , npctname$gene_name) == 0)
	stopifnot(intersect(pctname$gene_id , npctname$gene_id) == 0)

	resl  <-  list(	proteinCodingCanonical = list( gene_name = pctname , gene_id = pctgid ) ,
			nonProteinCodingCanonical = list(gene_name = npctname , gene_id = npctgid) )

	return(resl)

} # get_canonical_transcript_tables


















#' Construct a list with all transriptome information needed for computations
#'
#' This is a convenience function which loads a transcriptome gtf and performs
#' several pre-processing steps and extracting commonly used annotations, such as 
#' trasncript to gene tables and canonical transcript tables.
#'
#' @param  species	species
#'
#' @return	list with elements: gtf,trantable,canon,species
#'
#' @export
construct_transcriptome_gtf_data_list  <-  function(species , verbose) {

	if (verbose) { verb("construct_transcriptome_gtf_data_list.\n") }

	### gtf
	gtf  <-  load_transcriptome_gtf(species = species , verbose = verbose )
	gtf  <-  prepare_transcriptome_gtf(gtf = gtf)

	### transcript table
	trandf  <-  get_gene_transcript_table( gtf = gtf , species = species , verbose = verbose )

	### canonical
	canon  <-  get_canonical_transcript_tables( gtf = gtf , species = species , table = trandf , verbose = verbose )

	ret  <-  list( gtf = gtf , trantable = trandf , canon = canon , species = species )
	return(ret)

} # construct_transcriptome_gtf_data_list

























######################################################
######################################################
######################################################
######################################################
######################################################
############  enrichment preparation #################
######################################################
######################################################
######################################################
######################################################
######################################################






#' Subset the universe of genes and get canonical transcript id
#'
#' This function will return the subset of the provided universe of genes according to the indicated
#' subset critera and will provide only the canonical transcript id for the given level.
#'
#' @param  universe		character vector giving the initial universe of genes. Values can be gene names or gene ids.
#' @param  subset		If subset == "protein_coding", then the universe will be intersected with canonical transcripts for protein coding genes.
#'				If subset == "non_protein_coding", then the universe will be intersected with canonical transcripts for non-protein coding genes.
#'				If subset is NULL, then the universe will not be filtered by protein-coding or non-protein coding genes, and instead
#' 
#' @param  level		gene_name  to gene_id is one-to-many. Therefore, canonical transcripts can be identified at the gene_name or gene_id level.
#'				level must be "gene_name" or "gene_id" to indicate at which level a unique canonical transcript shoudl be selected.
#' @param  canon		list returned by "get_canonical_transcript_tables"
#'
#' @return	character vector of restricted universe transcript_ids.
subset_and_get_canonical_transcript  <-  function( universe = NULL , subset , level , canon , verbose = FALSE ) {

	stopifnot(level %in% c("gene_name","gene_id"))
	stopifnot(subset %in% c("protein_coding" , "non_protein_coding") | is.null(subset))


	# subset and universe to trasncript id
	if (subset == "protein_coding") {
		if (verbose) { verb("\t\tsubset protein coding.\n") }

		if (is.null(universe)) {
			universe  <-  unique(canon$proteinCodingCanonical[[level]][[level]])
		}
		
		universe  <-  intersect( universe , canon$proteinCodingCanonical[[level]][[level]] )
		universe  <-  unique(canon$proteinCodingCanonical[[level]]$transcript_id[canon$proteinCodingCanonical[[level]][[level]] %in% universe])

	} else if (subset == "non_protein_coding") {
		if (verbose) { verb("\t\tsubset non-protein coding.\n") }

		if (is.null(universe)) {
			universe  <-  unique(canon$nonProteinCodingCanonical[[level]][[level]])
		}

		universe  <-  intersect( universe , canon$nonProteinCodingCanonical[[level]][[level]] )
		universe  <-  unique(canon$nonProteinCodingCanonical[[level]]$transcript_id[canon$nonProteinCodingCanonical[[level]][[level]] %in% universe])

	} else if (is.null(subset)) {
		if (verbose) { verb("\t\tnull subset.\n") }


		rdf  <-  rbind( canon$proteinCodingCanonical[[level]] , canon$nonProteinCodingCanonical[[level]] )
		rdf  <-  rdf[ !duplicated(rdf[[level]]) ,,drop=FALSE]

		rdf  <-  unique(rdf[, c(level,"transcript_id") ,drop=FALSE])

		if (is.null(universe)) {
			universe  <-  unique(rdf$transcript_id)
		} else {
			universe  <-  unique(rdf$transcript_id[rdf[[level]]  %in%  universe])
		}
	} else {
		verb("\n\n\nERROR!  invalid subset value!!\n\n")
		stop()
	} # protein_coding


	return(universe)

} # subset_and_get_canonical_transcript




















#' Restrict gene sets and universe according to subset criteria and input universe
#'
#' This function will restrict the input universe according to the specified subset criteria (e.g. protein coding only),
#' and then will restrict the gene sets to the restricted universe.
#'
#' @param  sets			data frame with colname matching the 'level' parameter.
#' @param  universe		character vector giving the initial universe of genes. Values are as specified by the 'level' parameter.
#' @param  subset		If subset == "protein_coding", then the universe will be intersected with canonical transcripts for protein coding genes.
#'				If subset == "non_protein_coding", then the universe will be intersected with canonical transcripts for non-protein coding genes.
#'				If subset is NULL, then the universe will not be filtered by protein-coding or non-protein coding genes.
#' 
#' @param  level		"gene_name" or "gene_id" or "transcript_id", specifying which kind of IDs are provided in 'universe' and 'sets'.
#' @param  canon		list returned by "get_canonical_transcript_tables"
#'
#' @return  list with fields 'sets' and 'universe', giving the restricted versions of these inputs.
#'
#' @export
subset_sets_and_universe  <-  function( sets , universe , level , subset , canon = NULL , gtf = NULL , species = NULL , verbose = FALSE ) {

	if (is.null(canon)) {
		canon  <-  get_canonical_transcript_tables( gtf = gtf , species = species , verbose = verbose )
	} # proteinCodingTranscripts

	### subset universe
	if (is.null(subset)) {
		if (verbose) { verb("\t\tnull subset.\n") }

		universe  <-  unique(c(canon$proteinCodingCanonical[[level]][[level]] , canon$nonProteinCodingCanonical[[level]][[level]]))

	} else if (subset == "protein_coding") {
		if (verbose) { verb("\t\tsubset protein coding.\n") }

		if (is.null(universe)) {
			universe  <-  unique(canon$proteinCodingCanonical[[level]][[level]])
		}

		universe  <-  intersect( universe , canon$proteinCodingCanonical[[level]][[level]] )

	} else if (subset == "non_protein_coding") {
		if (is.null(universe)) {
			universe  <-  unique(canon$nonProteinCodingCanonical[[level]][[level]])
		}

		universe  <-  intersect( universe , canon$nonProteinCodingCanonical[[level]][[level]] )

	} else {
		verb("\n\n\nERROR!  invalid subset value!!\n\n")
		stop()
	} # protein_coding

	stopifnot(length(universe)>0)

	### subset sets
	sets  <-  sets[ sets[[level]]  %in%  universe ,,drop=FALSE]
	stopifnot(nrow(sets)>0)

	resl  <-  list(sets = sets , universe = universe)
	return(resl)

} # subset_sets_and_universe












































#' @param  sets		data frame containing gene sets of interest.  column names must include "group" and "gene_name" or "transcript_id" or "gene_id"
#' @param  level	gene_name  to gene_id is one-to-many. Therefore, canonical transcripts can be identified at the gene_name or gene_id level.
#'			level must be "gene_name" or "gene_id" to indicate at which level a unique canonical transcript shoudl be selected.
#' @param  universe	character vector of transcript_ids defining the universe.
#' @param  transcripts	data frame returned by "get_gene_transcript_table" 
#'
#'
#' @return	sets with transcript_id column
#'
#' @export
convert_sets_to_transcript_id_and_restrict_to_universe  <-  function( sets , level , universe , transcripts , verbose = FALSE ) {

	stopifnot(all(universe %in% transcripts$transcript_id)  &&  all(!is.na(universe)))


	### convert gene set to transcript id and restrict
	if (verbose) { verb("\t\tuniverse subset.\n") }

	subdf  <-  transcripts[ transcripts$transcript_id %in% universe ,,drop=FALSE]


	### get transcript ids
	if (!("transcript_id" %in% colnames(sets))) {
		if (verbose) { verb("\t\tmerging to get transcript_id.\n") }
		if (verbose) { show(head(subdf)) ; show(head(sets)) }
		

		nbefore  <-  length(unique(sets[[level]]))
		sets  <-  merge( sets , subdf , by = level )
		nafter  <-  length(unique(sets[[level]]))

		if (nbefore != nafter) {
			verb("WARNING!  %d genes removed by universe restriction.\n", nbefore - nafter )
		}
	} else {
		stopifnot(all(sets$transcript_id  %in%  transcripts$transcript_id))
	} # gene name

	### restrict by universe
	if (verbose) { verb("\t\trestrict.\n") }

	nbefore  <-  length(unique(sets[[level]]))
	sets  <-  sets[ sets$transcript_id %in% universe ,,drop=FALSE]
	nafter  <-  length(unique(sets[[level]]))

	if (nbefore != nafter) {
		verb("WARNING!  %d genes removed by universe restriction.\n" , nbefore - nafter )
	}


	if (!all(!is.na(sets$transcript_id))) {
		if (verbose) { verb("\t\tchecking for NA transcript ids.\n") }

		verb("\n\nshowing bad:\n")
		show(dim(sets))
		show(sum(!complete.cases(sets)))
		#show(sets[!complete.cases(sets),,drop=FALSE])
	}

	stopifnot(all(!is.na(sets$transcript_id)))

	return(sets)

} # convert_sets_to_transcript_id_and_restrict_to_universe







































































################## LEGACY USAGE ONLY
get_transcriptome_features_quant  <-  function(	gtf ,
						canonTrans ,
						feats.fa = NULL ,
						ucsc.fold = NULL ,
						rare.codons = c("TTA","CTA","TCG","GCG") ,
						usage.codon = NULL ) {

	geneFeatQuant  <-  list()

	###### gene level
	geneFeatQuant$gene_id  <-  canonTrans[, c("gene_id","transcript_id")]


	### number of transcripts
	uniq.df  <-  gtf[, c("transcript_id","gene_id")]
	tab.nt  <-  table(uniq.df$gene_id)
	num.tran.df  <-  data.frame( gene_id = names(tab.nt) , num.transcripts = as.vector(tab.nt) )
	geneFeatQuant$gene_id  <-  merge( geneFeatQuant$gene_id , num.tran.df , by = "gene_id" , all.x = TRUE )


	if ( !is.null(rare.codons)  &&  !is.null(usage.codon) ) {
	### rare codon
		rarc  <-  rowSums(usage.codon[, rare.codons])
		rar.df  <-  data.frame( transcript_id = names(rarc) , rareCodonCount = as.vector(rarc) )
		rar.df  <-  merge( rar.df , canonTrans , by = "transcript_id" , all.y = TRUE )
		
		geneFeatQuant$gene_id  <-  merge( geneFeatQuant$gene_id , rar.df[,c("gene_id","rareCodonCount")] , by = "gene_id" , all.x = TRUE )
	} # not NULL rare.codons

	
	for (featx  in  c("five_prime_utr","three_prime_utr","CDS"))  {
		verb("\t\t\t%s\n", featx)
	
		### init
		geneFeatQuant[[featx]]  <-  canonTrans[, c("gene_id","transcript_id")]
	
		### subset
		logi.feat  <-  gtf$feature == featx
		logi.canon  <-  gtf$transcript_id  %in%  canonTrans$transcript_id
		sub.df  <-  gtf[ logi.canon & logi.feat ,]
	
	
		### length
		sub.df$length  <-  sub.df$end - sub.df$start + 1
		sub.df$length[is.na(sub.df$length)]  <-  0
	
		geneFeatQuant[[featx]]  <-  merge( geneFeatQuant[[featx]] , sub.df[,c("gene_id","length")] , by = "gene_id" , all.x = TRUE )
	
		### log2Length
		geneFeatQuant[[featx]]$log2Length  <-  log2(geneFeatQuant[[featx]]$length)
		logi.zero  <-  geneFeatQuant[[featx]]$length == 0
		geneFeatQuant[[featx]]$log2Length[logi.zero]  <-  NA
	
		### exists
		geneFeatQuant[[featx]]$exists  <-  factor(geneFeatQuant[[featx]]$length > 0 , levels = c(TRUE,FALSE))
	
	
		### GC.content
		if (!is.null(feats.fa)) {

			gc.mat  <-  alphabetFrequency( feats.fa[[featx]] , as.prob = TRUE , baseOnly = TRUE )
			gc.mat  <-  as.data.frame(gc.mat)
			rownames(gc.mat)  <-  names(feats.fa[[featx]] )
			gc.df  <-  data.frame( GC = gc.mat$G + gc.mat$C , transcript_id = rownames(gc.mat) )
			gc.df  <-  merge( gc.df ,  canonTrans , by = "transcript_id" )
		
			geneFeatQuant[[featx]]  <-  merge( geneFeatQuant[[featx]] , gc.df[,c("gene_id","GC")] , by = "gene_id" , all.x = TRUE )
			
			### log2GC
			geneFeatQuant[[featx]]$log2GC  <-  log2(geneFeatQuant[[featx]]$GC)
			logi.zero  <-  geneFeatQuant[[featx]]$GC == 0
			geneFeatQuant[[featx]]$log2GC[logi.zero]  <-  NA

		} # feats.fa

		
		if (!is.null(ucsc.fold)) {	
			### MFE
			if (featx  %in%  names(ucsc.fold)) {
				### MFE
				geneFeatQuant[[featx]]  <-  merge( geneFeatQuant[[featx]] , ucsc.fold[[featx]][, c("transcript_id","MFE")] , by = "transcript_id" , all.x = TRUE )
		
				### log2MFE
				geneFeatQuant[[featx]]$log2MFE  <-  log2(-geneFeatQuant[[featx]]$MFE)
				logi.zero  <-  geneFeatQuant[[featx]]$length == 0
				geneFeatQuant[[featx]]$log2MFE[logi.zero]  <-  NA
			} # ucsc.fold
		} # not NULL ucsc.fold
	} # featx

	return(geneFeatQuant)
	
} # get_transcriptome_features_quant






































######################################################
######################################################
######################################################
######################################################
######################################################
##################  feature test  ####################
######################################################
######################################################
######################################################
######################################################
######################################################






prepare_transcript_quantitative_feature_matrix  <-  function( gtf = NULL , fasta = NULL , species = NULL , verbose = FALSE ) {

	if (is.null(gtf)) {
		if (verbose) { verb("\t\tloading transcriptome gtf.\n") }
		gtf  <-  load_transcriptome_gtf( species = species )
		gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
	} # gtf


	if (is.null(fasta)) {
		if (verbose) { verb("\t\tloading transcriptome feature fasta.\n") }
		fasta  <-  load_all_transcriptome_gtf_feature_fastas( species = species , gtf = gtf )
	} # fasta

	### MFE
	if (verbose) { verb("\t\tloading ucsc mfe.\n") }
	mfe  <-  load_ucsc_mrna_folding_info( species = species )

	### length, GC of fasta
	fasta  <-  fasta[names(fasta) %in% c("five_prime_utr","three_prime_utr","CDS","transcripts")]
	
	alldf  <-  data.frame()
	for (featx  in  names(fasta)) {
		if (verbose) { verb("\t\t\t%s\n", featx) }

		gc.mat  <-  alphabetFrequency( fasta[[featx]] , as.prob = TRUE , baseOnly = TRUE )
		gc.mat  <-  as.data.frame(gc.mat)
		rownames(gc.mat)  <-  names(fasta[[featx]] )

		fdf  <-  data.frame( GC = gc.mat$G + gc.mat$C , transcript_id = rownames(gc.mat) , length = width(fasta[[featx]]) )

		### add MFE
		if (featx %in% c("five_prime_utr","three_prime_utr") ) {
			if (verbose) { verb("\t\t\t\tadd mfe.\n") }
			subdf  <-  mfe[ mfe$feature == featx , c("transcript_id","MFE"),drop=FALSE]
			fdf  <-  merge( fdf , subdf , by = c("transcript_id") , all.x = TRUE )
			fdf$MFE[is.na(fdf$MFE)  &  fdf$length == 0]  <-  0
		} # MFE

		### long format
		fdf  <-  melt(fdf , id.vars = "transcript_id" )
		fdf$feature  <-  featx

		alldf  <-  rbind( alldf , fdf )
	} # featx

	
	### to matrix
	if (verbose) { verb("\t\tto matrix.\n") }

	alldf$qfeat  <-  paste( alldf$feature , alldf$variable , sep = " " )
	fmat  <-  acast( data = alldf , formula = transcript_id ~ qfeat , value.var = "value" , fill = 0 )
	
	return(fmat)

} # prepare_transcript_quantitative_feature_matrix















#' Test for rank imbalance of quantative features across gene sets
#'
#' This function will perform the KPMT to determine if the gene sets of interest
#' have unusual distributions of variosu quantitative features of transcripts,
#' including length and GC content of 5' UTR, 3' UTR, and CDS, as well as
#' length and GC content of the CDS.
#'
#' @param  gtf		transcriptome gtf, as obtained from "load_transcriptome_gtf"
#' @param  fasta	transcriptome feature fasta, as obtained from "load_all_transcriptome_gtf_feature_fastas"
#' @param  species	character string giving the species being analyzed.  Only necessary if quant is NULL and gtf, info, or fasta is NULL.
#' @param  quant	numeric matrix, rownames = transcript ids, colnames = features being tested. AS obtained from "prepare_transcript_quantitative_feature_matrix"
#' @param  sets		data frame containing gene sets of interest.  column names must include "group" and "gene_name" or "transcript_id" or "gene_id"
#' @param  universe	character vector defining the universe. THe values of the character vector should correspond to the "level" argument.
#'			i.e. if level == "gene_name", then "universe" should consist of gene names.
#'			If NULL, the rownames of "quant" will be used.
#'			If subset is not NULL, then the universe will be subsetted according to 'subset' argument.
#' @param  subset	If subset == "protein_coding", then the universe will be intersected with canonical transcripts for protein coding genes.
#'			If subset == "non_protein_coding", then the universe will be intersected with canonical transcripts for non-protein coding genes.
#'			if subset is NULL, then the universe will not be filtered by protein-coding or non-protein coding genes.
#' 
#' @param  level	gene_name  to gene_id is one-to-many. Therefore, canonical transcripts can be identified at the gene_name or gene_id level.
#'			level must be "gene_name" or "gene_id" to indicate at which level a unique canonical transcript shoudl be selected.
#'
#'
#' @return	data frame of KPMT results as returned by "kpmt".
#'
#' @export
test_transcriptome_quantitative_features_enrichment  <-  function( sets , universe = NULL , subset = "protein_coding" , level = NULL , gtf = NULL , fasta = NULL , species = NULL , quant = NULL , verbose = FALSE ) {

	stopifnot(level %in% c("gene_name","gene_id"))
	stopifnot(subset %in% c("protein_coding" , "non_protein_coding") | is.null(subset))

	if (is.null(gtf)) {
		if (verbose) { verb("\t\tload gtf.\n") }
		gtf  <-  load_transcriptome_gtf( species = species )
		gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
	} # proteinCodingTranscripts


	if (is.null(quant)) { 
		quant  <-  prepare_transcript_quantitative_feature_matrix( gtf = gtf , fasta = fasta , species = species , verbose = verbose )
	} # quant

	# transcript and gene name tables.
	trandf  <-  get_gene_transcript_table( gtf = gtf , verbose = verbose )
	canon  <-  get_canonical_transcript_tables( gtf = gtf , verbose = verbose )


	### universe
	# subset and universe to trasncript id
	if (verbose) { verb("\t\tinferring universe.\n") }

	universe  <-  subset_and_get_canonical_transcript( universe = universe , subset = subset , level = level , canon = canon , verbose = verbose )

	if (verbose) { verb("\t\tshowing universe:\n") }
	show(head(universe))
	show(length(universe))


	### restrict
	if (verbose) { verb("\t\trestrict.\n") }
	quant  <-  quant[rownames(quant) %in% universe ,,drop=FALSE]
	universe  <-  rownames(quant)  # important for consistency


	### convert gene set to transcript id and restrict
	if (verbose) { verb("\t\trestrict sets.\n") }
	sets  <-  convert_sets_to_transcript_id_and_restrict_to_universe( sets = sets , level = level , universe = universe , transcripts = trandf , verbose = verbose ) 


	#### split
	if (verbose) { verb("\t\tsplit.\n") }
	groups  <-  split(sets$transcript_id , sets$group )


	### kpmt
	if (verbose) { verb("\t\tkpmt.\n") }

	#show(head(quant))
	#show(dim(quant))
	#show(head(groups))
	
	resdf  <-  kpmt( pop = quant , obs = groups , tail = "two-sided" , verbose = FALSE) # verbose )
	resl  <-  list(kpmt = resdf)

	return(resl)

} # test_transcriptome_quantitative_features_enrichment
















#' Test for overrepresentation of start and stop codon usage
#' 
#' This function will perform hypergeometric tests to determine over- and under-
#' representation of start codon usage and stop codon usage, separately, for given gene sets.
#'
#' @param  gtf		transcriptome gtf, as obtained from "load_transcriptome_gtf"
#' @param  species	character string giving the species being analyzed.  Only necessary if quant is NULL and gtf, info, or fasta is NULL.
#' @param  sets		data frame containing gene sets of interest.  column names must include "group" and "gene_name" or "transcript_id" or "gene_id"
#' @param  universe	character vector defining the universe. The values of the character vector should correspond to the "level" argument.
#'			i.e. if level == "gene_name", then "universe" should consist of gene names.
#'			If NULL, the rownames of "quant" will be used.
#' @param  level	gene_name  to gene_id is one-to-many. Therefore, canonical transcripts can be identified at the gene_name or gene_id level.
#'			level must be "gene_name" or "gene_id" to indicate at which level a unique canonical transcript shoudl be selected.
#' @return	data.frame of hypergeometric test enrichment results in format returned by "test_hypergeometric_enrichment"
#'
#' @export
test_start_stop_codon_enrichment  <-  function( gtf = NULL , species = NULL , sets , universe = NULL , level , verbose = FALSE ) {

	stopifnot(level %in% c("gene_name","gene_id"))

	if (is.null(gtf)) {
		if (verbose) { verb("\t\tload gtf.\n") }
		gtf  <-  load_transcriptome_gtf( species = species )
	} # proteinCodingTranscripts


	### prepare
	if (verbose) { verb("\tprepare.\n") }
		
	gtf  <-  prepare_transcriptome_gtf( gtf = gtf )	
	gtf$gene_biotype  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_biotype" )
	gtf$transcript_biotype  <-  get_gtf_attribute_field( gtf = gtf , field = "transcript_biotype" )
	#gtf$TSL  <-  get_gtf_attribute_field( gtf = gtf , field = "transcript_support_level" , num = TRUE )


	if (verbose) { verb("\tload start/stop codon fa.\n") }
	afa  <-  read.transcriptome.feature.fasta( species = species , feat = "start_codon" )
	ofa  <-  read.transcriptome.feature.fasta( species = species , feat = "stop_codon" )

	if (verbose) { verb("\tadd start/stop codon.\n") }
	adf  <-  data.frame( transcript_id = names(afa) , codon = as.character(afa) , type = "start_codon" )
	odf  <-  data.frame( transcript_id = names(ofa) , codon = as.character(ofa) , type = "stop_codon" )
	aodf  <-  rbind( adf , odf )
	aodf  <-  dcast( data = aodf , formula = transcript_id ~ type , value.var = "codon" , fill = NA )
	#show(aodf[!complete.cases(aodf),,drop=FALSE])
	#stopifnot(all(complete.cases(aodf)))

	gtf  <-  merge( gtf , aodf , by = "transcript_id" , all.x = TRUE )
	#gtf  <-  gtf[ gtf$transcript_biotype == "protein_coding"  &  gtf$gene_biotype == "protein_coding" ,,drop=FALSE]
	
	### tables
	if (verbose) { verb("\ttables.\n") }
	trandf  <-  get_gene_transcript_table( gtf = gtf , verbose = verbose )
	canon  <-  get_canonical_transcript_tables( gtf = gtf , species = species , verbose = verbose )


	### universe
	# subset and universe to trasncript id
	if (verbose) { verb("\t\tinferring universe.\n") }

	universe  <-  subset_and_get_canonical_transcript( universe = universe , subset = "protein_coding" , level = level , canon = canon , verbose = verbose )


	### convert gene set to transcript id and restrict
	if (verbose) { verb("\t\trestrict sets.\n") }

	sets  <-  convert_sets_to_transcript_id_and_restrict_to_universe( sets = sets , level = level , universe = universe , transcripts = trandf , verbose = verbose ) 



	### test
	if (verbose) { verb("\ttest.\n") }

	groups  <-  data.frame( group = sets$group , element = sets$transcript_id )

	startdf  <-  unique(gtf[ gtf$transcript_id %in% universe  ,c("transcript_id","start_codon") ,drop=FALSE])
	stopdf  <-  unique(gtf[ gtf$transcript_id %in% universe  ,c("transcript_id","stop_codon") ,drop=FALSE])

	startdf$codon  <-  paste( "start" , startdf$start_codon , sep = " " )
	stopdf$codon  <-  paste( "stop" , stopdf$stop_codon , sep = " " )

	startdf  <-  startdf[,c("transcript_id","codon"),drop=FALSE]
	stopdf  <-  stopdf[,c("transcript_id","codon"),drop=FALSE]
	
	codondf  <-  rbind( startdf , stopdf )
	codondf$element  <-  codondf$transcript_id
	codondf$group  <-  codondf$codon
	codondf  <-  unique(codondf[,c("element","group"),drop=FALSE])
	
	resdf  <-  test_hypergeometric_enrichment( groups = groups , universe = universe , classes = codondf )
	resl  <-  list(enrich = resdf)

	return(resl)

} # test_start_stop_codon_enrichment







































#' Make a transcriptome-level GTF annotation file from Ensembl genome GTF annotation file
#'
#' This function will take a GTF file from Ensembl and create a transcriptome-level GTF file.
#' This means that the "seqname" field will correspond to individual transcripts, and the "start" and "end"
#' coordinate fields will be relative to the transcript, not the genome. Thus, all entries will have 
#' strand = "+" since the transcripts are all taken to be RNA, i.e. single-stranded.
#'
#' @param  speecies		species whose transcriptome GTF should be built.
#' @param  gtf			filename of a gtf file. can be gzipped. inferred from "species" if missing.
#' @param  outfbase		output filename prefix.  inferred if missing
#'
#' @return	nothing.  result is written to file.
#'
#' @export
make_eukaryotic_transcript_gtf  <-  function( species = NULL , gtf = NULL , outfbase , debug = FALSE , verbose = TRUE ) {


	######################################## DEBUG
	if (debug) {
		gtf  <-  "/home/map2085/riboprof/data/human/Homo_sapiens.GRCh38.84.gtf.gz"
		outfbase  <-  "/home/map2085/riboprof/work/gtf.info.out.debug"
	}
	######################################## DEBUG
	
	############## NOTE
	# 5' & 3' UTRs are nested inside of the first & last exons.
	
	##########  packages
	library(stringr)
	library(plyr)
	library(dplyr)
	library(ggplot2)
	
	
	
	
	######################## load
	verb("\tload.\n")
	
	myData.gtf  <-  read_gtf( file = gtf )
	
	myData.colnames  <-  colnames(myData.gtf)
	
	myData.tranGTF  <-  data.frame()
	
	
	# DEBUG
	if (debug) {
		myData.SAVE  <-  myData.gtf  
	} # debug
	
	
	
	
	### feature per feature
	verb("\tfeature per feature.\n")
	
	myData.gtf$transcript_id  <-  get_gtf_attribute_field( gtf = myData.gtf , field = "transcript_id" )
	myData.gtf$exon_number  <-  get_gtf_attribute_field( gtf = myData.gtf , field = "exon_number" , num = TRUE )

	stopifnot(!any(is.na(myData.gtf$transcript_id)))
	stopifnot(all(!is.na(myData.gtf$exon_number) & myData.gtf$feature == "exon"))
	
	
	### feat width
	myData.gtf$width  <-  myData.gtf$end - myData.gtf$start + 1
	stopifnot(all(!is.na(myData.gtf$width)))
	stopifnot(all(myData.gtf$width > 0))
	
	
	
	
	###### transcript length
	verb("\ttranscript length.\n")
	
	# Note that "exon" is anything other than intron, i.e. anything transcribed, even pseudogenes.  So the exons = 5' UTR & CDS & 3' UTR.
	# Thus, summing up exon lengths will give you the full transcript length.
	sub.feat  <-  myData.gtf[ myData.gtf$feature %in% c("exon") , ,drop=FALSE]
	sub.feat  <-  sub.feat  %>%  group_by(transcript_id)  %>%  summarise(transcript_length = sum(width) )
	sub.feat  <-  as.data.frame(sub.feat)
	
	myData.gtf  <-  merge( myData.gtf , sub.feat , by = "transcript_id" , all.x = TRUE)
	stopifnot(!any(is.na(myData.gtf$transcript_length)))
	stopifnot(all(myData.gtf$transcript_length > 0))
	
	
	
	
	
#	#### basic filter
#	verb("\tbasic filter.\n")
#	
#	#myData.gtf  <-  myData.gtf[ !is.na(myData.gtf$transcript_id)  &  !is.na(myData.gtf$transcript_length)  &  myData.gtf$transcript_length > 0 ,]
#	#myData.gtf  <-  myData.gtf[ !is.na(myData.gtf$transcript_id) , ,drop=FALSE]
#	stopifnot(!any(is.na(myData.gtf$transcript_id)))
#	
#	if ( any(is.na(myData.gtf$transcript_length))  ||  any(myData.gtf$transcript_length == 0) ) {
#		verb("\n\n\nERROR!  NA or 0 transcript lengths!!!!\n")
#	
#		verb("NA lengths:\n")
#		show(head(myData.gtf[ is.na(myData.gtf$transcript_length) , ,drop=FALSE ] ))
#	
#		verb("\n\n\n0 lengths:\n")
#		show(head(myData.gtf[ myData.gtf$transcript_length == 0 , ,drop=FALSE ] ))
#	
#		stop()
#	} # check
	
	
	# DEBUG
	if (debug) {
		myData.SAVE2  <-  myData.gtf  
	} # debug
	
	
	
	
	
	
	
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
			sub.utr  <-  sub.utr[,c("transcript_id","width"),drop=FALSE]
			sub.utr  <-  sub.utr  %>%  group_by(transcript_id)  %>%  summarise(utr.width = sum(width))
			sub.utr  <-  as.data.frame(sub.utr)
			sub.utr[[widfield]]  <-  sub.utr$utr.width
			sub.utr$utr.width  <-  NULL
	
			myData.gtf  <-  merge( myData.gtf , sub.utr , by = "transcript_id" , all.x = TRUE)
			logi.na  <-  is.na(myData.gtf[[widfield]])
			myData.gtf[[widfield]][logi.na]  <-  0
		} # has this UTR
	} # utrx
	
		
	
	
	
	############ cumulative relative position
	verb("\tcumulative relative position.\n")
	
	sub.exon  <-  myData.gtf[ myData.gtf$feature == "exon" , ]
	sub.exon  <-  sub.exon[ order(sub.exon$transcript_id , sub.exon$exon_number) , ,drop=FALSE]
	
	### orient start
	verb("\t\torient start.\n")
	
	sub.exon$exon.orient.start  <-  ifelse( sub.exon$strand == "+" , sub.exon$start , sub.exon$end )
	
	
	### cumsum
	verb("\t\tcumsum.\n")
	
	sub.exon  <-  sub.exon  %>%  group_by(transcript_id)  %>%  mutate(cs.pos = cumsum(width) - width)
	sub.exon  <-  as.data.frame(sub.exon)
	
	#logi.na  <-  is.na(sub.exon$five_prime_utr.width)
	#sub.exon$cs.pos[!logi.na]  <-  sub.exon$cs.pos[!logi.na] + sub.exon$five_prime_utr.width[!logi.na]
	
	myData.gtf  <-  merge( sub.exon[,c("transcript_id","exon_number","cs.pos","exon.orient.start"),drop=FALSE] , myData.gtf , by = c("transcript_id","exon_number") , all.y = TRUE , fill = 0 )
	
	
	
	
	
	### 3' UTR
	verb("\t\t3.\n")
	
	logi.utr  <-  myData.gtf$feature == "three_prime_utr"
	if (any(logi.utr)) {
		myData.gtf$cs.pos[logi.utr]  <-  myData.gtf$transcript_length[logi.utr] - myData.gtf$three_prime_utr.width[logi.utr]
	} # exists
	
	# 5' UTR # is already = 0 (by fill)
	
	
	
	
	
	
	
	#### add info to attributes
	verb("\tadd info to attributes.\n")
	
	p.3  <-  paste( 'three_prime_utr_length "' , myData.gtf$three_prime_utr.width , '"; ' , sep="")
	p.5  <-  paste( 'five_prime_utr_length "' , myData.gtf$five_prime_utr.width , '"; ' , sep="")
	p.L  <-  paste( 'transcript_length "' , myData.gtf$transcript_length , '"; ' , sep="")
	
	myData.gtf$attribute  <-  paste( myData.gtf$attribute , " " , p.3 , p.5 , p.L , sep = "")
	
	
	
	
	
	
	# feature = "transcript","three_prime_utr","five_prime_utr" have "NA" for exon_number and hence
	# will lack exon.cs.   This is easy to infer though
	
	
	
	########## to transcript coords.
	verb("\tto transcript coords.\n")
	
	
	### 5' UTR
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
		sub.feat$start  <-  sub.feat$transcript_length - sub.feat$three_prime_utr.width + 1
		sub.feat$end  <-  sub.feat$transcript_length
		
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
	sub.feat$start  <-  sub.feat$transcript_length - sub.feat$three_prime_utr.width - 2
	sub.feat$end  <-  sub.feat$transcript_length - sub.feat$three_prime_utr.width
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
	sub.feat$end  <-  sub.feat$transcript_length - sub.feat$three_prime_utr.width
	
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
	
	
	
	###### stop codon
	
	
	
	
	
	###### start/stop codon
	#verb("\t\tstart-stop codon.\n")
	#
	#sub.feat  <-   myData.gtf[ myData.gtf$feature %in% c("start_codon","stop_codon") , ,drop=FALSE]
	#if (nrow(sub.feat) > 0) {
	#	sub.feat$seqname  <-  sub.feat$transcript_id
	#	sub.feat$strand  <-  "+"
	#	sub.feat$start  <-  sub.feat$cs.pos + pmin(abs(sub.feat$end - sub.feat$exon.orient.start) , abs(sub.feat$start - sub.feat$exon.orient.start) ) + 1
	#	sub.feat$end  <-  sub.feat$start + sub.feat$width - 1
	#	
	#	# Sometimes, a start/stop codon may be split over entries (bizarrely).
	#	# Now that they are oriented, just take the earliest one and extend it to three bases.
	#	
	#	sub.feat  <-  sub.feat[ order(sub.feat$transcript_id , sub.feat$start , sub.feat$end)  ,]
	#	sub.feat  <-  sub.feat[ !duplicated(sub.feat[,c("feature","transcript_id")])  ,]
	#	sub.feat$end  <-  sub.feat$start + 2
	#	
	#	myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
	#} # has start/stop
	#
	#
	#
	#### CDS
	#verb("\t\tCDS.\n")
	#
	#sub.feat  <-  myData.gtf[ myData.gtf$feature == "exon" , ,drop=FALSE]
	#if (nrow(sub.feat) > 0) {
	#	logi.dup  <-  duplicated(sub.feat[, c("transcript_id")])
	#	sub.feat  <-  sub.feat[ !logi.dup ,]
	#	
	#	sub.feat$feature  <-  "CDS"
	#	sub.feat$seqname  <-  sub.feat$transcript_id
	#	sub.feat$strand  <-  "+"
	#	
	#	sub.feat$start  <-  sub.feat$five_prime_utr.width + 1
	#	sub.feat$start[is.na(sub.feat$start)]  <-  1
	#	
	#	sub.feat$end  <-  sub.feat$transcript_length - sub.feat$three_prime_utr.width
	#	sub.feat$end[is.na(sub.feat$end)]  <-  sub.feat$transcript_length[is.na(sub.feat$end)]
	#	
	#	myData.tranGTF  <-  rbind(myData.tranGTF , sub.feat[, myData.colnames] )
	#} # has exon
	
	
	
	
	
	
	
	
	
	####################### sanity checks
	verb("\t\tsanity check.\n")
	
	tdf  <-  myData.tranGTF  
	tdf$transcript_length  <-  get_gtf_attribute_field( gtf = tdf , field = "transcript_length" , num = TRUE )
	tdf$five_prime_utr_length  <-  get_gtf_attribute_field( gtf = tdf , field = "five_prime_utr_length" , num = TRUE )
	tdf$three_prime_utr_length  <-  get_gtf_attribute_field( gtf = tdf , field = "three_prime_utr_length" , num = TRUE )
	tdf$inferred.CDS.length  <-  tdf$transcript_length - (tdf$five_prime_utr_length + tdf$three_prime_utr_length)
	tdf$CDS_length  <-  tdf$end - tdf$start + 1
	
	
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
	
	fname  <-  sprintf("%s.pdf", outfbase )
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

	
	##### save
	verb("\tsave.\n")
	
	fname  <-  sprintf("%s", outfbase)
	write.table( myData.tranGTF , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
	
} # make_eukaryotic_transcript_gtf












