
# # # 			source("/home/map2085/riboprof/src/functions_rsem.R")




##########  packages
library(methods)
library(reshape2)









######################################################################
######################################################################
######################################################################
######################################################################
###########################  filenames  ##############################
######################################################################
######################################################################
######################################################################
######################################################################






#' Get filename for rsem results file(s)
#'
#' This function will return the filename of an rsem results file for a given sample.
#' This is vectorized if 'project' is a data frame.
#'
#' @param  project	character vector giving project name. must be length 1.
#' @param  individual	character vector giving individual. must be length 1.
#' @param  runid	character vector giving runid. must be length 1.
#' @param  type		character vector giving type of analysis: rrna-pop,mrna, etc.
#' @param  level	character vector: 'genes' or 'isoforms'. This is the level of RSEM output. usually you want 'genes'.
#'
#' @return	character vector of length 1 giving filename
rsem_analysis_type_to_filename_prefix  <-  function( project , individual = NULL, runid = NULL , type , level ,  verbose = FALSE ) {

	if (verbose) { verb("\trsem_analysis_type_to_filename_prefix.\n") }

	#### check
	stopifnot((is.null(individual)  &&  is.null(runid)  &&  is.data.frame(project))  ||  (!is.null(individual)  &&  !is.null(runid)  &&  is.character(project)  &&  is.character(individual)  && is.character(runid) ) )

	#### data frame to strings
	if (is.data.frame(project)) {
		stopifnot(is.null(individual)  && is.null(runid))
		runid  <-  project$runid
		individual  <-  project$individual
		project  <-  project$project 
	} # data frame

	if (type == "rrna-pop") {
		bdir  <-  "/home/map2085/riboprof/out/rsem-rnaseq-rrna-population"
		btag  <-  "rsem-rnaseq-rrna-population"
	} else {
		verb("\n\n\nERROR!  unrecognized type=[%s]!!!\n", type)
		stop()
	} # type


	# vectorized
	the.dir  <-  paste( bdir , project , individual , runid , sep = "/" )
	the.tag  <-  paste( btag , project , individual , runid , level , sep = "." ) 

	fname  <-  paste( the.dir , the.tag , sep = "/" )
	fname  <-  paste( fname , "results.gz" , sep = "." )
	
	if (verbose) { show(the.dir) ; show(the.tag) ; show(fname) } 

	return(fname)

} # rsem_analysis_type_to_filename_prefix









get_rsem_transcriptome_bam_filenames  <-  function( project , type ) {

        bdir  <-  "/home/map2085/riboprof/out/rsem-rnaseq-annotated-genome"

        stopifnot(is.data.frame(project))

        fnames  <-  c()
        for (rx  in  1:nrow(project)) {
                projx  <-  project$project[rx]
                indivx  <-  project$individual[rx]
                runx  <-  project$runid[rx]

                if (type == "mrna") {
                        fname  <-  sprintf("%s/%s/%s/%s/rsemRnaseq-annot.%s.%s.%s.transcript.bam" , bdir , projx , indivx , runx , projx , indivx , runx )
                        fnames  <-  c(fnames , fname)
                } else {
                        verb("\n\n\nERROR!  unrecognized type=[%s]!!!\n", type)
                        stop()
                } # type
        } # rx

        return(fnames)

} # get_rsem_transcriptome_bam_filenames



























######################################################################
######################################################################
######################################################################
######################################################################
###############################   load  ##############################
######################################################################
######################################################################
######################################################################
######################################################################






#' Load RSEM results
#'
#' This function loads RSEM results files and concatenates them, identified by runid
#'
#' @param  project      character vector giving project name or data frame of project.
#' @param  level	character vector: genes,isoforms. must be length 1
#' @param  type		character vector giving type of analysis: rrna-pop,mrna, etc.
#' @param  matrix       if TRUE, then will cast to a matrix. rownames = gene/transcript names, colnames = runids.
#'			if FALSE, then will be in long format as a data frame.
#'
#' @return      matrix or data frame, depending on "matrix" argument.
#'
#' @export
load_rsem_counts  <-  function( project , level , type , matrix , verbose = FALSE ) {

	stopifnot(level  %in%  c("genes","isoforms"))
	stopifnot(length(level) == 1)

	project  <-  process_project_input( project = project )

	rsem.fnames  <-  rsem_analysis_type_to_filename_prefix( project = project , type = type , level = level , verbose = verbose )

	gendf  <-  data.frame()
	for (rx  in  1:nrow(project)) {
		projx  <-  project$project[rx]
		indivx  <-  project$individual[rx]
		runx  <-  project$runid[rx]

		if (verbose) { verb("\t\t%s\n", runx) }

		fname  <-  rsem.fnames[rx]

		indf  <-  read.table( file = fname , header = TRUE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "" )
		indf$runid  <-  runx
		indf$project  <-  projx
		indf$individual  <-  indivx

		gendf  <-  rbind( gendf , indf )
	} # rx

	gendf$expected_count  <-  as.numeric(gendf$expected_count)


	### matrix
	if (matrix) {
		if (level == "genes") {
			cmat  <-  acast( data = gendf , formula = gene_id ~ runid , value.var = "expected_count" , fill = 0 )	
		} else if (level == "isoforms") {
			cmat  <-  acast( data = gendf , formula = transcript_id ~ runid , value.var = "expected_count" , fill = 0 )	
		} # level

		stopifnot(all(colnames(cmat) %in% project$runid))
		stopifnot(all(project$runid %in% colnames(cmat)))

		return(cmat)
	} else {
		return(gendf)
	}

} # load_rrna_population_expression















#' Load combined rsem counts for a project
#'
#' Load combined rsem counts as created by 'combine_rsem_transcriptome_for_project.sh'.
#'
load_combined_rsem_counts_for_project  <-  function( project , type , verbose = FALSE ) {

	stopifnot(type %in% c("gene","isoform"))
	bdir  <-  "/home/map2085/riboprof/out/combined-rsem-transcriptome"

	if (is.vector(project)  &&  is.character(project)) {
		project  <-  data.table(project = project)
	} # project

	stopifnot(is.data.frame(project))
	stopifnot("project" %in% colnames(project))

	alldt  <-  data.table()
	for (projx  in  unique(project$project)) {
		if (verbose) { verb("\t%s\n", projx) }
			
		fdir  <-  sprintf("%s/%s", bdir, projx)
		ftag  <-  sprintf("combined-rsem-transcriptome.%s.%s.gz" , projx , type)
		fname  <-  sprintf("%s/%s", fdir, ftag)

		indt  <-  fread( input = sprintf("zcat  %s", fname) , header = TRUE , sep="\t" , stringsAsFactors = FALSE , quote = "" )
		alldt  <-  rbind( alldt, indt)
	} # projx

	return(alldt)

} # load_combined_rsem_counts_for_project

















######################################################################
######################################################################
######################################################################
######################################################################
############################   normalize  ############################
######################################################################
######################################################################
######################################################################
######################################################################



#' Normalize gene counts by groups of samples
#'
#' For each group of samples, the gene counts will be normalized by the specified method
#' 
#' @param  counts	data.table with columns: count,gene_id,sample
#' @param  group	character vector of length = nrows(counts) giving the group for which samples will be normalized.
#' @param  norm		character string indicating the normalization method.
#'
#' @return	DGEList from limma
#'
#' @export
normalize_counts_by_group  <-  function( counts , group , norm , ... ) {

	sgdt  <-  unique(data.table( sample = counts$sample , group = group ))
	stopifnot(all(table(sgdt$sample) == 1))
	remove(sgdt)

	counts  <-  as.data.table(counts)
	counts$group  <-  group

	spl  <-  split( x = counts , by = "group" , keep.by = TRUE , use.names=TRUE)
	dgel  <-  list()
	for (splx  in  spl) {
		gx  <-  splx$group[1]
	
		cmat  <-  acast( splx , gene_id ~ sample , value.var = "count" , fill = 0)
		dge  <-  DGEList( counts = cmat , genes = rownames(cmat))
		dge  <-  calcNormFactors( dge , method = norm , ... )
		nl  <-  list(dge)
		names(nl)  <-  gx
		dgel  <-  c( dgel , nl )
	} # splx

	return(dgel)

} # normalize_counts_by_group


















