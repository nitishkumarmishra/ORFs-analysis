
#	  source("/home/map2085/riboprof/src/functions_orthologous_gene_maps.R")





########################################## libraries

library(stringr)
library(reshape2)

source("/home/map2085/riboprof/src/functions_entrez.R")




############################################ stringsAsFactors

options(stringsAsFactors = FALSE)




######################################### verb

### use this function for stderr output
verb <- function(...) cat(sprintf(...), sep='', file=stderr())












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

	orth.dir  <-  "/home/map2085/riboprof/data/mgi"
	orth.name  <-  "HOM_MouseHumanSequence.rpt.gz"
	
	fname  <-  sprintf("%s/%s", orth.dir , orth.name )

	indf  <-  read.table( file = fname , header = TRUE , sep = "\t" , row.names = NULL , stringsAsFactors = FALSE , quote = "" , fill = TRUE )

	if (!is.null(from)  ||  !is.null(to))  {
		indf  <-  restrict_orthologous_genes_by_species( ortho = indf , from = from , to = to )
	} # from, to

	return(indf)

} # load_orthologous_gene_maps_MouseHuman













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
	if (from == "mm10") { from  <-  "mouse" }
	if (to == "mm10") { to  <-  "mouse" }

	ortho$species  <-  str_extract( ortho$Common.Organism.Name , "^[^,]+" )
	ortho  <-  subset( ortho , species %in% c(from , to) )

	return(ortho)

} # restrict_orthologous_genes_by_species


















#' Create an orthologous gene map table
#'
#' This is usually an internal function, but can be used externally.
#' From long format to wide format, with duplicate entry (rows) if necessary.
#' 
#' @param	ortho		data frame, ortho table
#' @param	from, to	character vector, selected species
#'
#' @return	wide format data frame, each row is a HomolGene.ID from the ortho table
#' @export
make_orthologous_gene_map_table  <-  function( ortho = NULL , from , to , verbose = FALSE ) {

	if (is.null(ortho)) {
		ortho  <-  load_orthologous_gene_maps_MouseHuman(from = from , to = to)
	} # ortho

	if (!is.null(from)  ||  !is.null(to))  {
		ortho  <-  restrict_orthologous_genes_by_species( ortho = ortho , from = from , to = to )
	} # from, to

	# get all species
	uspec  <-  unique(ortho$species)

	sdf  <-  subset( ortho , select = c(HomoloGene.ID , species , EntrezGene.ID) )

	# split by species so can do a merge.
	# then, must rename columns so that the symbols become the species names
	speclist  <-  split( sdf , sdf$species )   # species list
	elist  <-  list()
	for (spx  in  names(speclist)) {

		### get entrez table
		entrez  <-  get_entrez_to_alias_table( species = spx )
		# columns = entrez_id , gene_name
		
		gdf  <-  unique(speclist[[spx]][, c("HomoloGene.ID","EntrezGene.ID") ,drop=FALSE])
		gdf  <-  merge( entrez , gdf , by.x = "entrez_id" , by.y = "EntrezGene.ID" )
		gdf  <-  gdf[,c("HomoloGene.ID","gene_name"),drop=FALSE]
		gdf  <-  unique(gdf)
		colnames(gdf)[colnames(gdf)=="gene_name"]  <-  spx
		
		elist[[spx]]  <-  gdf
	} # spx

	# Reduce via merge
	fmerg  <-  function(a,b) { merge(a , b , by = "HomoloGene.ID" , all = TRUE , fill = NA ) }
	mdf  <-  Reduce( f = fmerg , x = elist )
	mdf  <-  mdf[ order(mdf$HomoloGene.ID) , ,drop=FALSE]

	return(mdf)

} # make_orthologous_gene_map_table











#' Map gene names from one species to another within a given data frame.
#'
#' Given a data frame with a "gene_name" column, this function will convert the gene name to that of another species.
#'
#' @param  df		a data frame with column "gene_name".
#' @param  from		character string giving the species being mapped from.
#'			Thus, the gene_name column of "df" contains gene names fro the "from" species.
#' @param  to		character string giving the species being mapped to.
#' @param  universe	universe of gene names. This is because there are often many synonymous gene names ("Symbols") for the same gene.
#' @param  table	output of "make_orthologous_gene_map_table".  optional. Automatically lodaded if NULL.
#' @param  keep.NA	if FALSE, then gene names in "df" that do not have an orthologoue in the species being mapped to will be removed.
#'
#' @return  data.frame identical to input "df" but with "gene_name" column giving the gene name of the orthologous gene in the "to" species.
#'		Note that multiple entries may be created if a gene has multiple synonyms or multiple mappings exist.
#'		Old gene names are stored in the column "gene_name.<from>"
#'
#' @export
map_to_orthologous_gene_name  <-  function( df , from , to , universe = NULL , table = NULL , verbose = FALSE , keep.NA = FALSE ) {

	if (from == to) {
		return(df)
	} # same

	if (is.null(table)) {
		table  <-  make_orthologous_gene_map_table( from = from , to = to , verbose = verbose )
	} # table
	
	stopifnot(from %in% colnames(table))
	stopifnot(to %in% colnames(table))


	nameold  <-  sprintf("gene_name.%s" , from )

	df[[nameold]]  <-  df$gene_name
	df$gene_name  <-  NULL

	colnames(table)[colnames(table)==to]  <-  "gene_name"
	table  <-  unique(table[,c(from,"gene_name"),drop=FALSE])

	df  <-  merge( df , table , by.x = nameold , by.y = from , all.x = TRUE )

	if (!keep.NA) {
		df  <-  df[ !is.na(df$gene_name) ,,drop=FALSE]
	} # df

	if (!is.null(universe)) {
		df  <-  df[ df$gene_name %in% universe | is.na(df$gene_name) ,,drop=FALSE]
	} # universe

	return(df)

} # map_to_orthologous_gene_name


























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
	if (from == "mm10") { from <- "mouse" }
	if (to == "mm10") { to <- "mouse" }

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











