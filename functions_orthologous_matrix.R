
# # # 			source("/home/map2085/riboprof/src/functions_orthologous_matrix.R")



##########  packages

library(data.table)
library(stringr)
library(methods)







######################################### verb

### use this function for stderr output
verb <- function(...) cat(sprintf(...), sep='', file=stderr())










######################################################
######################################################
######################################################
######################################################
#####################  load  #########################
######################################################
######################################################
######################################################
######################################################



load_orthologous_matrix_database  <-  function() {

	bdir  <-  "/home/map2085/riboprof/data/orthologous_matrix/release_Dec2017"

	oma  <-  list()

	### groups
	verb("\t\tgroups.\n")

	fname  <-  sprintf("%s/oma-groups.txt.gz", bdir)
	inl  <-  readLines( con = fname )
	inl  <-  inl[ !str_detect(inl , "^#") ] # remove header lines
	inl  <-  str_split( inl , "\t" )
	names(inl)  <-  sapply( inl , '[[' , 1)
	inl  <-  lapply( inl , FUN=function(x) x[-1] ) # remove number
	oma$groups  <-  inl
	
	
	### genbank mappings.
	verb("\t\tgenbank mappings.\n")

	fname  <-  sprintf("%s/oma-ncbi.txt.gz", bdir)
	indt  <-  fread(input=sprintf("zcat  %s", fname) , sep="\t" , header=F)
	colnames(indt)  <-  c("oma.id","genbank.id")
	oma$genbank  <-  copy(indt)


	### species info
	verb("\t\tspecies info.\n")

	fname  <-  sprintf("%s/oma-species.txt.gz", bdir)
	indt  <-  fread(input=sprintf("zcat  %s", fname) , sep="\t" , header=F)
	colnames(indt)  <-  c("oma.id","taxid","scientific.name","genome.source","version")
	oma$species  <-  copy(indt)

	return(oma)

} # load_orthologous_matrix_database






get_oma_species_id_from_oma_protein_id  <-  function(oma.id) {

	lets  <-  str_extract(oma.id , "^.....")

	return(lets)

} # get_oma_species_id_from_oma_protein_id









