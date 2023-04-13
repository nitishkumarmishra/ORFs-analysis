## Codon count table
library("coRdon")
library('Biostrings')


dna <- DNAStringSet(c("AAAACGAAGTGTACTGTAATTTGCACAGTACTTAAATGTAAGTAAAAA",
                      "ACGTCCGTACTGATCGATTCCGTGATT"))
cT <- codonTable(dna)
mat <- codonCounts(cT)



###### Codon frequency
library(codondiffR)


fastaFile <- system.file(
  "extdata", "example_viruses.fna", package="codondiffR"
)
virusSet <- readSeq(file = fastaFile)

class(virusSet)



## Create a codonFreq object from the DNAStringSet.
virusCF <- codonFreq(virusSet)

class(virusCF)
mat1 <- virusCF@freq
rownames(mat1) <- virusCF@seqID




##### Codon counts from fasta file
fastaFile <- system.file(
  "extdata", "example_viruses.fna", package="codondiffR"
)
virusSet <- readSeq(file = fastaFile)


cT <- codonTable(virusSet)
mat <- codonCounts(cT)
rownames(mat) <- names(virusSet)