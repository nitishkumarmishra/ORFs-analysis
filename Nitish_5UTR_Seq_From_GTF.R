## get UTR sequences
## https://www.r-bloggers.com/2011/03/the-easiest-way-to-get-utr-sequence/

require(biomaRt)
require(org.Mm.eg.db)
require(rtracklayer)


mart <- useEnsembl("ensembl", dataset="mmusculus_gene_ensembl")
eg <- getBM(c("external_gene_name","ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=mart)

genes <- eg$external_gene_name
genes <- genes[genes != ""]
genes <- genes[50:100]
## Lots of genes doesn't have known UTR sequences
## Best way is to read GTF file and select which have UTR sequences
my_file <- "C:/Users/nmishra/Desktop/Polysomes RNAseq EMT/jupyter codes/data/mouse/Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz"
show(my_file)
granges_gtf <- import(my_file)
granges_gtf <- as.data.frame(granges_gtf)
five_prime_utr_gtf <- granges_gtf[granges_gtf$type=="five_prime_utr",]
five_prime_utr_gtf <- as.data.frame(five_prime_utr_gtf)
#head(five_prime_utr_gtf, n = 3)
genes <- five_prime_utr_gtf$gene_name

seq = getSequence(id = genes, 
                  type = "external_gene_name", 
                  seqType = "5utr", 
                  mart = mart,
                  upstream = 100)
show(seq)


exportFASTA(seq, file= 'SEQ.fasta') 




#library(BSgenome.Mmusculus.UCSC.mm39)
## assign this to a short name for convenience
#mm39 <- BSgenome.Mmusculus.UCSC.mm39
## extract sequence between specific coordinates
#getSeq(x = mm39, names = "chrX", start = 100636100, end = 100636120)
