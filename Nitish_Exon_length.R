library(GenomicFeatures)
library(dplyr)
txdb <- makeTxDbFromGFF("Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then :: https://www.biostars.org/p/83901/
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
gene.exon.counts <- lengths(exons.list.per.gene)
gene.info <- cbind(exonic.gene.sizes, gene.exon.counts)
colnames(gene.info) <- c("exon_length", "exon_number")


gtf <- rtracklayer::import("Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz")
gencode.vM28.gtf <- as.data.frame(gtf)
gencode.vM28.gtf.selected <- gencode.vM28.gtf %>%
  filter(type=="gene") %>%
  rename(Chr= seqnames, Type= type, GeneType=gene_biotype, Symbol= gene_name, Status= gene_version, Havana=gene_source, Length=width, Start = start, End=end, Strand=strand) %>%
  dplyr::select(gene_id, Symbol,  Chr, Start, End, Strand,  GeneType,  Status, Havana, Length) %>%
  filter(!grepl("chrM", Chr))


gene.info <- merge(gencode.vM28.gtf.selected, gene.info, by.x="gene_id", by.y= "row.names")
head(gene.info, n=3)

