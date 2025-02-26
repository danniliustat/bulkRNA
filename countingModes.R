genePlot <- function(){
my.gr <- GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=1, end=1000), strand=Rle("*") )
my.seqinfo <- Seqinfo(seqnames="chr1", seqlengths=1000, isCircular=FALSE, genome="simulated")
seqinfo(my.gr) <- my.seqinfo
gene1 <- GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=201, end=350), strand=Rle("+"), feature="protein_coding", gene="gene1", exon="exon1", transcript="transcript1.1", symbol="gene1", seqinfo=my.seqinfo)
gene2 <- GRanges(seqnames=Rle(rep("chr1", 2)), ranges=IRanges(start=c(301, 601), end=c(350, 700)), strand=Rle(rep("+", 2)), feature=rep("protein_coding", 2), gene=rep("gene2", 2), exon=c("exon2.1", "exon2.2"), transcript=rep("transcript2.1", 2), symbol="gene2", seqinfo=my.seqinfo)
gene3 <-  GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=451, end=550), strand=Rle("+"), feature="protein_coding", gene="gene3", exon="exon3.1", transcript="transcript3.1", symbol="gene3", seqinfo=my.seqinfo)
gene4 <- GRanges(seqnames=Rle(rep("chr1", 2)), ranges=IRanges(start=c(101, 301), end=c(200, 400)), strand=Rle(rep("-", 2)), feature=rep("protein_coding", 2), gene=rep("gene4", 2), exon=c("exon4.1", "exon4.2"), transcript=rep("transcript4.1", 2), symbol="gene4", seqinfo=my.seqinfo)
my.genes <- GRangesList(gene1=gene1, gene2=gene2, gene3=gene3, gene4=gene4)
my.genes.df <- as.data.frame(my.genes)
colnames(my.genes.df)[3] <- "chromosome"
gene.track <-GeneRegionTrack(my.genes.df, fill="green4", arrowHeadWidth=50, arrowHeadMaxWidth=50, col="grey75", shape="arrow", stacking="squish", genome="simulated", names="Demo Genes", grid=TRUE, lty.grid=2, geneSymbols=TRUE, cex.group=1)

gene1.reads <- GRanges(seqnames=Rle(rep("chr1", 20)), ranges=IRanges(start=sample(200:325, 20, replace=T), width=25), strand=Rle(rep("+", 20)), seqinfo=my.seqinfo)
gene2.reads <- GRanges(seqnames=Rle(rep("chr1", 60)), ranges=IRanges(start=sample(c(300:325, 600:675), 60, replace=T), width=25), strand=Rle(rep("+", 60)), seqinfo=my.seqinfo)
gene3.reads <- GRanges(seqnames=Rle(rep("chr1", 10)), ranges=IRanges(start=sample(450:525, 10, replace=T), width=25), strand=Rle(rep("+", 10)), seqinfo=my.seqinfo)
gene4.reads <- GRanges(seqnames=Rle(rep("chr1", 20)), ranges=IRanges(start=sample(c(100:175, 300:375), 20, replace=T), width=25), strand=Rle(rep("-", 20)), seqinfo=my.seqinfo)
all.reads <- c(gene1.reads, gene2.reads, gene3.reads, gene4.reads)

read.track <- AnnotationTrack(all.reads, fill=c("red", "blue")[as.integer(strand(all.reads))], stacking="squish", genome="simulated", name="Demo Reads", lty.grid=2)

ax.track <- GenomeAxisTrack(GRanges(seqnames="chr1", ranges=IRanges(1, 1000)), genome="simulated")

plotTracks(list(gene.track, read.track, ax.track), from=1, to=700, grid=1, sizes=c(4, 4, 1))
count.mat <- matrix(0, nrow=4, ncol=5)
colnames(count.mat) <- c("actual", "countOverlaps", "Union", "IntersectionStrict", "IntersectionNotEmpty")
rownames(count.mat) <- c("gene1", "gene2", "gene3", "gene4")
count.mat[, 1] <- c(20, 60, 10, 20)
count.mat[, 2] <- countOverlaps(my.genes, all.reads)
count.mat[, 3] <- assays(summarizeOverlaps(my.genes, all.reads, "Union", ignore.strand=FALSE))$counts
count.mat[, 4] <- assays(summarizeOverlaps(my.genes, all.reads, "IntersectionStrict", ignore.strand=FALSE))$counts
count.mat[, 5] <- assays(summarizeOverlaps(my.genes, all.reads, "IntersectionNotEmpty", ignore.strand=FALSE))$counts
count.mat
}
