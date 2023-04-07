require(seqinr)
require(ape)

alignment <- read.alignment(file = "alignment.aln", format = "clustal")
alignment

aln_clean <- cleanAlignment(alignment, minpcnongap = 50, minpcid = 50)
m <- dist.alignment(aln_clean)


n <- nj(m)
# boot.phylo(n, n)

h <- hclust(m, method = "average")

plot.phylo(n, type = "u")

plot(h)

unrootedNJtree(aln_clean, type = "protein")

rootedNJtree(aln_clean, type = "protein")

# -------------------------------------------------------------------------------
aln <- read.alignment(file = "fox_protein.aln", format = "clustal")
aln

dist_mat <- dist.alignment(aln)
dist_mat

tree <- rootedNJtree(aln, theoutgroup = "tr|Q9VT99|Q9VT99_DROME", type = "protein")
write.tree(tree, file = "tree.tre")
