library(seqinr)

read_seq <- function(x) {
  obj <- read.fasta(file = x)
  
  getSequence(obj)
  
  obj[[1]][1:10]
  
  seq <- getSequence(obj)
  
  seq[[1]][1:10]
  
  seq_vec <- unlist(seq)
  
  return (seq_vec)  
}

den1 <- read_seq("den1.fasta")
den2 <- read_seq("den2.fasta")
den3 <- read_seq("den3.fasta")

min_length <- min(length(den1), length(den2), length(den3))

den1 <- den1[1:min_length]
den2 <- den2[1:min_length]
den3 <- den3[1:min_length]


den1[den1 != den2]


# Using the database directly
choosebank("refseqViruses")

q <- query("query1", "AC=NC_001477")
attributes(q)

seq1 <- getSequence(q$req[[1]])

seq1

x1 <- getName(q$req[[1]])  
x2 <- getAnnot(q$req[[1]])

x1
x2

seq1[1:50]

q2 <- query("query2", "SP=Dengue Virus AND M=RNA")
q2

attributes(q2)

dengueSeq <- getSequence(q2$req[[1]])


closebank()


getncbiseq <- function(accession) {
  require("seqinr")
  dbs <- c("genbank", "refseq", "refseqViruses", "bacterial")
  numdbs <- length(dbs)
  
  for (i in 1: numdbs) {
    db <- dbs[i]
    choosebank(db)
    
    resquery <- try(query(".tmpquery", paste("AC=",accession)), silent = TRUE)
    
    if (!(inherits(resquery, "try-error"))) {
      queryname <- "query2"
      thequery <- paste("AC=", accession, sep="")
      query2 <- query(queryname, thequery)
      seq <- getSequence(query2$req[1])
      closebank()
      return (seq)
    }
    closebank()
  }
  print(paste("ERROR: accession", accession, "was not found"))
}


x <- getncbiseq("NC_001477")



# Class work - 3rd Feb

d1 <- list(getSequence(q2$req[[1]]))
d2 <- list(getSequence(q2$req[[2]]))
d3 <- list(getSequence(q2$req[[3]]))
d4 <- list(getSequence(q2$req[[4]]))

dengue <- getSequence(q2)

write.fasta(dengue, c("M1", "M2", "M3", "M4"), file.out = "my_output.fasta")

temp <- d1[[1]]
length(temp)

bases <- table(d1)

barplot(bases)


slidingwindowplot <- function(windowsize, inputseq) {
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkGCs <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    chunkGCs[i] <- chunkGC
  }
  plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}


slidingwindowplot(2000, d1[[1]])


myTable = count(d1[[1]], 2)
myTable

class(myTable)
names(myTable)

pie(myTable, labels = names(myTable))


# repr <- function (seq, pattern) {
#   count(seq, 2)["ta"] / (count(seq, 1)["t"] * count(seq, 1)["a"])  
# }

count(d1[[1]], 2)["ta"] / (count(d1[[1]], 1)["t"] * count(d1[[1]], 1)["a"])

count(d1[[1]], 2)["cg"] / (count(d1[[1]], 1)["c"] * count(d1[[1]], 1)["g"])


l <- length(d1[[1]])
(count(d1[[1]], 2)["cg"] / l) / ((count(d1[[1]], 1)["c"] / l) * (count(d1[[1]], 1)["g"] / l))

s <- sample(c("a", "c", "g", "t"), 1000, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
s

prob <- table(s) / sum(table(s))
prob
