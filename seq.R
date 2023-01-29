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

q2$req[[4]]


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
      query(queryname, thequery)
      seq <- getSequence(query2$req[1])
      closebank()
      return (seq)
    }
    closebank()
  }
  print(paste("ERROR: accession", accession, "was not found"))
}


x <- getncbiseq("NC_001477")
