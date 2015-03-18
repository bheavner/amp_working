# Do we need to remove genes with no readcounts prior to normalizing? Look at Tau data to explore the question.

require("synapseClient")
require("R.utils")
require("edgeR")

synapseLogin()

mergedCountFile <- ('syn3160706')

originalCountFile <- synGet(mergedCountFile)

# unzip file and load for processing
localFilePath <- getFileLocation(originalCountFile)

if(!file.exists(substr(localFilePath, 1, nchar(localFilePath) - 3))) {
  gunzip(localFilePath)
}

localFilePath <- sub('.gz', '', localFilePath) #trim the .gz suffix

rawCounts <- read.table(localFilePath, header = TRUE)

# does this have headers, etc. that I need?
str(rawCounts)
# 'data.frame':  60 obs. of  39179 variables:
# $ ENSMUSG00000000001: int  1499 1569 1364 1393 1530 1470 1569 1368 1649 1391 ...
# $ ENSMUSG00000000003: int  0 0 0 0 0 0 0 0 0 0 ...
# $ ENSMUSG00000000028: int  20 15 29 21 26 19 25 22 21 26 ...
# $ ENSMUSG00000000031: int  2 17 12 11 12 10 10 5 10 6 ...

row.names(rawCounts)
# row.names looks like this:
# "4088"     "4095"     "4808"     "NA00-136" "NA00-175"
# I think that's the sample names

colnames(rawCounts)
# looks like this:
# "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457"
# I think that's the gene IDs

# begin processing - first transpose b/c James did it differently than DGEList expects - I'll want to change James' code to change htis later.
transposedCounts <- t(rawCounts)

str(transposedCounts)
#int [1:39179, 1:60] 1499 0 20 2 70 5 1825 1374 2073 1460 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:39179] "ENSMUSG00000000001" "ENSMUSG00000000003" "ENSMUSG00000000028" "ENSMUSG00000000031" ...
#..$ : chr [1:60] "LP62_4" "R-g2_3" "R697_7" "R706_1" ...

# numerics changed to character strings?
typeof(transposedCounts)
# [1] "integer" - no, okay at the moment.

row.names(transposedCounts)
# row.names looks like this:
#"ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457"
# I think that's the gene IDs

colnames(transposedCounts)
# looks like this:
# "4088"     "4095"     "4808"     "NA00-136" "NA00-175"
# I think that's the sample names

# how many rows have zero entries in them?
dim(transposedCounts[(apply(transposedCounts, 1, function(y) any(y == 0))),])
# 21543!

# how many have _all_ zero entries in them?
dim(transposedCounts) - dim(transposedCounts[rowSums(transposedCounts[, -1])>0, ])
# 10252 do! That seems like a lot!

#Or...
sum(as.numeric(rowSums(transposedCounts < 1) >= 60))
#[1] 10231 (I'm not sure why it's not 10252)

# I think I do want to remove these genes prior to normalization.