# working to fix annotation, sigfigs, etc.
# ad_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160436')

require("synapseClient")
require("R.utils")
require("edgeR")

synapseLogin()

mergedCountFile <- ('syn3160436')

message("Normalizing ", mergedCountFile)
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
# 96 obs of 64253 variables - looks like:
#$ ENSG00000000003: int  530 2395 448 815 310 222 552 1125 1107 928 ...
#$ ENSG00000000005: int  5 4 3 5 4 3 4 5 4 0 ...
#$ ENSG00000000419: int  776 488 713 496 486 526 1081 934 559 620 ...

row.names(rawCounts)
# row.names looks like this:
# "4088"     "4095"     "4808"     "NA00-136" "NA00-175"
# I think that's the sample names

colnames(rawCounts)
# looks like this:
# "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457"
# I think that's the gene IDs

# begin processing - first transpose b/c James did it differently than DGEList expects
transposedCounts <- t(rawCounts)

str(transposedCounts)
#int [1:64253, 1:96] 530 5 776 381 277 600 387 785 2333 1370 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:64253] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" ...
#..$ : chr [1:96] "4088" "4095" "4808" "NA00-136" ...

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

# make smaller toy example to work with - all samples, 5 genes
smaller <- transposedCounts[1:5,]

expr <- DGEList(smaller, group = rep(1, ncol(smaller)))

# calculate normalization factors
normFactors <- calcNormFactors(expr, method = ("TMM"))

# use normaliztion factors to calculate cpm -
# per https://www.biostars.org/p/84087/, that's calculated as
# count / (library size * normalization factor))

normalizedCpm <- cpm(normFactors)

# trim to 7 sig figs
normalizedCpm <- format(normalizedCpm, digits = 7)

write.table(normalizedCpm, "testfile.txt", quote = FALSE, sep = "\t", row.names = TRUE)

# resulting .txt looks right to me - header is sample IDs,
