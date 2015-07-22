# hongdong noticed quotes on row names in files. I don't want those.

# syn4650257 - counts by gene ID was syn3818244
# syn4650265 - normalized counts by gene ID was syn3818429
# syn4650258 - counts by transcript ID was syn3818611
# syn4650430 - normalized counts by transcript ID was syn3818905

genePath = "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_JulyRerun.txt"

geneCount <- read.table(genePath, header = TRUE, check.names = FALSE)
write.table(geneCount, genePath, quote = FALSE)
system(paste("gzip", genePath))

# upload to synapse
## clear environment here for safety
rm(list = ls())

genePath = "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_Normalized_JulyRerun.txt"

geneCount <- read.table(genePath, header = TRUE, check.names = FALSE)
write.table(geneCount, genePath, quote = FALSE)
system(paste("gzip", genePath))

## clear environment here for safety
rm(list = ls())

genePath = "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts_JulyRerun.txt"

geneCount <- read.table(genePath, header = TRUE, check.names = FALSE)
write.table(geneCount, genePath, quote = FALSE)
system(paste("gzip", genePath))

## clear environment here for safety
rm(list = ls())

genePath = "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts_Normalized_JulyRerun.txt"

geneCount <- read.table(genePath, header = TRUE, check.names = FALSE)
write.table(geneCount, genePath, quote = FALSE)
system(paste("gzip", genePath))

## clear environment here for safety
rm(list = ls())

# then reannotate/provenance