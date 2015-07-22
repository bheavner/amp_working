# fix IDs in Mayo BB TCX countfiles
library(gdata) # to read .xlsx file -  install.packages("gdata")

genePath = "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_JulyRerun_gene_id_counts.txt"

geneCounts <- read.table(genePath, header = TRUE, check.names = FALSE)

keyPath  <- "/Users/bheavner/Desktop/TCXIDsKey.xlsx"

keyHash <- read.xls(keyPath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

# replace colnames(geneCounts) with corresponding keyHash$IlluminaSampleID values
colnames(geneCounts) <- keyHash$IlluminaSampleID[match(colnames(geneCounts), keyHash$Path_ID)]

# save and rezip the geneCounts table.
write.table(geneCounts, genePath, quote = FALSE)

test <- read.table(genePath, header = TRUE, check.names = FALSE)
colnames(test)

system(paste("gzip", genePath))

# then rename file to current conventions and upload it to synapse


transcriptPath = "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_JulyRerun_transcript_id_counts.txt"

transcriptCounts <- read.table(transcriptPath, header = TRUE, check.names = FALSE)

# replace colnames(transcriptCounts) with corresponding keyHash$IlluminaSampleID values
colnames(transcriptCounts) <- keyHash$IlluminaSampleID[match(colnames(transcriptCounts), keyHash$Path_ID)]

# save and rezip the geneCounts table.
write.table(transcriptCounts, genePath, quote = FALSE)

test <- read.table(transcriptPath, header = TRUE, check.names = FALSE)
colnames(test)

system(paste("gzip", transcriptPath))

# then rename file to current conventions and upload it to synapse
