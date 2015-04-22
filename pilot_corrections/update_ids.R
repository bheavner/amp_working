# fix IDs in PSP countfiles
library(gdata) # to read .xlsx file -  install.packages("gdata")

genePath = "/Users/bheavner/Desktop/psp/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_gene_id_counts_UpdatedID.txt"

geneCounts <- read.table(genePath, header = TRUE)

# colnames need to change to match the Path_ID
# remove leading X
colnames(geneCounts)  <- sub('X', '', colnames(geneCounts))
#replace . with -
colnames(geneCounts)  <- sub('[.:punct:]', '-', colnames(geneCounts))

keyPath  <- "/Users/bheavner/Desktop/Data\ errors/ToSend/PSP_MA_040915_IDsKey.xlsx"

keyHash <- read.xls(keyPath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

# replace colnames(geneCounts) with corresponding keyHash$IlluminaSampleID values
colnames(geneCounts) <- keyHash$IlluminaSampleID[match(colnames(geneCounts), keyHash$Path_ID)]

# save and rezip the geneCounts table.
write.table(geneCounts, genePath) #, quote = FALSE, sep = " ", row.names = TRUE)

#test <- read.table(genePath, header = TRUE)

system(paste("gzip", genePath))

# then upload it to synapse


transcriptPath = "/Users/bheavner/Desktop/psp/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_transcript_id_counts_UpdatedID.txt"

transcriptCounts <- read.table(transcriptPath, header = TRUE)

# colnames need to change to match the Path_ID
# remove leading X
colnames(transcriptCounts)  <- sub('X', '', colnames(transcriptCounts))
#replace . with -
colnames(transcriptCounts)  <- sub('[.:punct:]', '-', colnames(transcriptCounts))


# replace colnames(transcriptCounts) with corresponding keyHash$IlluminaSampleID values
colnames(transcriptCounts) <- keyHash$IlluminaSampleID[match(colnames(transcriptCounts), keyHash$Path_ID)]

# save and rezip the geneCounts table.
write.table(transcriptCounts, genePath, quote = FALSE, sep = " ", row.names = TRUE)

#test <- read.table(genePath, header = TRUE)

system(paste("gzip", transcriptPath))
