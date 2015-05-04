
# gene counts #
path <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_Alzheimer_gene_id_counts.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

# re-id by stepping through update_ids code


write.table(geneCounts, path, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))

### normalized genecounts ##

path <- "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_Normalized.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

write.table(format(geneCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))

### transcript counts ##

path <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_Alzheimer_TranscriptCounts_UpdatedID.txt"

transcriptCounts <- read.table(path, header = TRUE, check.names = FALSE)

# re-id by stepping through update_ids code

write.table(transcriptCounts, path, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))

## normalized transcript counts ##

### transcript counts ##

path <- "/Users/bheavner/Desktop/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts_Normalized.txt"

transcriptCounts <- read.table(path, header = TRUE, check.names = FALSE)

# re-id by stepping through update_ids code

write.table(format(transcriptCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))


### PSP:::

path <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_GeneCounts_UpdatedID.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

write.table(geneCounts, path, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

#test <- read.table(path, header = TRUE, check.names = FALSE, quote = " ")

system(paste("gzip", path, sep = " "))


## transcript counts

path <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_Transcript_Counts_UpdatedID.txt"

transcriptCounts <- read.table(path, header = TRUE, check.names = FALSE)

# fix IDs

write.table(transcriptCounts, path, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

#test <- read.table(path, header = TRUE, check.names = FALSE)

system(paste("gzip", path, sep = " "))



path <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_TranscriptCounts_Normalized_UpdatedID.txt"

transcriptCounts <- read.table(path, header = TRUE, check.names = FALSE)

# re-id by stepping through update_ids code

write.table(format(transcriptCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

#test <- read.table(path, header = TRUE, check.names = FALSE)

system(paste("gzip", path, sep = " "))
