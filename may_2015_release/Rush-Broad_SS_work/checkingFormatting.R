# covariates #
#look ok

# gene counts #
path <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_gene_id_counts.txt"
newpath  <- "/Users/bheavner/Desktop/AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_GeneCounts.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

write.table(geneCounts, newpath, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

system(paste("gzip", newpath, sep = " "))


### normalized genecounts ##

path <- "/Users/bheavner/Desktop/AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_geneCounts_normalized.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

#fix colnames to agree with covariates
colnames(geneCounts) <- sub('X', '', colnames(geneCounts))

write.table(format(geneCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))

### transcript counts ##

path2 <- "/Users/bheavner/Desktop/AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_TranscriptCounts.txt"

transcriptCounts <- read.table(path2, header = TRUE, check.names = FALSE)

write.table(transcriptCounts, path2, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

system(paste("gzip", path2, sep = " "))

## normalized transcript counts ##

path <- "/Users/bheavner/Desktop/AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_TranscriptCounts_Normalized.txt"

transcriptCounts <- read.table(path, header = TRUE, check.names = FALSE)

#fix colnames to agree with covariates
colnames(transcriptCounts) <- sub('X', '', colnames(transcriptCounts))

write.table(format(transcriptCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))
