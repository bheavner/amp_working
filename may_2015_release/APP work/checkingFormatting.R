# gene counts #
path <- "/Users/bheavner/Desktop/AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

write.table(geneCounts, path, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))

### normalized genecounts ##

path <- "/Users/bheavner/Desktop/AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts_Normalized.txt"

geneCounts <- read.table(path, header = TRUE, check.names = FALSE)

#fix colnames to agree with covariates
colnames(geneCounts) <- sub('X', '', colnames(geneCounts))

colnames(geneCounts) <- sub('[/.]', '-', colnames(geneCounts))

write.table(format(geneCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))

### transcript counts ##

path2 <- "/Users/bheavner/Desktop/AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts.txt"

transcriptCounts <- read.table(path2, header = TRUE, check.names = FALSE)

write.table(transcriptCounts, path2, quote = FALSE, sep = " ", row.names = TRUE, col.names = NA)

system(paste("gzip", path2, sep = " "))

## normalized transcript counts ##

### transcript counts ##

path <- "/Users/bheavner/Desktop/AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts_Normalized.txt"

transcriptCounts <- read.table(path, header = TRUE, check.names = FALSE)

#fix colnames to agree with covariates
colnames(transcriptCounts) <- sub('X', '', colnames(transcriptCounts))

colnames(transcriptCounts) <- sub('[/.]', '-', colnames(transcriptCounts))

write.table(format(transcriptCounts, scientific = FALSE, digits = 5), 
            path, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

system(paste("gzip", path, sep = " "))
