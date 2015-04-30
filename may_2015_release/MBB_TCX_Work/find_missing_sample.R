localFilePath  <- '/tmp/AMP-AD_MSBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_gene_id_counts.txt'
sampleListPath  <- '/tmp/TCXSamples.csv'

Counts <- read.table(localFilePath, header = TRUE)
Samples <- read.table(sampleListPath, header = TRUE)

## modify colnames(Counts) to match Samples format
# Remove leading X
colnames(Counts) <- sub("X", "", colnames(Counts))
colnames(Counts) <- sub("[[:punct:]]", "-", colnames(Counts))

setdiff(colnames(Counts), Samples[,1]) # none
setdiff(Samples[,1], colnames(Counts)) # NA07-139 is missing from Counts.


# and transcript files

localFilePath  <- '/tmp/AMP-AD_MSBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_transcript_id_counts.txt'
sampleListPath  <- '/tmp/TCXSamples.csv'

Counts <- read.table(localFilePath, header = TRUE)
Samples <- read.table(sampleListPath, header = TRUE)

## modify colnames(Counts) to match Samples format
# Remove leading X
colnames(Counts) <- sub("X", "", colnames(Counts))
colnames(Counts) <- sub("[[:punct:]]", "-", colnames(Counts))

setdiff(colnames(Counts), Samples[,1]) # none
setdiff(Samples[,1], colnames(Counts)) # NA08-162 is missing from Counts.
