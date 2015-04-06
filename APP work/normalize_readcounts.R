# This script normalizes the merged count_files produced by running merge_file_counts.R
# on SNAPR output. It does so by dowloading the merged read count files, processing them,
# and uploading new synapse objects with appropriate provenance.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/TODO")

# The files for this batch are:

# TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts ('syn3483519')
# TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts ('syn3483625')

countFileSynapseIDs <- c('syn3483519','syn3483625') 

for (mergedCountFile in countFileSynapseIDs) {
    message("Normalizing ", mergedCountFile)

    # Download file from Synapse
    originalCountFile <- synGet(mergedCountFile)

    # unzip file and load for processing
    localFilePath <- getFileLocation(originalCountFile)

    if(!file.exists(substr(localFilePath, 1, nchar(localFilePath) - 3))) {
        gunzip(localFilePath)
    }

    localFilePath <- sub('.gz', '', localFilePath) #trim the .gz suffix

    Counts <- read.table(localFilePath, header = TRUE)

    # make DGEList object
    expr <- DGEList(Counts, group = rep(1, ncol(Counts)))

    # calculate normalization factors
    normFactors <- calcNormFactors(expr, method = ("TMM"))

    # use normaliztion factors to calculate cpm -
    # per https://www.biostars.org/p/84087/, that's calculated as
    # count / (library size * normalization factor))

    normalizedCpm <- cpm(normFactors)

    # write the data to local dir

    newFileName <- sub('_transposed.txt.gz', '', originalCountFile$properties$name) #legacy?
    newFileName <- paste0(newFileName, "_normalized.txt", sep="")

    write.table(normalizedCpm, newFileName, quote = FALSE, sep = "\t", row.names = TRUE)

    # package it up, then create a Synapse object for the output file and upload with provenance

    gzip(newFileName)

    newFileName <- paste0(newFileName, ".gz", sep="")

    parentId <- originalCountFile$properties$parentId

    normalizedCountFile <- File(newFileName, parentId = parentId)

    normalizedCountFile <- synStore(normalizedCountFile,
                                    activityName="CPM (using TMM) from edgeR normalization",
                                    used=list(list(name = "normalize_readcounts.R",
                                                   url = codeFile, wasExecuted = T),
                                              list(entity=originalCountFile,
                                                   wasExecuted=F)))
}
