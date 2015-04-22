# MAY NO LONGER NEED SINCE CHANGES TO MERGE_PSP_PILOT_RNASEQ_COUNTS.R

# This script downloads the merged count_files produced by running merge_file_counts.R on SNAPR output, transposes the
# dataframe so genes are rows and samples are columns, removes samples to be omitted for QC reasons, and reuploads
# them to synapse with appropriate provenance that this was run.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/PriceLab/AMP-ad/pilotCorrections/reformat_merged_readcounts.R")

# The files to process are:

# AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_gene_id_counts_UpdatedID.txt.gz ('syn3578797')
# AMP-AD_MayoPilot_UFL-Mayo-ISB_IlluminaHiSeq2000_TemporalCortex_ProgressiveSupranuclearPalsy_TranscriptCounts_UpdatedID.txt.gz ('syn3582433')

pspCountFileSynapseIDs <- c('syn3578797', 'syn3582433')

# samples to exclude for QC reasons:

pspOmit <- c('1811024561_B', '1811024502_B', '1811024560_B', '1811024331_B')

for (dataset in list(pspCountFileSynapseIDs)) {
    omit <- pspOmit
    
    for (mergedCountFile in dataset) {
        message("Processing SynID: ", mergedCountFile)

        # Download file from Synapse
        originalCountFile <- synGet(mergedCountFile)

        # unzip file and load for processing
        localFilePath <- getFileLocation(originalCountFile)

        if(!file.exists(sub('.gz', '', localFilePath))) {
            gunzip(localFilePath)
        }

        localFilePath <- sub('.gz', '', localFilePath) #trim the .gz suffix

        rawCounts <- read.table(localFilePath, header = TRUE)

        # begin processing - first transpose
        transposedCounts <- t(rawCounts)

        # remove samples that don't pass QC
        keepSamples <- setdiff(colnames(transposedCounts), omit) #probably not the best way..
        transposedCounts <- transposedCounts[, keepSamples]

        # write the data to local dir

        newFileName <- sub('.txt.gz', '', originalCountFile$properties$name)
        newFileName <- paste0(newFileName, "_transposed.txt", sep="")

        write.table(transposedCounts, newFileName, quote = FALSE, sep = "\t", row.names = TRUE)

        # package it up, then create a Synapse object for the output file and upload with provenance

        gzip(newFileName)

        newFileName <- paste0(newFileName, ".gz", sep="")

        parentId <- originalCountFile$properties$parentId

        transposedCountFile <- File(newFileName, parentId = parentId)

        transposedCountFile <- synStore(transposedCountFile,
                                        activityName="Transposed merged readcount file",
                                        used=list(list(name = "reformat_merged_readcounts.R",
                                                       url = codeFile, wasExecuted = T),
                                                  list(entity=originalCountFile,
                                                       wasExecuted=F)))
    }
}
