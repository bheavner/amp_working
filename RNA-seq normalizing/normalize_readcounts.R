# This script normalizes the merged count_files produced by running merge_file_counts.R on SNAPR output. It does so by dowloading the merged read count files, processing them, and uploading new synapse objects with appropriate provenance.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/bheavner/ampSynapseProjects/blob/master/rnaseqAnalysis/normalize_readcounts.R")

# The files for the first batch are:

# ad_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160436')
# ad_pilot_rnaseq_gene_name_counts.txt.gz ('syn3160433')
# ad_pilot_rnaseq_transcript_id_counts.txt.gz ('syn3160437')

# psp_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160443')
# psp_pilot_rnaseq_gene_name_counts.txt.gz ('syn3160442')
# psp_pilot_rnaseq_transcript_id_counts.txt.gz ('syn3160444')

# mouse_tau_rnaseq_gene_id_counts.txt.gz ('syn3160706')
# mouse_tau_rnaseq_gene_name_counts.txt.gz ('syn3160705')
# mouse_tau_rnaseq_transcript_id_counts.txt.gz ('syn3160709')

countFileSynapseIDs <- c('syn3160436', 'syn3160433', 'syn3160437', 'syn3160443', 'syn3160442', 'syn3160444', 'syn3160706', 'syn3160705','syn3160709' )

for (mergedCountFile in countFileSynapseIDs) {

  message("Normalizing ", mergedCountFile)

  # Download file from Synapse
  originalCountFile <- synGet(mergedCountFile)

  # unzip file and load for processing
  localFilePath <- getFileLocation(originalCountFile)

  if(!file.exists(substr(localFilePath, 1, nchar(localFilePath)-3))) {
    gunzip(localFilePath)
  }

  localFilePath <- substr(localFilePath, 1, nchar(localFilePath)-3) #trim the .gz suffix

  rawCounts <- read.table(localFilePath, header = TRUE)

  # begin processing - first transpose b/c James did it differently than DGEList expects
  transposedCounts <- t(rawCounts)

  # make DGEList object
  expr <- DGEList(transposedCounts, group = rep(1,ncol(transposedCounts)))

  # calculate normalization factors
  normFactors <- calcNormFactors(expr, method=("TMM"))

  # use normaliztion factors to calculate cpm -
  # per https://www.biostars.org/p/84087/, that's calculated as
  # count / (library size * normalization factor))

  normalizedCpm <- cpm(normFactors)

  # make it a data frame and write the table to local dir
  normalizedCounts <- as.data.frame(t(normalizedCpm))

  newFileName <- substr(originalCountFile$properties$name, 1, nchar(originalCountFile$properties$name)-7) # strip .txt.gz
  newFileName <- paste(newFileName, "_normalized.txt", sep="")

  write.table(normalizedCounts, newFileName, quote = FALSE, row.names = FALSE)

  # package it up, then create a Synapse object for the output file and upload with provenance

  gzip(newFileName)

  newFileName <- paste(newFileName, ".gz", sep="")

  parentId <- originalCountFile$properties$parentId

  normalizedCountFile <- File(newFileName, parentId = parentId)

  normalizedCountFile <- synStore(normalizedCountFile, activityName="CPM (using TMM) from edgeR normalization", used=list(list(name = "normalize_readcounts.R", url = codeFile, wasExecuted = T), list(entity=originalCountFile, wasExecuted=F)))
}
