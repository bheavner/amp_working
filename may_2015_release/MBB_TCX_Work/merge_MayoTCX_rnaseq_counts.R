
# This script uses the merging function saved in merge_file_counts.R to combine
# SNAPR count files for reprocessed pilot AD RNAseq data from Mayo

library(synapseClient)
library(tools)
source("merge_count_files.R")

# Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# Define paths for required Synapse objects
mayoTCX_rnaseq_counts <- "syn3632530" # all count files in a zipped directory

# Download files from Synapse
mayoTCX_count_files <- synGet(mayoTCX_rnaseq_counts)
mayoTCX_files_path <- getFileLocation(mayoTCX_count_files)

# Get name of temporary directory to store unzipped files (same as name of
# original compressed directory)
fileDir <- file_path_sans_ext(basename(mayoTCX_files_path))
tmpDir <- tempdir()
unzip(mayoTCX_files_path, exdir = tmpDir)

inputDir <- paste(tmpDir, '/tmp/tcx', sep = "")

#inputDir <- file.path(tmpDir, fileDir)

prefix <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX"

countTypes <- c("gene_id", "transcript_id") #  

for (countType in countTypes) {
  message(paste("Merging", prefix, "files of count type", countType, "..."))
  
  # Create the merged file and store the output file path
  merged_file <- create_merged_file(inputDir, countType, prefix)
  
  # Create a Synapse object for the output file and upload
  merged_file_object <- File(path = merged_file, 
                             parentId = mayoTCX_count_files$properties$parentId)
  merged_file_object <- synStore(merged_file_object)
}
# unlink(inputDir, recursive = TRUE)
