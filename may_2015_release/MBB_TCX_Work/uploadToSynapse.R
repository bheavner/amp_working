# This is a reusable function to upload covariates file to Synapse with annotation and provenance as described at syn3382222.
# it assumes that there's already a synapse ID created by a previous upload, but none of the existing data is necessarily correct.

# Example usage:
# 
# source("uploadToSynapse.R")
#
# file <- "AMP-AD_MSBB_UFL-Mayo-ISB_IlluminaHiSeq2000_MSBBTCX-Covariates.csv"
# synID <- syn3509411"
# 
# fileObject <- uploadToSynapse(file, synID)

file <- "AMP-AD_MSBB_UFL-Mayo-ISB_IlluminaHiSeq2000_MSBBTCX-Covariates.csv"
synID <- "syn3509411"

uploadToSynapse <- function(file, synID) {
  
  # get info about current synapse object to modify
  metaDataOnly <- synGet(synID, downloadFile=F)
  
  
  # to clean up
  metaDataOnly@filePath = '/Users/bheavner/foo/some_file_path.foo'
# also change name in properties
  synStore(metaDataOnly) # will version it this way, which you must to update file
  
  # make a new object with the file for editing and upload
  objectAnnotations <- synGetAnnotations(metaDataOnly)
  objectProperties <- synGetProperties(metaDataOnly)
  objectActivity  <- synGetActivity(metaDataOnly)
  
  fileObject <- File(path = file,
                     parentId = objectProperties$parentId,
                     annotations = objectAnnotations,
                     generatedBy = objectActivity,
                     name = "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates"
  )
  
  # Add all required annotation (described at syn3382222)
  objectAnnotations$consortium <- "AMP-AD"
  objectAnnotations$study <- "MSBB"
  objectAnnotations$center <- "UFL-Mayo-ISB"
  objectAnnotations$dataType <- "Covariates"
  objectAnnotations$disease <- c("Progressive Supranuclear Palsy", "Alzheimers Disease")
  objectAnnotations$platform <- "IlluminaHiSeq2000"
  # objectAnnotations$other # NA
  # objectAnnotations$mouseModel #NA
  # objectAnnotations$imputationReference # NA
  objectAnnotations$tissueTypeAbrv <- "TCX"
  objectAnnotations$tissueType <- "Temporal Cortex"
  # objectAnnotations$dataSubType  <- "patient" # NOT SURE IF RIGHT
  objectAnnotations$fileType <- "csv"
  objectAnnotations$organism <- "Homo sapiens"
  objectAnnotations$dataContact <- "Nilufer Ertekin-Taner, M.D., Ph.D. Taner.Nilufer@mayo.edu"
  
  synSetAnnotations(fileObject) <- objectAnnotations
  
  # Add provenance
  
  # code used to make this covariates file
  #codeFile <- ("https://github.com/PriceLab/AMP-ad/APP_work/tree/v0.1-alpha/makeAPPCovariates.R")    
  
  #  act <- Activity(name='HBTRC Reference Data Migration',
  #                  used=list(list(entity=emoryTable@values$originalSynapseId[i],wasExecuted=F)),
  #                  executed=list("https://github.com/Sage-Bionetworks/ampAdScripts/blob/master/Emory/migrateEmoryFeb2015.R"))
  
  #  generatedBy(hbtrcSyn) <- act
  
  # Upload to synapse without incrementing version
  fileObject <- synStore(fileObject, forceVersion = T)
  
  return(fileObject)
}