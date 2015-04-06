# This is a reusable function to upload covariates file to Synapse with annotation and provenance as described at syn3382222.
# it assumes that there's already a synapse ID created by a previous upload, but none of the existing data is necessarily correct.

# Example usage:
# 
# source("uploadToSynapse.R")
#
# file <- "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates.csv"
# synID <- syn3483880"
# 
# fileObject <- uploadToSynapse(file, synID)

file <- "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates.csv"
synID <- "syn3483880"

uploadToSynapse <- function(file, synID) {

  # get info about current synapse object to modify
  metaDataOnly <- synGet(synID, downloadFile=F)
  
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
  
  # Add all required annotation
  objectAnnotations$consortium <- "AMP-AD"
  objectAnnotations$study <- "TAUAPPms"
  objectAnnotations$center <- "UFL-Mayo-ISB"
  objectAnnotations$dataType <- "Covariates"
  objectAnnotations$disease <- "Alzheimers Disease"
  objectAnnotations$platform <- "IlluminaHiSeq2000"
  # objectAnnotations$other # NA
  objectAnnotations$mouseModel <- "APP"
  # objectAnnotations$imputationReference # NA
  objectAnnotations$tissueTypeAbrv <- "FB" # IS THIS OKAY?
  objectAnnotations$tissueType <- "Forebrain" # IS THIS OKAY? (from syn3157182)
  # objectAnnotations$dataSubType # I THINK NA
  objectAnnotations$fileType <- "csv"
  objectAnnotations$organism <- "Mus musculus"
  objectAnnotations$dataContact <- "Ben Heavner <ben.heavner@systemsbiology.org>"
  
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