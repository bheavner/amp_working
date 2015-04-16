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

#file <- "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates.csv"
#synID <- "syn3483880"

uploadToSynapse <- function(file, synID) {

  # get info about current synapse object to modify
  metaDataOnly <- synGet(synID, downloadFile=F)
  
  # make a new object with the file for editing and upload
  objectAnnotations <- synGetAnnotations(metaDataOnly)
  objectProperties <- synGetProperties(metaDataOnly)
  objectActivity  <- synGetActivity(metaDataOnly)
  
  fileObject <- File(path = file,
                     parentId = objectProperties$parentId,
                     name = "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates")
  
  # define annotations as a list
  metaDataAnnotations <- list(consortium = "AMP-AD",
                              study = "TAUAPPms",
                              center = "UFL-Mayo-ISB",
                              dataType = "Covariates",
                              disease = "Alzheimers Disease",
                              platform = "IlluminaHiSeq2000",
                              mouseModel = "APP",
                              tissueTypeAbrv = "FB",
                              tissueType = "Forebrain",
                              fileType = "csv",
                              organism = "Mus musculus",
                              dataContact = "Ben Heavner <ben.heavner@systemsbiology.org>")
  
  # set the annotations of the fileObject
  synSetAnnotations(fileObject) <- metaDataAnnotations
  
  # define activity
  
  # code used to make this covariates file
  codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.1/UFL/APP/makeAPPCovariates.R")    
  
  act <- Activity(name='Covariate file generation',
                       executed = as.list(codeFile))

  generatedBy(fileObject) = act
    
  # Upload to synapse without incrementing version
  fileObject <- synStore(fileObject, forceVersion = T)
  
  return(fileObject)
}