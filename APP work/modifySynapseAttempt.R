file <- "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates.csv"
synID <- "syn3483880"

# get info about current synapse object to modify
fileObject <- synGet(synID, downloadFile=F)

# define annotations (an example)
fileObjectAnnotations <- list(consortium = "AMP-AD",
                              study = "TAUAPPms")
synSetAnnotations(fileObject) <- fileObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.1/UFL/APP/makeAPPCovariates.R")   

act <- Activity(name='Covariate file generation',
                executed = as.list(codeFile))
generatedBy(fileObject) = act

# define properties
fileProperties <- list(path = file,
                       name = "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates")

synSetProperties(fileObject) <- fileProperties

fileObject <- synStore(fileObject, forceVersion = T)