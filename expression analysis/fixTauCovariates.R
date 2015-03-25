# Goal: fix Tau covariates file so sample IDs there are the same as in the
# RNA-seq data (that is, modify Mouse_IDs in syn2875343)

library(dplyr)
library(synapseClient)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# get the normalized readcount and covariates files from synapse
covariatesFile <- synGet('syn2875343')

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

# First, replace . with _
covariates$Mouse_ID <- sub("[.]", "_", covariates$Mouse_ID)

# add _ after "LP"
covariates$Mouse_ID <- sub("LP", "LP_", covariates$Mouse_ID)

# next, I want to replace - with . in all cases
covariates$Mouse_ID <- sub('-', '.', covariates$Mouse_ID)

# finally, I want to prepend all IDs that start with a number with "X". My regex-fu is weak:
covariates$Mouse_ID <- sub('^1', 'X1', covariates$Mouse_ID)
covariates$Mouse_ID <- sub('^3', 'X3', covariates$Mouse_ID)

# write the updated table
write.table(covariates, file="Mouse_samples.txt", quote=F)

# upload to synapse (provenance undefined)


## later (3/24/15), discover one more:
#in count data, sample “LP62_4”
#in covariates, sample “LP_62_4”
# fix: change covariates sample to match count data ID.
covariatesFile <- synGet('syn2875343')

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

# First, replace . with _
covariates$Mouse_ID <- sub("LP_62_4", "LP62_4", covariates$Mouse_ID)

# write the updated table
write.table(covariates, file="Mouse_samples.txt", quote=F)

# upload to synapse (provenance undefined)