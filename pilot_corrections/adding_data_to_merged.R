# make new readcounts file from snapr output
library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# in shell, use AWSCLI to do this:
# aws s3 cp s3://mayo-prelim-rnaseq/PSP_Samples/ /tmp/psp-rnaseq-counts_Updated/ --recursive --exclude "*.bam*" --exclude "*.fa" --exclude "*.gtf" --include "*.txt" --dryrun
# rezip with:
# zip -r psp-rnaseq-counts_Updated /tmp/psp/psp-rnaseq-counts_Updated

# upload synapse object, get ID, establish provenance (point to https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha)
codeFile <- ("https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha")
synapseId <- ('syn3578144')
metaDataOnly <- synGet('syn2875347', downloadFile=F)

localFilePath <- "/tmp/psp-rnaseq-counts_Updated.zip"

mergedFileObject <- File(path = localFilePath, 
                           parentId = metaDataOnly$properties$parentId,
                           activityName = "SNAPR processing", 
                           used = list(list(name = "SNAPR_tools", 
                                            url = codeFile,
                                            wasExecuted = T)))

mergedFileObject <- synStore(mergedFileObject)

# to fix provenance, redownload metadata, get ID, establish provenance (point to https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha)
synapseId <- ('syn3578144')
metaDataOnly <- synGet(synapseId, downloadFile=F)

attempt <- synStore(metaDataOnly,
                    activityName="SNAPR processing",
                    used=list(list(name = "snapr_tools",
                                   url = codeFile, wasExecuted = T)
                    ))