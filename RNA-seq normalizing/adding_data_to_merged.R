# adding reprocessed snapr output to readcounts
library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

tauZipFile <- synGet('syn2875347')
localFilePath <- getFileLocation(tauZipFile)
localFilePath

# then in shell, go there, unzip, use AWSCLI to do this:
# aws s3 cp s3://ufl-u01-rnaseq/ /tmp/ --recursive --exclude "*" --include "*370763*" --exclude "*.bam*"
# rezip with:
# zip -r *  mouse-tau-rnaseq-counts.zip

# reupload synapse object, get ID, establish provenance (point to https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha)
codeFile <- ("https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha")
synapseId <- ('syn2875347')
metaDataOnly <- synGet('syn2875347', downloadFile=F)

mergedFileObject <- File(path = localFilePath, 
                           parentId = metaDataOnly$properties$parentId,
                           activityName = "SNAPR processing", 
                           used = list(list(name = "SNAPR_tools", 
                                            url = codeFile,
                                            wasExecuted = T)))

mergedFileObject <- synStore(mergedFileObject)

# to fix provenance, redownload metadata, get ID, establish provenance (point to https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha)
synapseId <- ('syn2875347')
metaDataOnly <- synGet(synapseId, downloadFile=F)

attempt <- synStore(metaDataOnly,
                    activityName="SNAPR processing",
                    used=list(list(name = "snapr_tools",
                                   url = codeFile, wasExecuted = T)
                    ))