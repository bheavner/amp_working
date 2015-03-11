# I need to fix provenance as follows:
# syn2875349 - ad-rnaseq-counts.zip #should point to repository
# syn2875350 - psp-rnaseq-counts.zip #should point to repository

# syn3160709 - mouse_tau_rnaseq_transcript_id_counts.txt.gz # should be like syn3160444
# syn3160706 - mouse_tau_rnaseq_gene_id_counts.txt.gz # should be like syn3160444

# syn3192651 - mouse_tau_rnaseq_transcript_id_counts_transposed.txt.gz # should be like syn3192634

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# reupload synapse object, get ID, establish provenance (point to https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha)
# syn2875349 - ad-rnaseq-counts.zip #should point to repository
synapseId <- ('syn2875349')
codeFile <- ("https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha")

metadataOnly <- synGet(synapseId, downloadFile=F)

attempt <- synStore(metadataOnly,
                    activityName="SNAPR processing",
                    used=list(list(name = "snapr_tools",
                                   url = codeFile, wasExecuted = T)
                    ))

# syn2875350 - psp-rnaseq-counts.zip #should point to repository
synapseId <- ('syn2875350')
codeFile <- ("https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha")

metadataOnly <- synGet(synapseId, downloadFile=F)

attempt <- synStore(metadataOnly,
                    activityName="SNAPR processing",
                    used=list(list(name = "snapr_tools",
                                   url = codeFile, wasExecuted = T)
                    ))

# syn3160709 - mouse_tau_rnaseq_transcript_id_counts.txt.gz # should be like syn3160444
synapseId <- ('syn3160709')
codeFile <- ("https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha")

metaDataOnly <- synGet(synapseId, downloadFile=F)

activity_name <- "RNA-seq Count File Merging"
input_files <- c("syn2875347")
code_files <- list(list(name = "merge_mouse_tau_rnaseq_counts.R",
                        url = "https://github.com/jaeddy/ampSynapseProjects/blob/v0.1-alpha/dataEnablement/R/merge_mouse_tau_rnaseq_counts.R", 
                        wasExecuted = T),
                   list(name = "merge_count_files.R",
                        url = "https://github.com/jaeddy/ampSynapseProjects/blob/v0.1-alpha/dataEnablement/R/merge_count_files.R", 
                        wasExecuted = T))
output_files <- c("syn3160709")
description <- paste("To execute run:", 
                     "Rscript merge_mouse_tau_rnaseq_counts.R")

attempt <- synStore(metaDataOnly,
                    used = input_files,
                    executed = code_files,
                    activityName = activity_name,
                    activityDescription = description)


# syn3160706 - mouse_tau_rnaseq_gene_id_counts.txt.gz # should be like syn3160444
synapseId <- ('syn3160706')
codeFile <- ("https://github.com/jaeddy/snapr_tools/tree/v0.1-alpha")

metaDataOnly <- synGet(synapseId, downloadFile=F)

activity_name <- "RNA-seq Count File Merging"
input_files <- c("syn2875347")
code_files <- list(list(name = "merge_mouse_tau_rnaseq_counts.R",
                        url = "https://github.com/jaeddy/ampSynapseProjects/blob/v0.1-alpha/dataEnablement/R/merge_mouse_tau_rnaseq_counts.R", 
                        wasExecuted = T),
                    list(name = "merge_count_files.R",
                         url = "https://github.com/jaeddy/ampSynapseProjects/blob/v0.1-alpha/dataEnablement/R/merge_count_files.R", 
                         wasExecuted = T))
output_files <- c("syn3160706")
description <- paste("To execute run:", 
                      "Rscript merge_mouse_tau_rnaseq_counts.R")

attempt <- synStore(metaDataOnly,
                    used = input_files,
                    executed = code_files,
                    activityName = activity_name,
                    activityDescription = description)


# syn3192651 - mouse_tau_rnaseq_transcript_id_counts_transposed.txt.gz # should be like syn3192634