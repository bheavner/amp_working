#goal: build list of public and working synapse IDs to help track progress on project deliverables

require(synapseClient)
synapseLogin()

#query the knowledge portal for all the mayo/ufl/isb data
df <- synQuery('select name,id from file where projectId==\'syn2580853\' and center==\'UFL-Mayo-ISB\'')

#grab all the synapse objects for all mayo/ufl/isb data, but don't download the data
synapseEntity<-lapply(df$file.id,function(x){return(synGet(x,downloadFile = F))})

#extract internal Synapse id from provenance of all files, and add it to the df table
mayoProvenance <- lapply(synapseEntity,function(x){return(synGetActivity(x))})
file.oldId <- sapply(mayoProvenance,function(x){return(sapply(x$used,function(x){return(x$reference$targetId)}))})

names(file.oldId) <- df$file.id
