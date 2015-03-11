# goal: download the merged read count files from Synapse, normalize with edgeR, and upload the results to new Synapse objects (with provenance)

# dowload merged read count files. I _think_ the files are:
# ad_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160436')
# ad_pilot_rnaseq_gene_name_counts.txt.gz ('syn3160433')
# ad_pilot_rnaseq_transcript_id_counts.txt.gz ('syn3160437')

# psp_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160443')
# psp_pilot_rnaseq_gene_name_counts.txt.gz ('syn3160442')
# psp_pilot_rnaseq_transcript_id_counts.txt.gz ('syn3160444')

# mouse_tau_rnaseq_gene_id_counts.txt.gz ('syn3160706')
# mouse_tau_rnaseq_gene_name_counts.txt.gz ('syn3160705')
# mouse_tau_rnaseq_transcript_id_counts.txt.gz ('syn3160709')

# QUESTION FOR JAMES: Should there be a transcript_name_counts file, too? (he says no/optional)

setwd("~/Projects/UO1-AMP/RNA-seq normalizing")
require("synapseClient")
require("R.utils")
require("edgeR")

synapseLogin()

# I'll want to make this a loop, but for now, one at a time..
#input_files <- c("syn2875349", "syn2875350", "syn2875347")
# for (input_file in input_files) { }

ad_pilot_rnaseq_gene_id_counts <- synGet('syn3160436')
localFilePath <- ad_pilot_rnaseq_gene_id_counts@filePath

if(!file.exists(substr(localFilePath, 1, nchar(localFilePath)-3))) {
    gunzip(localFilePath)
}

localFilePath <- substr(localFilePath, 1, nchar(localFilePath)-3)

raw_counts <- read.table(localFilePath, header = TRUE)

## BH playing
str(raw_counts)
# 96 obs of 64253 variables - looks like:
#$ ENSG00000000003: int  530 2395 448 815 310 222 552 1125 1107 928 ...
#$ ENSG00000000005: int  5 4 3 5 4 3 4 5 4 0 ...
#$ ENSG00000000419: int  776 488 713 496 486 526 1081 934 559 620 ...

row.names(raw_counts)
# row.names looks like this:
# "4088"     "4095"     "4808"     "NA00-136" "NA00-175"
# I think that's the sample names

colnames(raw_counts)
# looks like this:
# "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457"
# I think that's the gene IDs

raw_counts[1]
# looks like this:
#ENSG00000000003
#4088                 530
#4095                2395

# I think this means that for geneID ENSG00000000003, sample 4088 has 530 reads, sample 4095 has 2395, sample NA00-136 has 815, etc. (total 96 samples, 64253 geneIDs).

raw_counts[1,1]
# 530
# that's the number of counts for geneID ENSG00000000003 in sample 4088

raw_counts[,1]
# 530 2395  448  815  ...
# that's the counts for geneID ENSG00000000003 in each of the 96 samples


# To use calcNormFactors, I think I want to transform this data frame to a DGEList object, which expects as a minimum "a table of counts (rows=features, columns=samples), group indicator for each column".

# For the table of counts, I want raw_counts(GENES=colnames, SAMPLES=row.names) with values = counts. So there should be 64253 rows, 96 columns, and 6168288 entries. Basically, just a transpose of the data frame.

# to practice, I'll try to make the first two rows of the table:
#         sample 1  sample 2
#gene 1    count    count
#gene 2    count    count

colnames(raw_counts[1:2])
row.names(test)
raw_counts[1:2,1]
raw_counts[1:2,2]

# so it should look like this:
#                  4088  4095
#ENSG00000000003   530   2395
#ENSG00000000005     5      4

# I'll do 3 genes to make it nonsemetric for testing
matrix(data = c(raw_counts[1:2, 1], raw_counts[1:2, 2], raw_counts[1:2, 3]), nrow = 3, ncol = 2, byrow = TRUE, dimnames = list(c(colnames(raw_counts[1:3])), c(row.names(raw_counts)[1:2])))

#                 4088 4095
# ENSG00000000003  530 2395
# ENSG00000000005    5    4
# ENSG00000000419  776  488
#
# woot! now try to make the big one - just transpose the data frame.

transposed_counts <- t(raw_counts)
counts[1:3,1:2]

#                4088 4095
# ENSG00000000003  530 2395
# ENSG00000000005    5    4
# ENSG00000000419  776  488

# Now make a DGEList object (I'm unsure of conditions for this, so assign all to 1 condition)

expr <- DGEList(transposed_counts, group = rep(1,ncol(transposed_counts)))

# normaliztion: CPM (using TMM) from edgeR - If object is a DGEList, then it is returned as output with the relative normalization factors in object$samples$norm.factors.

norm_factors <- calcNormFactors(expr, method=("TMM"))

# normalized counts are made by dividing the counts by the normalization factor. There are 96 normalization factors, so they're by sample, not by gene. So, for example, to normalize the counts for the first sample, I'd do:
counts[,1]/norm_factors$samples$norm.factors[1]

# next: normalize all 96 samples.

normalized_cpm <- cpm(norm_factors)

dim(normalized_cpm) # it's 64253 x 96

# now, repack normalized counts in data frame, export to delimited file, gzip, and upload to synapse with provenance pointing to pretty version of this code on github.

# I want it to look like the orignial file, so when loaded, a data frame with 96 obs of 64253 variables, row names are sample ids, column names are gene ids, and values are normalized counts.

# so normalized_counts[1]
# should look like this:
#ENSG00000000003
#4088               ###
#4095               ###

# I think this means that for geneID ENSG00000000003, sample 4088 has 530 reads, sample 4095 has 2395, sample NA00-136 has 815, etc. (total 96 samples, 64253 geneIDs).

normalized_counts <- as.data.frame(t(normalized_cpm))

normalized_counts[1]
#ENSG00000000003
#4088           11.452189
#4095           63.768737

# so.. that seems to work. Next, write this data frame to a delimited text file, then upload it to synapse.

# Save to file
out_path <- ("ad_pilot_rnaseq_gene_id_counts_normalized.txt")

write.table(normalized_counts, out_path, quote = FALSE, row.names = FALSE)

gzip(file.path)

out_path <- paste(out_path, ".gz", sep="")

# Create a Synapse object for the output file and upload with provenance
parentID = ("syn2875332")

cer_clinical_object <- File(path = out_path,
                            parentId = folder_id)

synSetAnnotations(cer_clinical_object) <- list(model="lm", response="random")

cer_clinical_object <- synStore(cer_clinical_object, activityName="CPM (using TMM) from edgeR normalization", used=list(entity=codeFile, wasExecuted=T), list(entity=exprFile, wasExecuted=F))
