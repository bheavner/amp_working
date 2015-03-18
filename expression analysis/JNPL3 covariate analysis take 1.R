# Goal: investigate which covariates most influence gene expression levels 
# (see https://www.biostars.org/p/97624/ and http://stats.stackexchange.com/questions/59879/logistic-model-what-is-more-important-anova-chi-sq-test-or-significance-of-coe)

library(dplyr) # for subsetting data
library(R.utils) # for unzipping data
library(synapseClient) # to download data
library(reshape2) # for melt for ggplot

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# get the transposed readcount file and covariates file from synapse
countFile <- synGet('syn3192634') 
covariatesFile <- synGet('syn2875343') # the working dir copy

# unzip count file and load for processing
localCountFilePath <- getFileLocation(countFile)

if(!file.exists(sub('.gz', '', localCountFilePath))) {
  gunzip(localCountFilePath)
}

localCountFilePath <- sub('.gz', '', localCountFilePath) #trim the .gz suffix

counts <- read.table(localCountFilePath, header = TRUE, stringsAsFactors = FALSE)

#PROBLEM: - an rTGMinus sample
#in count data, sample “LP62_4”
#in covariates, sample “LP_62_4”
# so change column name in counts
colnames(counts)[1] <- "LP_62_4"

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), header = TRUE, stringsAsFactors = FALSE)

# subset counts and covariates for JNPL3 samples

JNPL3Samples <- dplyr::select(filter(covariates, Experiment == "MAPT_P301L"), Mouse_ID)

JNPL3Counts <- dplyr::select(counts, one_of(JNPL3Samples$Mouse_ID))

rownames(covariates) <- covariates$Mouse_ID

JNPL3Covariates <- covariates[JNPL3Samples$Mouse_ID,]

# make DGEList object of counts for normalization and filtering

JNPL3 <- DGEList(counts = JNPL3Counts)

# filter data - require minimum of 100 counts per million for at least 2 samples
d.full <- JNPL3 # keep the old one in case we mess up
dim(d.full)
#39179    24

keep <- rowSums(cpm(JNPL3)>100) >= 2
d <- JNPL3[keep,]
dim(d)
#3133   24

# reset library sizes after filtering
d$samples$lib.size <- colSums(d$counts)

# normalize the data using TMM
d <- calcNormFactors(d)

# Now, I want to make a general linear model using covariates, and see which are most significant for predicting expression

# select covariates
colnames(JNPL3Covariates)
#[1] "Mouse_ID"             "Experiment"           "RIN"                  "Genotype"             "Sex"                 
#[6] "Age_months"           "RLIMS.ID"             "Seq.Run.ID"           "Lane.Number"          "Clusters"            
#[11] "Raw_RNAseq_file_name"

# I'd like to test RIN, Genotype, Sex, Age_months, Lane.number, and Clusters
data <- data.frame(t(d$counts))

dim(data)
#[1]   24 3133

dim(JNPL3Covariates)
#[1] 24 11

data$Mouse_ID <- JNPL3Covariates[, "Mouse_ID"]
data$RIN <- JNPL3Covariates[, "RIN"]
data$Genotype <- JNPL3Covariates[, "Genotype"]
data$Sex <- JNPL3Covariates[, "Sex"]
data$Age_months <- JNPL3Covariates[, "Age_months"]
data$Lane.Number <- JNPL3Covariates[, "Lane.Number"]
data$Clusters <- JNPL3Covariates[, "Clusters"]

mdata <- melt(data, id.vars = c("Mouse_ID", "Genotype")

# the dream:
#fit <- glm(data$counts ~ JNPL3Covariates[,"Genotype"], data = mdata)


head(melt(t(d$counts)), id.vars = c(""))

#Var1               Var2 value
#1 X112.10 ENSMUSG00000000088  6036
#2 X112.11 ENSMUSG00000000088  6285
#3 X186848 ENSMUSG00000000088  7493
#4 X307983 ENSMUSG00000000088  5843
#5 X307987 ENSMUSG00000000088  5737
#6 X309581 ENSMUSG00000000088  7515

test <- melt(t(d$counts))
colnames(test) <- c("Mouse_ID", "Gene_Name", "Count")


head(melt(JNPL3Covariates[,c("Mouse_ID", factors)], id.vars = c("Mouse_ID", "Genotype", "Sex", "Age_months", "Lane.Number", "Clusters", "RIN")))

#Mouse_ID Genotype Sex Age_months Lane.Number Clusters RIN
#X112.10  X112.10        -   F         12           6      748 7.8
#X112.11  X112.11        -   F         12           3      761 8.1
#X186848  X186848        -   F         12           6      748 7.5
#X307983  X307983        +   F         12           4      766 7.8
#X307987  X307987        +   F         12           7      701 8.0
#X309581  X309581        +   F         12           7      701 7.1

#... I'm beginning to think I might have to do this for each gene... Try it for just the first one

fit <- glm(ENSMUSG00000000088 ~ RIN + Genotype + Sex + Age_months + Clusters + Lane.Number + Mouse_ID, data = data)
