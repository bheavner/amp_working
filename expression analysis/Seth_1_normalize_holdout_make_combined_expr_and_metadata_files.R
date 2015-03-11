setwd("/proj/price1/sament/hdlux")
options(stringsAsFactors=F)
library(sva)
library(WGCNA)
allowWGCNAThreads()
library(flashClust)

load("/proj/hdlux/data/DataCut.13-11-12/data_setaside.expr_proc_32.RData")
expr = x

load("/proj/hdlux/data/DataCut.13-11-12/metadata_setaside.expr_proc_32.RData")
traitData = x

B6 = which( traitData$strain == "B6" )

table( traitData[B6,"date_hyb"] , traitData[ B6,"allele_nominal"] )

table( traitData[B6,"date_hyb"] , traitData[ B6,"week"] )

table( traitData[,"date_hyb"] , traitData[ ,"allele_nominal"]

#all samples
datExpr0 = t(expr)
sampleTree = flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(10,7)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#cex.axis = 1.5, cex.main = 2)

strain = as.numeric( as.factor(traitData$strain ) )
allele_nominal = as.numeric( as.factor(traitData$allele_nominal))
sex = as.numeric(as.factor(traitData$sex ))
week = traitData$week
date_hyb = as.numeric(as.factor(traitData$date_hyb))
date_harvest = as.numeric(as.factor(traitData$date_harvest))
slide_id = as.numeric(as.factor(traitData$slide_id))
datTraits = data.frame( strain , allele_nominal , sex , week , date_hyb , slide_id , date_harvest)
traitColors = numbers2colors(datTraits, signed = FALSE);
pdf("hdlux_holdout_sampletree_with_trait_heatmap_2015-2-10.pdf")
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
dev.off()

#several samples appear to be outliers.
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0 # [keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#traitData2 = traitData[ keepSamples , ]
traitData2 = traitData
batch = traitData2$date_hyb
n1 = names( which( table(batch) == 1 ) )
batch[ batch %in% n1 ] = "other"
mod = model.matrix(~ as.factor( traitData2$strain ) + as.factor( traitData2$allele_nominal ) + as.factor( traitData2$week ) )
ComBatExpr = ComBat(dat = t(datExpr) , batch = batch , mod = mod , par.prior = TRUE , prior.plots = TRUE )

sampleTree2 = flashClust(dist(t(ComBatExpr)), method = "average");
strain = as.numeric( as.factor(traitData2$strain ) )
allele_nominal = as.numeric( as.factor(traitData2$allele_nominal))
sex = as.numeric(as.factor(traitData2$sex ))
week = traitData2$week
date_hyb = as.numeric(as.factor(traitData2$date_hyb))
date_harvest = as.numeric(as.factor(traitData2$date_harvest))
slide_id = as.numeric(as.factor(traitData2$slide_id))
datTraits2 = data.frame( allele_nominal , sex , week , date_hyb , slide_id , date_harvest)
traitColors2 = numbers2colors(datTraits2, signed = FALSE);
pdf("hdlux_combat_holdout_sampletree_with_trait_heatmap_2015-2-15.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors2, groupLabels = names(datTraits2), main = " ComBat-corrected, no-outlier-removed samples and traits")
dev.off()

save(sampleTree, sampleTree2 , file="sample_clustering_before_and_after_ComBat.RData")
save(ComBatExpr , file="hdlux.holdout.ComBat_corrected_expr.RData")
save(traitData2 , file="metadata.holdout.after_outlier_removal.RData")

expr.holdout = ComBatExpr
meta.holdout = traitData2

#expression data
load("hdlux_ComBat_corrected_expr.RData")
expr.train = ComBatExpr

#column metadata
load("metadata_after_outlier_removal.RData")
meta.train = traitData2

expr.full = cbind( expr.train , expr.holdout )
meta.full = rbind( meta.train , meta.holdout )

system( "mkdir holdout_emats_2015-2-16" )
setwd( "holdout_emats_2015-2-16" )

holdout_file = meta.full[ which( meta.full$setaside == 1 ) , 2 ]
write.table( holdout_file , file = "holdout_file" , col.names=F , row.names=F , quote=F )

strains = unique( meta.full$strain )
for( i in strains ) {
	cat( i )
	cat(": ")
	e = expr.full[ , which( meta.full$strain == i ) ]
	m = meta.full[ which( meta.full$strain == i ) , ]
	colnames(e) = m[,2]
	alleles = setdiff( unique( m$allele_nominal ) , "WT" )
	for( j in alleles ) {
		cat( j )
		cat(", ")
		x = which( m$allele_nominal %in% c("WT",j) )
		e2 = e[,x]
		m2 = m[x,]
		efile = paste( "expr_full." , i , "." , j , "vWT.txt" , sep = "" )
		mfile = paste( "meta_full." , i , "." , j , "vWT.txt" , sep = "" )
		write.table( e2 , file = efile , sep = "\t" , quote=F , col.names=T , row.names=F )
		write.table( m2 , file = mfile , sep = "\t" , quote=F , col.names=T , row.names=F )
	}
	cat("\n")
}

