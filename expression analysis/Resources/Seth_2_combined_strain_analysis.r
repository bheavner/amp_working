#preliminaries
#set working directory
setwd("/proj/price1/sament/hdlux")
options(stringsAsFactors=F)
#load libraries needed in this analysis
library(limma)
library(WGCNA)
library(sva)
library(splines)
allowWGCNAThreads()

#expression data
load("hdlux_ComBat_corrected_expr.RData")

#column metadata
load("metadata_after_outlier_removal.RData")

#probe annotation
anno = read.delim("/proj/price1/sament/hdlux/HDLux_agilent_gene_list-2.txt")
matchProbes = match( rownames(ComBatExpr) , anno$ProbeID )
anno = anno[ matchProbes , ]

#which samples go into the current analysis?
useSamples = which(traitData2$allele_nominal %in% c("Q111" , "WT") )

expr.strain = ComBatExpr[,useSamples]
meta.strain = traitData2[useSamples,]
attach(meta.strain)
table( allele_nominal , week , sex )

strain = factor( meta.strain$strain )

#recode variables from the metadata file to 'factors', needed by limma for the linear model
#subset the samples by age, according to bins
age.quartiles = quantile( meta.strain$age )
age.bins = vector(length = nrow(meta.strain) )
age.bins[ meta.strain$age <= age.quartiles[3] ] = "young"
age.bins[ meta.strain$age > age.quartiles[3] ] = "old"

allele = as.factor(meta.strain$allele_nominal)
allele = relevel(allele , ref = "WT")

sex = factor(meta.strain$sex)

group = paste( strain , allele , sex , age.bins , sep="." )
group = factor( group , levels = c( 
	"B6.WT.F.young" , "B6.Q111.F.young" ,
	"B6.WT.M.young" , "B6.Q111.M.young" ,
	"B6.WT.F.old"  , "B6.Q111.F.old" ,
	"B6.WT.M.old"  , "B6.Q111.M.old" ,
	"CD1.WT.F.young" , "CD1.Q111.F.young" ,
	"CD1.WT.M.young" , "CD1.Q111.M.young" ,
	"CD1.WT.F.old"  , "CD1.Q111.F.old" ,
	"CD1.WT.M.old"  , "CD1.Q111.M.old" ,
	"FVB.WT.F.young" , "FVB.Q111.F.young" ,
	"FVB.WT.M.young" , "FVB.Q111.M.young" ,
	"FVB.WT.F.old"  , "FVB.Q111.F.old" ,
	"FVB.WT.M.old"  , "FVB.Q111.M.old" ,
	"129.WT.F.young" , "129.Q111.F.young" ,
	"129.WT.M.young" , "129.Q111.M.young" ,
	"129.WT.F.old"  , "129.Q111.F.old" ,
	"129.WT.M.old"  , "129.Q111.M.old" )
	 )

#now run a linear model using the variables that we have created
#the design matrix contains the levels of each sample for each variable
design = model.matrix( ~ 0 + group )
colnames(design) = c( 
	"B6.WT.F.young" , "B6.Q111.F.young" ,
	"B6.WT.M.young" , "B6.Q111.M.young" ,
	"B6.WT.F.old"  , "B6.Q111.F.old" ,
	"B6.WT.M.old"  , "B6.Q111.M.old" ,
	"CD1.WT.F.young" , "CD1.Q111.F.young" ,
	"CD1.WT.M.young" , "CD1.Q111.M.young" ,
	"CD1.WT.F.old"  , "CD1.Q111.F.old" ,
	"CD1.WT.M.old"  , "CD1.Q111.M.old" ,
	"FVB.WT.F.young" , "FVB.Q111.F.young" ,
	"FVB.WT.M.young" , "FVB.Q111.M.young" ,
	"FVB.WT.F.old"  , "FVB.Q111.F.old" ,
	"FVB.WT.M.old"  , "FVB.Q111.M.old" ,
	"s129.WT.F.young" , "s129.Q111.F.young" ,
	"s129.WT.M.young" , "s129.Q111.M.young" ,
	"s129.WT.F.old"  , "s129.Q111.F.old" ,
	"s129.WT.M.old"  , "s129.Q111.M.old" )

#fit an initial linear model with the variables in design
fit5 = lmFit( expr.strain , design )

#to test the significance of particular contrasts we care about, we create a contrast matrix, designating comparisons between the levels of the variables
contr.mat = makeContrasts( 
	Q111vsWT = (B6.Q111.F.young+CD1.Q111.F.young+s129.Q111.F.young+FVB.Q111.F.young+B6.Q111.M.young+CD1.Q111.M.young+s129.Q111.M.young+FVB.Q111.M.young+B6.Q111.F.old+CD1.Q111.F.old+s129.Q111.F.old+FVB.Q111.F.old+B6.Q111.M.old+CD1.Q111.M.old+s129.Q111.M.old+FVB.Q111.M.old)/16 - (B6.WT.F.young + CD1.WT.F.young + s129.WT.F.young + FVB.WT.F.young + B6.WT.M.young+CD1.WT.M.young+s129.WT.M.young+FVB.WT.M.young+B6.WT.F.old+CD1.WT.F.old+s129.WT.F.old+FVB.WT.F.old+B6.WT.M.old+CD1.WT.M.old+s129.WT.M.old+FVB.WT.M.old)/16 ,
	Q111vsWT.young = (B6.Q111.F.young+CD1.Q111.F.young+s129.Q111.F.young+FVB.Q111.F.young+B6.Q111.M.young+ CD1.Q111.M.young+s129.Q111.M.young+FVB.Q111.M.young)/8 - (B6.WT.F.young + CD1.WT.F.young + s129.WT.F.young+FVB.WT.F.young+B6.WT.M.young+CD1.WT.M.young+s129.WT.M.young+FVB.WT.M.young)/8 ,
	Q111vsWT.old = (B6.Q111.F.old+ CD1.Q111.F.old+ s129.Q111.F.old+ FVB.Q111.F.old+ B6.Q111.M.old+ CD1.Q111.M.old+ s129.Q111.M.old+ FVB.Q111.M.old)/8 - (B6.WT.F.old+CD1.WT.F.old+s129.WT.F.old+FVB.WT.F.old+B6.WT.M.old+CD1.WT.M.old+s129.WT.M.old+FVB.WT.M.old)/8 , 
	Q111vsWT.F.young = (B6.Q111.F.young+CD1.Q111.F.young+s129.Q111.F.young+FVB.Q111.F.young)/4 - (B6.WT.F.young + CD1.WT.F.young + s129.WT.F.young + FVB.WT.F.young)/4 ,
	Q111vsWT.M.young = (B6.Q111.M.young+CD1.Q111.M.young+s129.Q111.M.young+FVB.Q111.M.young)/4 - (B6.WT.M.young+CD1.WT.M.young+s129.WT.M.young+FVB.WT.M.young)/4 , 
	Q111vsWT.F.old = (B6.Q111.F.old+CD1.Q111.F.old+s129.Q111.F.old+FVB.Q111.F.old)/4 - (B6.WT.F.old+CD1.WT.F.old+s129.WT.F.old+FVB.WT.F.old)/4 ,
	Q111vsWT.M.old = (B6.Q111.M.old+CD1.Q111.M.old+s129.Q111.M.old+FVB.Q111.M.old)/4 - (B6.WT.M.old+CD1.WT.M.old+s129.WT.M.old+FVB.WT.M.old)/4 ,
	Q111vsWT.B6.old = (B6.Q111.M.old+B6.Q111.F.old)/2 - (B6.WT.M.old+B6.WT.F.old)/2 ,
	Q111vsWT.CD1.old = (CD1.Q111.M.old+CD1.Q111.F.old)/2 - (CD1.WT.M.old+CD1.WT.F.old)/2 ,
	Q111vsWT.FVB.old = (FVB.Q111.M.old+FVB.Q111.F.old)/2 - (FVB.WT.M.old+FVB.WT.F.old)/2 ,
	Q111vsWT.s129.old = (s129.Q111.M.old+s129.Q111.F.old)/2 - (s129.WT.M.old+s129.WT.F.old)/2,
	Q111vsWT.B6 = (B6.Q111.M.old+B6.Q111.F.old+B6.Q111.M.young+B6.Q111.F.young)/4 - (B6.WT.M.old+B6.WT.F.old+B6.WT.M.young+B6.WT.F.young)/4 ,
	Q111vsWT.CD1 = (CD1.Q111.M.old+CD1.Q111.F.old+CD1.Q111.M.young+CD1.Q111.F.young)/4 - (CD1.WT.M.old+CD1.WT.F.old+CD1.WT.M.young+CD1.WT.F.young)/4 ,
	Q111vsWT.FVB = (FVB.Q111.M.old+FVB.Q111.F.old+FVB.Q111.M.young+FVB.Q111.F.young)/4 - (FVB.WT.M.old+FVB.WT.F.old+FVB.WT.M.young+FVB.WT.F.young)/4 ,
	Q111vsWT.s129 = (s129.Q111.M.old+s129.Q111.F.old+s129.Q111.M.young+s129.Q111.F.young)/4 - (s129.WT.M.old+s129.WT.F.old+s129.WT.M.young+s129.WT.F.young)/4 ,
	sex = (B6.WT.F.young + B6.Q111.F.young + B6.WT.F.old + B6.Q111.F.old +
		CD1.WT.F.young + CD1.Q111.F.young + CD1.WT.F.old  + CD1.Q111.F.old +
		FVB.WT.F.young + FVB.Q111.F.young + FVB.WT.F.old  + FVB.Q111.F.old +
		s129.WT.F.young + s129.Q111.F.young + s129.WT.F.old  + s129.Q111.F.old )/16 -
		(B6.WT.M.young + B6.Q111.M.young + B6.WT.M.old + B6.Q111.M.old +
		CD1.WT.M.young + CD1.Q111.M.young + CD1.WT.M.old  + CD1.Q111.M.old +
		FVB.WT.M.young + FVB.Q111.M.young + FVB.WT.M.old  + FVB.Q111.M.old +
		s129.WT.M.young + s129.Q111.M.young + s129.WT.M.old  + s129.Q111.M.old )/16 ,
	age = (B6.WT.F.young + B6.Q111.F.young + CD1.WT.F.young + CD1.Q111.F.young +
		FVB.WT.F.young + FVB.Q111.F.young +s129.WT.F.young + s129.Q111.F.young +
		B6.WT.M.young + B6.Q111.M.young + CD1.WT.M.young + CD1.Q111.M.young + 
		FVB.WT.M.young + FVB.Q111.M.young + s129.WT.M.young + s129.Q111.M.young)/16 -
		(B6.WT.F.old + B6.Q111.F.old + CD1.WT.F.old + CD1.Q111.F.old +
		FVB.WT.F.old + FVB.Q111.F.old +s129.WT.F.old + s129.Q111.F.old +
		B6.WT.M.old + B6.Q111.M.old + CD1.WT.M.old + CD1.Q111.M.old + 
		FVB.WT.M.old + FVB.Q111.M.old + s129.WT.M.old + s129.Q111.M.old)/16 ,
	levels = design
)

#now we update the linear model with these contrasts.
fit6 = contrasts.fit( fit5 , contr.mat )
#and we use an empirical Bayesian modeler to test for the significance of each probe.
fit6 = eBayes(fit6)


nProbes = nrow(results)
combined.results = data.frame( probeID = rownames(results) )
for( i in 1:ncol(results) ) {
	top = topTable( fit6 , coef=i , number = nProbes , sort.by = "none" , confint=T)
	colnames(top) = paste( colnames(results)[i] , colnames(top) , sep=".")
	combined.results = cbind(combined.results , top[ , c(2,3,4,7,8) ] )
}
combined.results = cbind( anno[ , c(1,3,6) ] , combined.results[,-1] )

write.table( combined.results , file="hdlux_all_strains_limma_fit_pvals.txt" , row.names=F  , quote=F , sep="\t")
save(combined.results , file="hdlux_all_strains_limma_fit_pvals.RData")

sig.b6 = topTable( fit6 , coef=8 , number=nrow(results) )
sig.cd1 = topTable( fit6 , coef=9 , number=nrow(results) )
sig.combined = topTable( fit6 , coef = 8:11 , number = nrow(results ) )

a = which( sig.combined$adj.P < 0.05 )

pdf("Comparison of Fold Changes Between Strains")
par(mar=c(2,2,1,1) , bty="l", oma=c(5,5,5,5))
layout( matrix( c( 1 , 2 , 3 , 0 , 4 , 5 , 0 , 0 , 6 ) , 3 , 3 ) , ) 
plot( sig.combined[a,2] , sig.combined[a,3] , xlim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ),  ylim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ) , ylab = "CD1")
plot( sig.combined[a,2] , sig.combined[a,4] , xlim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ),  ylim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ) , ylab = "FVB")
plot( sig.combined[a,2] , sig.combined[a,5] , xlim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ),  ylim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ) )
plot( sig.combined[a,3] , sig.combined[a,4] , xlim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ),  ylim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ))
plot( sig.combined[a,3] , sig.combined[a,5] , xlim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ),  ylim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ))
plot( sig.combined[a,4] , sig.combined[a,5] , xlim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ),  ylim=c( min( sig.combined[,2:5] ) , max( sig.combined[,2:5] ) ) )
mtext(outer=T , adj = 0.13 , side=1 , "C57BL/6" , line=1.5)
mtext(outer=T , adj = 0.52 , side=1 , "CD1" , line=1.5)
mtext(outer=T , adj = 0.87 , side=1 , "FVB" , line=1.5)
mtext(outer=T , adj = 0.14 , side=2 , "129" , line=1.5)
mtext(outer=T , adj = 0.52 , side=2 , "FVB" , line=1.5)
mtext(outer=T , adj = 0.87 , side=2 , "CD1" , line=1.5)
mtext(outer=T , side=3 , cex=.8 , line=0.5 , "1551 Probes , FDR < 0.05 in Old Mice")
mtext(outer=T , side=3 , cex=1.5 , line=2 , "Comparison of Fold Changes Between Strains")
dev.off()

probes = which( results[,4] == 1 )

sort(unique(anno[ probe , 3 ]))

#probes = sig[1:100,1]
sig = topTable( fit6 , number=nrow(results) , coef=1)
x = merge( sig[1:100 , ] , anno , by.x = 1 , by.y = 1 , all.x = T)

x = x[ order(x[,5]) ,  ]

#probe = x[ x$GeneSymbol == "Htt" , 1 ]
#probe = x[ 25 , 1 ]

anno[which(anno[,3] == "Klf16") , ]
probe = "A_52_P326354"
probe = "A_55_P1992676"
probe = "A_51_P215097"
pdf("Expression Patterns of Probes with 100 Lowest Q111vsWT p-values.pdf" , width=11 , height = 7)
for( probe in x[,1] ) {
	par( oma=c(2,2,4,0) , mar=c(10,2,1,0) )
	layout( matrix(c(1,1,1,1,2,3,4,5),2,4, byrow=T) )
	boxplot( expr.strain[ probe , ] ~ group , las=3 , cex=0.8 , col=c("blue" , "orange"))
	mtext( anno[ which(anno[,1] == probe)  , 3 ] , cex=1.5 , font=4 , line=1 , outer=T , las=1)
	mtext( paste( x[ x[,1]==probe, 12 ] , "; P =" , signif( x[ x[,1] == probe ,5] ,  digits=2) ) , cex=1 , line = -0.5 , outer=T , font=3 , las=1)
	for( i in c("B6" , "CD1" , "FVB" , "129" ) ) {
		par(mar=c(3,3,1,0))
		plot( meta.strain$age , expr.strain[probe,] , type="n" , bty="l" , xlab="Age (weeks)" ,ylab= paste( "Expression in" , i ) ) 
		wt = which( meta.strain$allele_nominal == "WT" & meta.strain$strain == i)
		points( meta.strain$age[wt] , expr.strain[ probe , wt ] , col="blue" , pch=24)
		lines( lowess( x =meta.strain$age[wt] , y = expr.strain[ probe , wt ] ), col="blue" , lty=1 , lwd=2)
		q111 = which(meta.strain$allele_nominal == "Q111" & meta.strain$strain == i)
		points( meta.strain$age[q111] , expr.strain[probe,q111]  , pch=19 , col="orange")
		lines( lowess( x =meta.strain$age[q111] , y = expr.strain[probe,q111] ), col="orange" , lty=1, lwd=2)
		legend( x=min( meta.strain$age) , y=max( expr.strain[probe,]), legend = c("WT" , "Q111") , col=c("blue" , "orange") , pch=c(24,19) , lty = c(1,1))
		mtext( side=3 , line=0 , paste( "Strain" , i ) , las=1)
	}
	mtext( side=1 , "Age (weeks)" , outer=T , las=1)
	mtext( side=2 , "log2(Expression)" , outer=T , las=3)
}
dev.off()


contr.mat2 = makeContrasts( 
	Q111vsWT.B6 = (B6.Q111.M.old+B6.Q111.F.old+B6.Q111.M.young+B6.Q111.F.young)/4 - (B6.WT.M.old+B6.WT.F.old+B6.WT.M.young+B6.WT.F.young)/4 ,
	Q111vsWT.CD1 = (CD1.Q111.M.old+CD1.Q111.F.old+CD1.Q111.M.young+CD1.Q111.F.young)/4 - (CD1.WT.M.old+CD1.WT.F.old+CD1.WT.M.young+CD1.WT.F.young)/4 ,
	Q111vsWT.FVB = (FVB.Q111.M.old+FVB.Q111.F.old+FVB.Q111.M.young+FVB.Q111.F.young)/4 - (FVB.WT.M.old+FVB.WT.F.old+FVB.WT.M.young+FVB.WT.F.young)/4 ,
	Q111vsWT.s129 = (s129.Q111.M.old+s129.Q111.F.old+s129.Q111.M.young+s129.Q111.F.young)/4 - (s129.WT.M.old+s129.WT.F.old+s129.WT.M.young+s129.WT.F.young)/4 ,
	levels = design )
fit2 = contrasts.fit( fit5 , contr.mat2 )
fit2 = eBayes(fit2)


#these functions summarize the results of the linear model.
results = decideTests( fit2 , method = "global" , adjust.method = "BH" , p.value = 0.1 )
nSig = colSums(abs(results))
pdf("Venn_Diagram_FDR05_by_strain.pdf")
par(mfrow=c(1,1) , oma=c(0,0,0,0))
vennDiagram( results , names=c("C57BL/6" , "CD1" , "FVB" , "129") )
dev.off()





