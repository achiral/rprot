#WGCNA
#install.packages("BiocManager") 
#BiocManager::install("WGCNA") 
#install.packages("flashClust")
#install.packages("lavaan")
#install.packages("sem")
setwd("/Users/user/Dropbox/0_Work/R/WGCNA")

library(WGCNA)
library(cluster)
library(flashClust)
library(lavaan) #sem function
options(stringsAsFactors = FALSE)

#############################################################
#Section 12.1 Constructing a sample network for outlier detection
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv")
dim(femData)
names(femData)
#Remove gene information and transpose the expression data
datExprFemale=as.data.frame(t(femData[, -c(1:8)]))
names(datExprFemale)=femData$substanceBXH
rownames(datExprFemale)=names(femData)[-c(1:8)]
# Now we read in the physiological trait data
traitData = read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)
# use only a subset of the columns
allTraits=traitData[,c(2, 11:15, 17:30, 32:38)]
names(allTraits)
# Order the rows of allTraits so that
# they match those of datExprFemale:
Mice=rownames(datExprFemale)
traitRows = match(Mice, allTraits$Mice)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
# show that row names agree
table(rownames(datTraits)==rownames(datExprFemale))
# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(t(datExprFemale),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits),"C",sep="")
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")
# Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
# the following 2 lines differ from what is written in the book
datExprFemale=datExprFemale[!remove.samples,]
datTraits=datTraits[!remove.samples,]
# Recompute the sample network among the remaining samples
A=adjacency(t(datExprFemale),type="distance")
# Let's recompute the Z.k values of outlyingness
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)

#############################################################
#Section 12.2: Co-expression modules in female mouse livers 
#12.2.1: Choosing the soft threshold beta via scale free topology
# Choose a set of soft thresholding powers
powers=c(1:10) # in practice this should include powers up to 20.
# choose power based on SFT criterion
sft=pickSoftThreshold(datExprFemale,powerVector=powers)

#Digression: if you want to pick a soft threshold for a signed network write
#sft=pickSoftThreshold(datExprFemale,powerVector=powers, networkType = "signed")
# but then you should consider higher powers. Default beta=12.

# Plot the results:
par(mfrow=c(1,2))
# SFT index as a function of different powers
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red")
#############################################################
#12.2.3 Automatic module detection via dynamic tree cutting
mergingThresh = 0.25
net = blockwiseModules(datExprFemale,corType="pearson",
                       maxBlockSize=5000,networkType="unsigned",power=7,minModuleSize=30,
                       mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
                       pamRespectsDendro=FALSE,saveTOMFileBase="femaleMouseTOM")
moduleLabelsAutomatic=net$colors
# Convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)

# A data frame with module eigengenes can be obtained as follows
MEsAutomatic=net$MEs

#this is the body weight
weight = as.data.frame(datTraits$weight_g)
names(weight)="weight"
# Next use this trait to define a gene significance variable
GS.weight=as.numeric(cor(datExprFemale,weight,use="p"))
# This translates the numeric values into colors
GS.weightColor=numbers2colors(GS.weight,signed=T)
blocknumber=1
datColors=data.frame(moduleColorsAutomatic,GS.weightColor)[net$blockGenes[[blocknumber]],]

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[blocknumber]],colors=datColors,
                    groupLabels=c("Module colors","GS.weight"),dendroLabels=FALSE,
                    hang=0.03,addGuide=TRUE,guideHang=0.05)

#############################################################
#12.2.3 Blockwise module detection for large networks
bwnet = blockwiseModules(datExprFemale,corType="pearson",
                         maxBlockSize=2000,networkType="unsigned",power=7,minModuleSize=30,
                         mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
                         pamRespectsDendro=FALSE,saveTOMFileBase="femaleMouseTOM-blockwise",verbose=3)
# Relabel blockwise modules so that their labels
# match those from our previous analysis
moduleLabelsBlockwise=matchLabels(bwnet$colors,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsBlockwise=labels2colors(moduleLabelsBlockwise)

# measure agreement with single block automatic procedure
mean(moduleLabelsBlockwise==moduleLabelsAutomatic) 
blockNumber=2
# Plot the dendrogram for the chosen block
plotDendroAndColors(bwnet$dendrograms[[blockNumber]],
                    moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],"Module colors",
                    main=paste("Dendrogram and module colors in block",blockNumber),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

#############################################################
#12.2.4 Manual, stepwise module detection
# We now calculate the weighted adjacency matrix, using the power 6:
A = adjacency(datExprFemale, power = 7)
# Digression: to define a signed network choose
#A = adjacency(datExprFemale, power = 12, type="signed")

#define a dissimilarity based on the topological overlap
dissTOM =TOMdist(A)
#hierarchical clustering
geneTree = flashClust(as.dist(dissTOM),method="average")
# here we define the modules by cutting branches
moduleLabelsManual1=cutreeDynamic(dendro=geneTree,distM=dissTOM,
                                  method="hybrid",deepSplit=2,pamRespectsDendro=F,minClusterSize=30)

# Relabel the manual modules so that their labels
# match those from our previous analysis
moduleLabelsManual2=
  matchLabels(moduleLabelsManual1,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsManual2=labels2colors(moduleLabelsManual2)

# Calculate eigengenes
MEList=moduleEigengenes(datExprFemale,colors=moduleColorsManual2)
MEs = MEList$eigengenes
# Add the weight to existing module eigengenes
MET=orderMEs(cbind(MEs,weight))
# Plot the relationships among the eigengenes and the trait
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

# automatically merge highly correlated modules
merge=mergeCloseModules(datExprFemale,moduleColorsManual2,
                        cutHeight=mergingThresh)
# resulting merged module colors
moduleColorsManual3 = merge$colors
# eigengenes of the newly merged modules:
MEsManual = merge$newMEs

# Show the effect of module merging by plotting the
# original and merged module colors below the tree
datColors=data.frame(moduleColorsManual3,moduleColorsAutomatic, moduleColorsBlockwise,GS.weightColor)
plotDendroAndColors(geneTree,colors=datColors,
                    groupLabels=c("manual hybrid","single block","2 block","GS.weight"),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
# check the agreement between manual and automatic module labels
mean(moduleColorsManual3==moduleColorsAutomatic)

#############################################################
#12.2.5  Relating modules to physiological traits
# Choose a module assignment
moduleColorsFemale=moduleColorsAutomatic
# Define numbers of genes and samples
nGenes = ncol(datExprFemale)
nSamples = nrow(datExprFemale)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprFemale,moduleColorsFemale)$eigengenes
MEsFemale = orderMEs(MEs0)
modTraitCor = cor(MEsFemale, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
#Since we have a moderately large number of modules and traits,
#a suitable graphical representation will help in reading
#the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",
                   signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits),
               yLabels = names(MEsFemale), ySymbols = names(MEsFemale), 
               colorLabels =FALSE,colors=greenWhiteRed(50),textMatrix=textMatrix,
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# calculate the module membership values
# (aka. module eigengene based connectivity kME):
datKME=signedKME(datExprFemale, MEsFemale)

colorOfColumn=substring(names(datKME),4)
par(mfrow = c(2,2))
selectModules=c("blue","brown","black","grey")
par(mfrow=c(2,length(selectModules)/2))
for (module in selectModules) {
  column = match(module,colorOfColumn)
  restModule=moduleColorsFemale==module
  verboseScatterplot(datKME[restModule,column],GS.weight[restModule],
                     xlab=paste("Module Membership ",module,"module"),ylab="GS.weight",
                     main=paste("kME.",module,"vs. GS"),col=module)}

#############################################################
#12.2.6 Output file for gene ontology analysis
# Read in the probe annotation
GeneAnnotation=read.csv(file="GeneAnnotation.csv")
# Match probes in the data set to those of the annotation file
probes = names(datExprFemale)
probes2annot = match(probes,GeneAnnotation$substanceBXH)
# data frame with gene significances (cor with the traits)
datGS.Traits=data.frame(cor(datExprFemale,datTraits,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datOutput=data.frame(ProbeID=names(datExprFemale),
                     GeneAnnotation[probes2annot,],moduleColorsFemale,datKME,datGS.Traits)
# save the results in a comma delimited file
write.table(datOutput,"FemaleMouseResults.csv",row.names=F,sep=",")

#Copy and paste the following R code from section 11.8.1 into the R session.
library(sem)
# Now we specify the 5 single anchor models
# Single anchor model 1: M1->A->B
CausalModel1=specifyModel()
A -> B, betaAtoB, NA
M1 -> A, gammaM1toA, NA
A <-> A, sigmaAA, NA
B <-> B, sigmaBB, NA
M1 <-> M1, sigmaM1M1,NA

# Single anchor model 2: M->B->A
CausalModel2=specifyModel()
B -> A, betaBtoA, NA
M1 -> B, gammaM1toB, NA
A <-> A, sigmaAA, NA
B <-> B, sigmaBB, NA
M1 <-> M1, sigmaM1M1,NA

# Single anchor model 3:  A<-M->B
CausalModel3=specifyModel()
M1 -> A, gammaM1toA, NA
M1 -> B, gammaM1toB, NA
A <-> A, sigmaAA, NA
B <-> B, sigmaBB, NA
M1 <-> M1, sigmaM1M1,NA

# Single anchor model 4:  M->A<-B
CausalModel4=specifyModel()
M1 -> A, gammaM1toA, NA
B -> A, gammaBtoA, NA
A <-> A, sigmaAA, NA
B <-> B, sigmaBB, NA
M1 <-> M1, sigmaM1M1,NA

# Single anchor model 5:  M->B<-A
CausalModel5=specifyModel()
M1 -> B, gammaM1toB, NA
A -> B, gammaAtoB, NA
A <-> A, sigmaAA, NA
B <-> B, sigmaBB, NA
M1 <-> M1, sigmaM1M1,NA

# Function for obtaining the model fitting p-value
# from an sem object:
ModelFittingPvalue=function(semObject){
  ModelChisq=summary(semObject)$chisq
  df= summary(semObject)$df
  ModelFittingPvalue=1-pchisq(ModelChisq,df)
  ModelFittingPvalue}

# this function calculates the single anchor score
LEO.SingleAnchor=function(M1,A,B){
  datObsVariables= data.frame(M1,A,B)
  # this is the observed correlation matrix
  S.SingleAnchor=cor(datObsVariables,use="p")
  m=dim(datObsVariables)[[1]]
  
  semModel1 =sem(CausalModel1,S=S.SingleAnchor,N=m)
  semModel2 =sem(CausalModel2,S=S.SingleAnchor,N=m)
  semModel3 =sem(CausalModel3,S=S.SingleAnchor,N=m)
  semModel4 =sem(CausalModel4,S=S.SingleAnchor,N=m)
  semModel5 =sem(CausalModel5,S=S.SingleAnchor,N=m)
  
  # Model fitting p-values for each model
  p1=ModelFittingPvalue(semModel1)
  p2=ModelFittingPvalue(semModel2)
  p3=ModelFittingPvalue(semModel3)
  p4=ModelFittingPvalue(semModel4)
  p5=ModelFittingPvalue(semModel5)
  LEO.SingleAnchor=log10(p1/max(p2,p3,p4,p5))
  
  data.frame(LEO.SingleAnchor,p1,p2,p3,p4,p5)
} # end of function

# Get SNP markers which encode module QTLs 
SNPdataFemale=read.csv("SNPandModuleQTLFemaleMice.csv")
# find matching row numbers to line up the SNP data
snpRows=match(dimnames(datExprFemale)[[1]],SNPdataFemale$Mice)
# define a data frame whose rows correspond to those of datExpr
datSNP = SNPdataFemale[snpRows,-1]
rownames(datSNP)=SNPdataFemale$Mice[snpRows]
# show that row names agree
table(rownames(datSNP)==rownames(datExprFemale))

SNP=as.numeric(datSNP$mQTL19.047)
weight=as.numeric(datTraits$weight_g)
MEblue=MEsFemale$MEblue
MEbrown=MEsFemale$MEbrown

# evaluate the relative causal fit
# of SNP -> ME -> weight
LEO.SingleAnchor(M1=SNP,A=MEblue,B=weight)[[1]]
LEO.SingleAnchor(M1=SNP,A=MEbrown,B=weight)[[1]]

whichmodule="blue"
restGenes=moduleColorsFemale==whichmodule
datExprModule=datExprFemale[,restGenes]
attach(datExprModule)
LEO.SingleAnchorWeight=rep(NA, dim(datExprModule)[[2]] )
for (i in c(1:dim(datExprModule)[[2]]) ){
  printFlush(i)
  LEO.SingleAnchorWeight[i]=
    LEO.SingleAnchor(M1=SNP,A=datExprModule[,i],B=weight)[[1]]}

# Find gens with a very significant LEO score
restCausal=LEO.SingleAnchorWeight>3 # could be lowered to 1
names(datExprModule)[restCausal]

# this provides us the corresponding gene symbols
data.frame(datOutput[restGenes,c(1,6)], LEO.SingleAnchorWeight= LEO.SingleAnchorWeight )[restCausal,]

#############################################################
#12.4 Visualizing the network
# WARNING: On some computers, this code can take a while to run (20 minutes??). 
# I suggest you skip it.
# Set the diagonal of the TOM disscimilarity to NA
diag(dissTOM) = NA
# Transform dissTOM with a power to enhance visibility
TOMplot(dissim=dissTOM^7,dendro=geneTree,colors=moduleColorsFemale, main = "Network heatmap plot, all genes")
cmd1=cmdscale(as.dist(dissTOM),2)
par(mfrow=c(1,1))
plot(cmd1,col=moduleColorsFemale,main="MDS plot", xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")

#############################################################
#VisANT plot and software
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExprFemale, power=7)
module = "blue"
# Select module probes
probes = names(datExprFemale)
inModule = (moduleColorsFemale==module)
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
# Because the module is rather large,
# we focus on the 30 top intramodular hub genes
nTopHubs = 30
# intramodular connectivity
kIN = softConnectivity(datExprFemale[, modProbes])
selectHubs = (rank (-kIN) <= nTopHubs)
vis = exportNetworkToVisANT(modTOM[selectHubs,selectHubs],
                            file=paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted=TRUE,threshold = 0, probeToGene=
                              data.frame(GeneAnnotation$substanceBXH,GeneAnnotation$gene_symbol))

#############################################################
#Cytoscape and Pajek software
# select modules
modules = c("blue","brown")
# Select module probes
inModule=is.finite(match(moduleColorsFemale,modules))
modProbes=probes[inModule]
match1=match(modProbes,GeneAnnotation$substanceBXH)
modGenes=GeneAnnotation$gene_symbol[match1]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
                               nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
                               weighted = TRUE, threshold = 0.02,nodeNames=modProbes,
                               altNodeNames = modGenes, nodeAttr = moduleColorsFemale[inModule])

#############################################################
#12.5 Module preservation between female and male mice
# WARNING: On some computers, this code can take a while to run (30 minutes??). 
# I suggest you skip it.
maleData = read.csv("LiverMaleFromLiverFemale3600.csv")
names(maleData)
datExprMale=data.frame(t(maleData[,-c(1:8)]))
names(datExprMale)=maleData$substanceBXH

# We now set up the multi-set expression data
# and corresponding module colors:
setLabels = c("Female", "Male")
multiExpr=list(Female=list(data=datExprFemale),
               Male=list(data=datExprMale))
moduleColorsFemale=moduleColorsAutomatic
multiColor=list(Female=moduleColorsFemale)

# The number of permutations drives the computation time
# of the module preservation function. For a publication use 200 permutations.
# But for brevity, let's use a small number
nPermutations1=20
# Set it to a low number (e.g. 3) if only the medianRank statistic
# and other observed statistics are needed.
# Permutations are only needed for calculating Zsummary
# and other permutation test statistics.
# set the random seed of the permutation test analysis
set.seed(1)
system.time({
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1, nPermutations = nPermutations1,
                          randomSeed = 1, quickCor = 0, verbose = 3)
})
# Save the results of the module preservation analysis
save(mp, file = "modulePreservation.RData")
# If needed, reload the data:
load(file = "modulePreservation.RData")

# specify the reference and the test networks
ref=1; test = 2

Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats

# Z statistics from the permutation test analysis
Z.PreservationStats

# Let us now visualize the data.
modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modColors[selectModules]

#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres

par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules],pch = 21,main="medianRank Preservation",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)

# plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Zsummary Preservation",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)

whichmodule="lightgreen"
Eigengene=MEsFemale$MElightgreen
datExprModule=datExprFemale[,moduleColorsFemale==whichmodule]
# set the margins of the graphics window
par(mfrow=c(1,1),mar=c(0.3, 5.5, 3, 2))
# create a heatmap whose columns correspond to the arrays
# and whose rows correspond to genes
plotMat(t(scale(datExprModule)),cex.axis=2,nrgcols=30,rlabels=F,
        rcols=whichmodule,main=paste("heatmap",whichmodule,"module"))

#scatter plot between eigengene and sample network connectivity
par(mfrow=c(1,1))
verboseScatterplot(Eigengene,Z.k,xlab=paste("ME",whichmodule,sep=""))
abline(h=-2,col="red",lwd=2)

# This calculates the quality of the modules in the female data
Obs.QualityStatsFemale= mp$quality$observed[[1]][[2]]
Z.QualityStatsFemale=mp$quality$Z[[1]][[2]]
QualityStats=data.frame(Obs.QualityStatsFemale[,1:2],
                        Zsummary.qual=Z.QualityStatsFemale$Zsummary.qual)
QualityStats["lightgreen",]

#############################################################
#12.6 Consensus modules between female and male liver tissues
# number of networks used in the consensus network analysis:
nSets=2
# Vector with descriptive names of the two sets.
setLabels=c("Female liver", "Male liver")
shortLabels=c("Female", "Male")
#Define a list whose components contain the data
multiExpr=vector(mode="list",length=nSets)
multiExpr[[1]] = list(data = datExprFemale)
names(multiExpr[[1]]$data) = names(datExprFemale)
rownames(multiExpr[[1]]$data) = dimnames(datExprFemale)[[1]]
multiExpr[[2]] = list(data = datExprMale)
names(multiExpr[[2]]$data) = names(datExprMale)
rownames(multiExpr[[2]]$data) = dimnames(datExprMale)[[1]]
# Check that the data has the correct format:
exprSize = checkSets(multiExpr)

# The variable exprSize contains useful information
# about the sizes of all of the data sets
# now we run automatic module detection procedure
netConsensus=blockwiseConsensusModules(multiExpr,maxBlockSize=5000,
                                       power=7,minModuleSize=30,deepSplit=2,pamRespectsDendro=FALSE,
                                       mergeCutHeight=0.25,numericLabels=TRUE,
                                       minKMEtoStay=0,saveTOMs=TRUE)
# list of consensus module eigengenes
consMEs = netConsensus$multiMEs
# module labels
modLabCons0 = netConsensus$colors
# Relabel the consensus modules so that their labels match those
# from the automatic analysis in female mice only
modLabCons = matchLabels(modLabCons0,moduleLabelsAutomatic)
# check the agreement between consensus- and females-only modules
mean(modLabCons==moduleLabelsAutomatic) 

# Convert the numeric labels to color labels
moduleColorsConsensus = labels2colors(modLabCons)

blocknumber=1
consTree = netConsensus$dendrograms[[blocknumber]]
# plot the consensus dendrogram and module colors
datColors=data.frame(moduleColorsConsensus,moduleColorsFemale)[netConsensus$blockGenes[[blocknumber]],]
plotDendroAndColors(consTree,datColors,c("Cons module","female module"),dendroLabels=FALSE,
                    hang=0.03,addGuide=TRUE,guideHang=0.05,main="Consensus gene dendrogram and module colors")

#############################################################
#12.6.1 Relating consensus modules to the traits


