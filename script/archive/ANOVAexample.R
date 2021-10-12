#ANOVAexample
### CREATE TABLE
Fruit<-c("Apple","Apple","Orange","Orange","Banana","Banana","Banana","Avocado","Avocado")
Weight<-c(60,50,80,72,40,45,85,90,95)
horDiam<-c(60,70,75,70,30,35,50,60,70)
Price<-c(5,8,7,9,3,4,8,13,14)
Dummy<-c(4,6,2,8,1,2,3,8,6)
myData<-data.frame(Fruit=Fruit, Weight=Weight, horDiam=horDiam, Price=Price, Dummy=Dummy)
rownames(myData)<-c("Apple1","Apple2","Orange1","Orange2","Banana1","Banana2","Banana3","Avocado1","Avocado2")


### ANOVA
fit.aov<-list()
summaryAOV<-list()
for (i in 1:3){
  fit.aov[[i]]<-aov(myData[,i+1]~myData[,1])
  summaryAOV[[i]]<-summary(fit.aov[[i]])
}

### TUKEY
par(mfrow=c(1,3))
testTukey<-list()
mainTukey<-c("Weight", "Horiz. Diameter", "Price")
for (i in 1:3){
  testTukey[[i]]<-TukeyHSD(fit.aov[[i]], conf.level = 0.95)
  plot(testTukey[[i]], main=mainTukey[i])
}

### CLUSTERING
plot( hclust(dist(myData), method="ward") )

### CLUSTERING WITH P-VALUE
#install.packages("pvclust")
#install.packages("pvrect") #not available (for R version 3.6.2)
#fit <- pvclust(t(myData[,-1]), method.hclust="ward", method.dist="euclidean")
#plot(fit)
#pvrect(fit, alpha=0.95)

### analysis of multivariate homogeneity of group dispersions
#install.packages("vegan")
require(vegan)
distance<-vegdist(myData[,2:5], method="euclidean")
model<-betadisper(distance, myData[,1])
permutest(model, pairwise = TRUE)


#install.packages("ade4")
#install.packages("randtest.discrimin") #not available (for R version 3.6.2)
#require(ade4)
#discr <- discrimin(dudi.pca(myData[,2:5], scan = FALSE), myData[,1], scan = FALSE)
#randtest.discrimin(discr)


