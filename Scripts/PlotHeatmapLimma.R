suppressPackageStartupMessages(require(gplots))
suppressPackageStartupMessages(library(RColorBrewer, quiet=TRUE))
library(limma)
library(edgeR)

rnaData= commandArgs()[7]
#testInFiles= commandArgs()[7]
classFilePath = commandArgs()[8]
#logTheData = commandArgs()[9] == "TRUE"
outFiletxt = commandArgs () [9]
#histS = commandArgs() [10]
outFilePath = commandArgs()[10]


#Reads in the RNA Seq Data
data = t(data.matrix(read.table(rnaData, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1, quote="\"", check.names=FALSE)))
print ("dim of RNA data")
#dim(data)
head(data)

#testData = replicate(360, rnorm(21, mean=100))
#dimnames(testData) <- list(rownames(testData, do.NULL =FALSE, prefix = "row"), colnames(testData,do.NULL =FALSE, prefix = "col"))
#write.table(testData, file="testData.txt" ,sep = "\t", quote=F, col.names=NA)
#testData = (data.matrix(read.table(testInFiles, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1, quote="\"", check.names=FALSE)))


classData = read.table(classFilePath, sep="\t", header=FALSE, check.names=FALSE, row.names=NULL, stringsAsFactors=FALSE)
classData = classData[classData[,1] %in% colnames(data),]
classData = data.frame(classData)
classData = classData[colnames(data),] # make sure classData rows are in same order as expression columns


#dge = DGEList(counts=data)
#head(dge)
#dge =  calcNormFactors(dge)

#Group = factor(classData$V2, levels=c("Serous","NonSerous"))
#design =  model.matrix(~0+Group)
#print ("dim of design matrix")
#colnames(design) = c("Serous","NonSerous")
#rownames(design) = classData$V1
#print(design)

#contrast.matrix = makeContrasts(Serous-NonSerous, levels=design)
#print(contrast.matrix)

#v <- voom(data,design,plot=TRUE,normalize="quantile")

#v <- voom(data,design,plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit,coef=ncol(design))
#Determines the foldchange
fit <- treat(fit,lfc=log2(1.2))
topTreat(fit,coef=ncol(design))

#write.table(fit, outFiletxt ,sep = "\t", col.names=NA)
#Seperates the Class Samples into Seperate Variables
class1 = sort(unique(classData[,2]))[1]
class2 = sort(unique(classData[,2]))[2]
class1Samples = classData[which(classData[,2]==class1),1]
class2Samples = classData[which(classData[,2]==class2),1]
print(class1Samples)
print(class2Samples)

NonSerous = testData[,class1Samples]

DUSP3_NS = (NonSerous[1,]) #all the NS gene counts for each sample
DUSP3_S = (Serous[1,])


pValue = t.test(DUSP3_NS, DUSP3_S)
print(pValue)

Averages all of the P-values for Serous vs. NonSerous
NonSerousMean = apply((NonSerous), 1, mean)
SerousMean = apply((Serous), 1, mean)

tTest= t.test(NonSerousMean,SerousMean)
tTest
stop()



pdf(histNS)
maxValue=max(c(DUSP3_NS, DUSP3_S))
plot(density(DUSP3_NS),col= "blue", xlim =c(0,10000),main= "Serous vs Non-Serous RNA-seq Values(RPS6KA1, p=0)")
lines(density(Serous), col="red")
legend("topright", inset=.05,
   c("NonSerous","Serous"),fill= c("blue", "red"), horiz=TRUE)
dev.off()

pdf(histS)
plot(SerousMean)
dev.off()


print (head(data1[,1]))
print("__________")
print (head(data2[,1]))

geneTTestPValues = NULL
print(data1[1:5,1:5])
print(data2[1:5,1:5])


for (i in 1:nrow(data1))
{
 pValue = t.test(data1[i,], data2[i,])$p.value
  geneTTestPValues = c(geneTTestPValues, pValue)
}

#print(geneTTestPValues)
print(dim(data1))
data1 = data1[which(geneTTestPValues<=0.0000000000000000000000000000001),]
data2 = data2[which(geneTTestPValues<=0.0000000000000000000000000000001),]
print(dim(data1))

#Averages all of the P-values for Serous vs. NonSerous
data1Means = apply(v1, 1, mean)
data2Means = apply(v2, 1, mean)


#Brings together the two means into one variabe?
dataMeans = rbind(data1Means, data2Means)

#Assigns Serous or Non Serous to the Columns
rownames(dataMeans) = c(class1, class2)

colors = brewer.pal(11,"RdGy")
colors = brewer.pal(11,"RdBu")
colors = colors[11:1]
colors = colors[round(quantile(1:length(colors), probs=seq(0, 1, 0.20)))]

key = TRUE
cexRow = 1
cexCol= 1.5
pdf(outFilePath, width=5, height=10)

#Creates the Heat Map
heatmap.2(t(dataMeans), Rowv=TRUE, col=colors, dendrogram="row", trace="none", key=TRUE,  symkey=FALSE,, density.info="none", scale="none", main="", cexRow=cexRow,cexCol=cexCol,  labRow=colnames(dataMeans), labCol=rownames(dataMeans),
margin=c(6,16), lmat=rbind(c(4,3),c(2,1)), lwid=c(2,4), lhei=c(2.0,14))

margin=FALSE, lmat=rbind(c(0,3),c(2,1),c(0,4)), lwid=c(1.5,4), lhei=c(1.5,4,1),srtCol=180)


