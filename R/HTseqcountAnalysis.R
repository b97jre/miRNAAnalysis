

#load.package("ggplot2")

library(ggplot2)
source('/Users/johanreimegard/git/miRNAAnalysis/R/miRNAAnalysisFunctions.R')




setwd("/Users/johanreimegard/Vetenskap/Data/Mouse_miRNAs/HTseqCount/test")

metaData = read.table("metaDataTable.txt", sep = "\t", header = TRUE )
metaData$PhenotypeColor = "Blue"
metaData$PhenotypeColor[ metaData$Phenotype =="wt" ] = "red"
metaData$DateSymbol = 1
metaData$DateSymbol[ metaData$Date =="2nd" ] = 2
metaData$lty = 1
metaData$lty[ metaData$Date =="2nd" ] = 3




metaDataPremiRNA = metaData[ metaData$Reference =="pre-miRNA",  ] 

length121 = read.table("up_121.lengthDistribution.tab.txt", sep = "\t", header = TRUE)
length154 = read.table("up_154.lengthDistribution.tab.txt", sep = "\t", header = TRUE)
AllLength = cbind(length121,length154)
AllLength = AllLength[,order(colnames(AllLength))]

lty = c(1:length(colnames(length121)),1:length(colnames(length154)))
col = c(rep("red",length(colnames(length121))),rep("blue", times = length(colnames(length154))))

start = 5
stop = 30
colnames(Alllength) <- metaDataPremiRNA$Name

pdf(paste ("miRNA_mapping",start,stop,"pdf", sep = "."))
plotLengthDistribution(Alllength, startLength = start, stopLength = stop, lty = metaDataPremiRNA$lty, col = metaDataPremiRNA$PhenotypeColor )
dev.off()

start = 18
stop = 25
pdf(paste ("miRNA_mapping",start,stop,"pdf", sep = "."))
plotLengthDistribution(AllLength, startLength = start, stopLength = stop, lty = metaDataPremiRNA$lty, col = metaDataPremiRNA$PhenotypeColor )
dev.off()




head(AllLength )

tempLength = AllLength[15:40,]
myPca <- prcomp(t(tempLength))

pdf('pca_plot.pdf')
plot(myPca$x[,1],myPca$x[,2],col=metaDataPremiRNA$PhenotypeColor,pch=metaDataPremiRNA$DateSymbol)
text(myPca$x[,1], myPca$x[,4], labels=metaDataPremiRNA$Name, cex= 0.7)

dev.off()


rownames(AllLength)
length

dev.off()






OneTwoOne = read.table("121.htseq.count.table.txt", sep = "\t")
colnames(OneTwoOne) <- paste("121",colnames(OneTwoOne),sep = "_")
OneTwoOne = OneTwoOne[1:1187, ]


test = read.table("test.txt", sep = "\t")[1:1187, ]
test = test[order(test$V1),]

dim(test)

test2 = read.table("test2.txt", sep = "\t")
test2 = test2[order(test2$V1),]
test4 = test2[!test2$V1 %in% test$V1,]
dim(test3)
head(test3)
test3$HT = test$V2
test3$diff = test3$V3 -test3$HT

test3 = test3[order(test3$diff),]


colnames(test3) <- c("Name","length","idxstats","notsure","HTseq","diff")

write.table(test3, "diff.txt", sep = "\t", quote = FALSE,row.names = FALSE,col.names = TRUE)


OneFiveFour = read.table("154.htseq.count.table.txt", sep = "\t", header = TRUE)
colnames(OneFiveFour) <- paste("154",colnames(OneFiveFour),sep = "_")
OneFiveFour = OneFiveFour[1:1187, ]


AllSamples = cbind(OneFiveFour,OneTwoOne)
AllSamples = AllSamples[,order(colnames(AllSamples))]

AllSamplesNorm = sweep((AllSamples+1), 2, colSums((AllSamples+1)), FUN="/")

head(AllSamplesNorm)
myPca <- prcomp(t(AllSamplesNorm))

pdf('pca_plot.pdf')
plot(myPca$x[,1],myPca$x[,4],col=metaDataPremiRNA$PhenotypeColor,pch=metaDataPremiRNA$DateSymbol)
text(myPca$x[,1], myPca$x[,4], labels=metaDataPremiRNA$Name, cex= 0.7)

dev.off()
dim(OneFiveFour)







