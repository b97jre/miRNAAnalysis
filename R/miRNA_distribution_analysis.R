#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)




library(ggplot2)
source('/Users/johanreimegard/git/miRNAAnalysis/R/miRNAAnalysisFunctions.R')


setwd("/Users/johanreimegard/Vetenskap/Data/Mouse_miRNAs/HTseqCount")


metaData = read.table("metaDataTable.txt", sep = "\t", header = TRUE )
metaData$PhenotypeColor = "Blue"
metaData$PhenotypeColor[ metaData$Phenotype =="wt" ] = "red"
metaData$DateSymbol = 1
metaData$DateSymbol[ metaData$Date =="2nd" ] = 2

metaDatamiRNA = metaData[ metaData$Reference == "pre-miRNA",  ]




miRNAs = read.table("Name_Stranded_18_24.htseq.count.table.txt", sep = "\t", header =  TRUE)
miRNAsBoth = read.table("htseq_id_Stranded_18_24.htseq.count.table.txt", sep = "\t", header =  TRUE)
premiRNAsStranded = read.table("premiRNA_Stranded_18_24.htseq.count.table.txt", sep = "\t", header =  TRUE)

premiRNAsNotStranded = read.table("premiRNA_Unstranded_18_24.htseq.count.table.txt", sep = "\t", header =  TRUE)




class(miRNAs)
head(premiRNAsStranded)
Allinfo =as.data.frame(rowSums(miRNAsBoth))  
Allinfo$premiRNAs = rowSums(premiRNAsStranded)
Allinfo$All = rowSums(premiRNAsNotStranded)

colnames( Allinfo ) <- c("miRNA" , "premiRNA" , "All")
dim(Allinfo)
head(Allinfo)


Allinfo = Allinfo[1:1187, ]

tail(Allinfo)
Allinfo$wrongStrand = Allinfo$All-Allinfo$premiRNA
Allinfo$notmiRNA = Allinfo$premiRNA-Allinfo$miRNA

nomiRNAs = Allinfo[order(Allinfo$notmiRNA,decreasing = TRUE),]
nomiRNAs = nomiRNAs[nomiRNAs$All>0 &(nomiRNAs$notmiRNA/nomiRNAs$All) > 0.25, ]

reverseStrand = Allinfo[order(Allinfo$wrongStrand,decreasing = TRUE),]
reverseStrand = reverseStrand[reverseStrand$All>0 &(reverseStrand$wrongStrand/reverseStrand$All) > 0.25, ]

inBoth = intersect(row.names(nomiRNAs), row.names(reverseStrand))

Allinfo$quality = "GOOD"
dim(Allinfo)
Allinfo$quality[row.names(Allinfo) %in% row.names(nomiRNAs)] = "ANNOTATION"
Allinfo$quality[row.names(Allinfo) %in% row.names(reverseStrand)] = "STRAND"
Allinfo$quality[row.names(Allinfo) %in% inBoth] = "BOTH"

Allinfo$quality[ Allinfo$All <36] = "LOW"

Allinfo$correct = 1 - (Allinfo$All- Allinfo$miRNA)/Allinfo$All
Allinfo$correct[is.nan(Allinfo$correct)] = 0


m <- ggplot(Allinfo, aes(x=correct))
m + geom_density(adjust=1/5)

ggsave("AllDistribution.pdf")

m <- ggplot(Allinfo, aes(x=correct, col = quality))
m + geom_density(adjust=1/5)+ facet_wrap(~ quality, ncol = 3,, scales = "free_y")

ggsave("AllDistribution_Facet.pdf")

c <- ggplot(Allinfo, aes(factor(quality), fill = quality))

c + geom_bar()

ggsave("miRNAreadsDistribution.pdf")


head(exp.data)
exp.data = miRNAsBoth[rownames(miRNAsBoth) %in% goodmiRNAs,metaDatamiRNA$colnames ]
colnames(exp.data) <- metaDatamiRNA$Name
exp.data <- exp.data + 1
lib.size <- apply(exp.data,2,sum)
scale.factors <- calcNormFactors(exp.data, method="TMM")

norm.data <- t(t(exp.data)/(scale.factors*lib.size))
norm.data <- log(norm.data)

mir.pca <- prcomp(t(norm.data),) 



pdf(file = "PCA_1_2.pdf")
plot(mir.pca$x[,1], mir.pca$x[,2],pch=metaDatamiRNA$DateSymbol,col= metaDatamiRNA$PhenotypeColor)  ## plot  PC1 and PC2
text(mir.pca$x[,1], mir.pca$x[,2], rownames(mir.pca$x), cex=0.7, pos=1, col= metaDatamiRNA$PhenotypeColor) 
dev.off()

pdf(file = "PCA_3_2.pdf")
plot(mir.pca$x[,3], mir.pca$x[,2],pch=metaDatamiRNA$DateSymbol,col= metaDatamiRNA$PhenotypeColor)  ## plot  PC1 and PC2
text(mir.pca$x[,3], mir.pca$x[,2], rownames(mir.pca$x), cex=0.7, pos=1, col= metaDatamiRNA$PhenotypeColor) 
dev.off()




norm.data[head(rownames(as.data.frame(sort(mir.pca$rotation[,2], decreasing=TRUE)))),]
norm.data[head(rownames(as.data.frame(sort(mir.pca$rotation[,2], decreasing=FALSE)))),]

mir.pca$rotation[,1]

mir.pca <- prcomp(t(exp.data)) 
mir.pca$sdev


dev.off()
pdf(file = "HeatMap.pdf")
heatmap(as.matrix(exp.data), scale="none", cexCol=0.7, cexRow=0.2)
dev.off()

class(mir.pca)
mir.pca$
