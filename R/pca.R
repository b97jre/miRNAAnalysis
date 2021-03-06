#install.packages("GGally")

library(plyr)
library(ggplot2)
library(GGally)



gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

gg_shapeValues  <- function(n){
  shapes = c(16,17,15,3,7,8)
  subshapes = rep(shapes,len = n)
}


runPC <- function(prcompResult, metaDataTable, n.comp = 5,scale = FALSE, center = TRUE,
                  sampleColName = "Name", colorComponent = NULL , pchComponent = NULL){
  
  PCinfo = list()
  pca = prcomp(t(rpkm), scale=scale, center=center)
  
  
  
  
}
  




# Create PCAplot requires that the colnames in the expressionTable are present in the metaDataTable[sampleColName, ] 
createMultipePCAplot <- function(prcompResult, metaDataTable, plot.comp = 1:5,
                                 sampleColName = "Name", colorComponent = NULL , pchComponent = NULL)
  
  # Nr of components
  plot.comp 
  n.comp = length(plot.comp)

  #Plot info
  pca.basis = prcompResult[['x']]
  pca.basis = pca.basis[order(rownames(pca.basis) ),]
  
  # PC componenet information info 
  e.var = (prcompResult[['sdev']]^2 / sum(prcompResult[['sdev']]^2))[plot.comp]
  e.var = as.data.frame( e.var)
  
  
  e.var$components = colnames(pca.basis)[plot.comp]
  e.var$info = sprintf('%s: %.1f%s',e.var$components,e.var$e.var*100, "%")
  
  
  
  
  
  
  #Table info 
  metaInfo = metaDatatable[metaDatatable[[sampleColName]] %in% rownames(pca.basis),  ]
  rownames(metaInfo) <- metaInfo[[sampleColName]]
  metaInfo = metaInfo[order(rownames(metaInfo)),]


  ggInfo = as.data.frame(pca.basis[, plot.comp])
  ggInfo[pchComponent] = metaInfo[pchComponent]
  ggInfo[colorComponent] = metaInfo[colorComponent]
  ggInfo[paste(pchComponent, colorComponent, sep = "_")] =as.factor(paste(ggInfo[[pchComponent]], ggInfo[[colorComponent]], sep = "_"))

  breakInfo =  as.data.frame(levels(ggInfo[[paste(pchComponent, colorComponent, sep = "_")]]))
  colnames(breakInfo) <- "Both"
  ll <- unlist(strsplit(as.character(breakInfo$Both), "_"))
  breakInfo$shapeLevels = as.factor(ll[seq(from = 1, to = length(ll)-1, by = 2)])
  breakInfo$colourLevels = as.factor(ll[seq(from = 2, to = length(ll), by = 2)])
  breakInfo$colourValues = "black"
  breakInfo$shapeValues = 16
  
  
  df = data.frame(colLevels = levels(breakInfo$colourLevels), colValues = gg_color_hue(length(levels(breakInfo$colourLevels))))
  for( i in 1:length(df$colLevels)){
    breakInfo$colourValues[which(breakInfo$colourLevels == df$colLevels[i])] = as.character(df$colValues[[i]])
  }
  
  df = data.frame(shapeLevels = levels(breakInfo$shapeLevels), shapeValues = gg_shapeValues(length(levels(breakInfo$shapeLevels))))
  for( i in 1:length(df$shapeLevels)){
    breakInfo$shapeValues[which(breakInfo$shapeLevels == df$shapeLevels[i])] = as.integer(df$shapeValues[[i]])
  }
  
   pcaPlot <- ggpairs(ggInfo,columns = c(1:(n.comp+1)) , params = c(adjust = 0.5),  colour=colorComponent,  shape = pchComponent, linetype = pchComponent,    
        lower = list(continuous = "points", combo = "blank", discrete = "blank"),
        upper = list(continuous = "blank", combo = 'box', discrete = "blank")
  )

   
 
   plot <- ggplot2::ggplot()
   plot <- plot +
     ggplot2::scale_y_discrete(name="",  limits = c("Line","Shape","Colour", "Head")) +
     ggplot2::scale_x_continuous(limits=c(0,1), breaks=NULL, name="") + 
     ggplot2::geom_text(data=NULL, mapping=aes(y='Line', x=0.05,label = pchComponent,hjust=0), size=4)+
     ggplot2::geom_text(data=NULL, mapping=aes(y='Shape', x=0.05,label = pchComponent,hjust=0), size=4)+
     ggplot2::geom_text(data=NULL, mapping=aes(y='Colour', x=0.05,label =colorComponent,hjust=0), size=4)+
     ggplot2::geom_text(data=NULL, mapping=aes(y='Head', x=0.05,label ="Legend info",hjust=0), size=4)
   pcaPlot <- putPlot(pcaPlot, plot, n.comp+1, 1)
   
   
  plot <- ggplot2::ggplot()
  plot <- plot +
    ggplot2::scale_y_discrete(name="",  limits = c( as.character(rev(breakInfo$Both)),"Legend") ) +
    ggplot2::scale_x_continuous(limits=c(0,1), breaks=NULL, name="") + 
    ggplot2::geom_point(data=breakInfo, mapping=aes(y=Both, x=0.10,colour=colourValues,  shape = shapeLevels, linetype = shapeLevels), size=3)+
    ggplot2::geom_text(data=breakInfo, mapping=aes(y=Both, x=0.25,label =Both,hjust=0), size=4)+
    ggplot2::geom_text(data=NULL, mapping=aes(y='Legend', x=0.5,label ="Legend"), size=4)
  pcaPlot <- putPlot(pcaPlot, plot, n.comp+1, 2)
  
  plot <- ggplot2::ggplot()
  plot <- plot +
    ggplot2::scale_y_discrete(name="",  limits = c( as.character(rev(e.var$components)),"PCs") ) +
    ggplot2::scale_x_continuous(limits=c(0,1), breaks=NULL, name="") + 
    ggplot2::geom_text(data=e.var, mapping=aes(y=components, x=0.05,label =info,hjust=0), size=4)+
    ggplot2::geom_text(data=NULL, mapping=aes(y='PCs', x=0.5,label ="PCs"), size=4)
  pcaPlot <- putPlot(pcaPlot, plot, 1, 2)


  personal_plot <- ggally_text(
    "ggpairs allows you\nto put in your\nown plot.\nLike that one.\n <---"
  )
  pcaPlot <- putPlot(pcaPlot, personal_plot, , 1)
  
  
  
  plot <- plot +
    ggplot2::scale_y_discrete(name="") +
    ggplot2::scale_x_continuous(limits=c(0,1), breaks=NULL, name="") + 
    ggplot2::geom_point(data=breakInfo, mapping=aes(y=Both, x=0.10,colour=colourValues,  shape = shapeLevels, linetype = shapeLevels), size=3)+
    ggplot2::geom_text(data=breakInfo, mapping=aes(y=Both, x=0.25,label =Both,hjust=0), size=4)

  pcaPlot <- putPlot(pcaPlot, plot, n.comp+1, 1)
  
      
  for(i in 1:(n.comp-1)){
    for(j in  (i+1):n.comp){
      plot <- ggplot2::ggplot()
      plot <- plot +
        ggplot2::scale_y_discrete(name="") +
        ggplot2::scale_x_continuous(limits=c(0,1), breaks=NULL, name="") + 
        ggplot2::geom_point(data=breakInfo, mapping=aes(y=Both, x=0.10,colour=colourValues,  shape = shapeLevels, linetype = shapeLevels), size=3)+
        ggplot2::geom_text(data=breakInfo, mapping=aes(y=Both, x=0.25,label =Both,hjust=0), size=4)
      
        pcaPlot <- putPlot(pcaPlot, plot, n.comp+1, j)
    }
  }

  
  
  pcaPlot
  
clNames = colnames(ggInfo)

m = ggplot(ggInfo, aes(x = PC1))
m <- ggplot(ggInfo2, aes(x=PC1, colour=Timepoint, linetype=Phenotype))
m + geom_density(fill=NA,adjust=1/3)



# Create PCAplot requires that the colnames in the expressionTable are present in the metaDataTable[sampleColName, ] 

createTwoPCAplot <- function(expressionTable, metaDataTable, xComp = 1,YComp = 2 ,
                                 sampleColName = "Name", colorComponent = "Phenotype" , pchComponent = "Date"
                                 ,PCAcolors =  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")
                                 ,pchs = 16:40,
                                 scale = TRUE, center = TRUE){
  ######################
  #PCA
  ######################
  
  
  metaInfo = metaDatatable[metaDatatable[[sampleColName]] %in% colnames(expressionTable),  ]
  rownames(metaInfo) <- metaInfo[[sampleColName]]
  
  #components
  plot.comp = 1:n.comp
  
  
  #pca
  pca = prcomp(t(expressionTable), scale=scale, center=center)
  pca.basis = pca[['x']]
  e.var = pca[['sdev']]^2 / sum(pca[['sdev']]^2)
  names(e.var) = colnames(pca.basis)
  
  #color maps
  ColorVariables = as.factor(metaInfo[,colorComponent])
  SampleColors = PCAcolors[1:length(levels(ColorVariables))]
  names(SampleColors) =  levels(ColorVariables)
  metaInfo['Color'] = as.character(revalue(ColorVariables,SampleColors))
  
  #pch variables maps
  pchVariables = as.factor(metaInfo[,pchComponent])
  SamplePCHs = pchs[1:length(levels(pchVariables))]
  names(SamplePCHs) =  levels(pchVariables)
  metaInfo['pchs'] = as.numeric(revalue(pchVariables,SamplePCHs))
  
  
  
  samples = rownames(pca.basis)
  PCs = colnames(pca.basis)
  
  
  stage.color = as.character(metaInfo[samples, 'Color'])
  stage.pchs = metaInfo[samples, 'pchs']
  metaInfo$Both = paste( metaInfo[[colorComponent]],metaInfo[[pchComponent]], sep = " ") 
  stage2color.map = unique(metaInfo[samples, c('Both','Color','pchs')])
  
  
  
  #plot
  ##    pdf(pca.pairs.preqc.pdf)
  
  p = pairs(pca.basis[, plot.comp], upper.panel = point.text.panel, lower.panel = point.panel, diag.panel = diag.panel, text.panel = text.panel)
  
  ##    dev.off()
  
  return (p)
}

# ggInfo =  p  + geom_point() + 
#   scale_colour_manual(name = paste(pchComponent, colorComponent, sep = "_"),
#                       breaks = breakInfo$Both,
#                       values = breakInfo$colourValues) +   
#   scale_shape_manual(name = paste(pchComponent, colorComponent, sep = "_"),
#                      breaks = breakInfo$Both,
#                      values = breakInfo$shapeValues) 
# )



# point.text.panel <- function(x, y){
#   xmin = min(y)
#   xmax = max(y)
#   xbuffer  = as.integer((xmax-xmin)/10 +2 )
#   ymin = min(x)
#   ymax = max(x)  
#   ybuffer  = as.integer((ymax-ymin)/10 +2 )
#   par("usr" = c(xmin-xbuffer,xmax+xbuffer, ymin-ybuffer, ymax+ybuffer))
#   
#   text(y, x, labels = samples, col = stage.color)  
#   
# }
# 
# point.panel <- function(x, y){
#   xmin = min(x)
#   xmax = max(x)
#   xbuffer  = as.integer((xmax-xmin)/10 +2 )
#   ymin = min(y)
#   ymax = max(y)  
#   ybuffer  = as.integer((ymax-ymin)/10 +2 )
#   par("usr" = c(xmin-xbuffer,xmax+xbuffer, ymin-ybuffer, ymax+ybuffer))
#   
#   points(x, y, col = stage.color, pch = stage.pchs)
# }
# 
# diag.panel <- function(x, y, labels, ...){
#   pu <- par("usr")
#   d <- density(x,bw = 0.5, ...)
#   par("usr" = c(pu[1:2], 0, max(d$y)*1.5))
#   
#   lines(d)
# }
# 
# text.panel <- function(x, y, labels, ...){
#   
#   if(labels == 'PC1'){
#     text(0.5, 0.3, sprintf('%s: %.1f%s', labels, e.var[grep(paste('^', labels, '$', sep = ''), colnames(pca.basis))] * 100, '%'));
#     legend('topleft', legend = stage2color.map[, 'Both'], col = stage2color.map[, 'Color'], pch = stage2color.map[, 'pchs'])
#   }
#   else{
#     text(0.5, 0.8, sprintf('%s: %.1f%s', labels, e.var[grep(paste('^', labels, '$', sep = ''), colnames(pca.basis))] * 100, '%'));
#   }
# }
# 
# createPCAplot(expressionTable = expressionTable, metaDataTable = metaDataPremiRNA, n.comp = 5 , sampleColName = "Name", colorComponent = "Date", pchComponent = "Phenotype")
# 
# 
# # Create PCAplot requires that the colnames in the expressionTable are present in the metaDataTable[sampleColName, ] 
# createMultipePCAplot <- function(expressionTable, metaDataTable, n.comp = 4,
#                                  sampleColName = "Name", colorComponent = "Phenotype" , pchComponent = "Date"
#                                  ,PCAcolors =  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000")
#                                  ,pchs = 16:40,
#                                  scale = TRUE, center = TRUE){
#   ######################
#   #PCA
#   ######################
#   
#   
#   metaInfo = metaDatatable[metaDatatable[[sampleColName]] %in% colnames(expressionTable),  ]
#   rownames(metaInfo) <- metaInfo[[sampleColName]]
#   
#   #components
#   plot.comp = 1:n.comp
#   
#   pca = prcomp(t(expressionTable),scale=FALSE, center=TRUE )  
#   
#   #pca
#   pca = prcomp(t(expressionTable), scale=scale, center=center)
#   pca.basis = pca[['x']]
#   e.var = pca[['sdev']]^2 / sum(pca[['sdev']]^2)
#   names(e.var) = colnames(pca.basis)
#   
#   #color maps
#   ColorVariables = as.factor(metaInfo[,colorComponent])
#   SampleColors = PCAcolors[1:length(levels(ColorVariables))]
#   names(SampleColors) =  levels(ColorVariables)
#   metaInfo['Color'] = as.character(revalue(ColorVariables,SampleColors))
#   
#   #pch variables maps
#   pchVariables = as.factor(metaInfo[,pchComponent])
#   SamplePCHs = pchs[1:length(levels(pchVariables))]
#   names(SamplePCHs) =  levels(pchVariables)
#   metaInfo['pchs'] = as.numeric(revalue(pchVariables,SamplePCHs))
#   
#   
#   
#   samples = rownames(pca.basis)
#   PCs = colnames(pca.basis)
#   
#   
#   stage.color = as.character(metaInfo[samples, 'Color'])
#   stage.pchs = metaInfo[samples, 'pchs']
#   metaInfo$Both = paste( metaInfo[[colorComponent]],metaInfo[[pchComponent]], sep = " ") 
#   stage2color.map = unique(metaInfo[samples, c('Both','Color','pchs')])
#   
#   
#   
#   #plot
#   ##    pdf(pca.pairs.preqc.pdf)
#   
#   p = pairs(pca.basis[, plot.comp], upper.panel = point.text.panel, lower.panel = point.panel, diag.panel = diag.panel, text.panel = text.panel)
#   
#   ##    dev.off()
#   
#   return (p)
# }
# 
# 

info2 = ggInfo[,"PC1"]
metaInfo2 = ggInfo[,"Date"]

compareValue(info2, metaInfo2){
 Levels = levels(metaInfo2)
 lowestPvalue = 1;
 for(i = 1:(length(Levels)-1)){
   for(j =  (i+1):(length(Levels))){
     A = info2[which(metaInfo2 == Levels[i])]
     B = info2[which(metaInfo2 == Levels[j])]
     KS.out = ks.test(A,B)
     if(KS.out$p.value < lowestPvalue){
       lowestPvalue
     }
     KS.out$data.name
   }
 } 
 
  
}


