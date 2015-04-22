
plotLengthDistribution <- function(AllLength, startLength = 1, stopLength = dim(AllLength)[1], lty = 1, col = 1 ){
  tempLength = AllLength[startLength:stopLength,]
  matplot(rownames(tempLength),y = tempLength,type = "l",xlab = "nt Length", ylab = "raw counts",lty =  lty , col = col)
  legend("topright", c(colnames(tempLength)),cex = 0.6, lty = lty , col = col)
}


