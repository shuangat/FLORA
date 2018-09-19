library(gProfileR)
library(data.table)

calculateCorrelation <- function(file) {
  if (file.exists(file) == TRUE) {
    data <- read.table(file, sep=",", header = T, row.names = 1)
    
    n.coding <- nrow(subset(data,info == "coding")) #i 
    n.lnc <- nrow(subset(data, info == "lnc"))   #j
    
    data <- data[,-1]
    data<-data.frame(cbind(data[,1],data.frame( log2(as.matrix(data[,-1]) +0.1) )))
    names(data)[1] <- "gene"
    
    coding <- data[1:n.coding,]
    lnc <- data[(n.coding+1):(n.coding+n.lnc),]
    
    pearson <- cor(t(data[,-1]))
    pearson <- pearson[1:n.coding,(n.coding+1):(n.coding+n.lnc)]
    pearson <- data.frame(pearson)
    
    correlation <- cbind( coding[,1], data.frame(abs(pearson)) )
    names(correlation)[1] <- "gene"
    
    coding.info <- cbind(row.names(coding),data.frame(coding[,1]) ) 
    colnames(coding.info) <- c("id","name")
    
    lnc.info <- cbind(row.names(lnc),data.frame(lnc[,1]) )
    colnames(lnc.info) <- c("id","name")
    
    return(list(correlation=correlation, coding.info=coding.info, lnc.info=lnc.info))
    
  } else {
    return("can not open the file!")
  }
}


functionalPrediction <- function(lnc.name, file_correlation, file_lnc.info) {
  if ( (file.exists(file_correlation) == TRUE) &(file.exists(file_lnc.info) == TRUE)) {
    correlation <- fread(file_correlation,sep="\t")
    lnc.info <- read.table(file_lnc.info, sep = "\t", head=T, row.names=1)
    
    j <- which(lnc.info$name == lnc.name)
    if (length(j) > 0) {
      k <- j+2
      od <- data.frame(order(correlation[,..k],decreasing = T))
      correlation.top <- correlation[od[1:30,],c(2,..k)]
      gene <- data.frame(correlation.top[,1])
      
      go <- gprofiler(gene, organism = "hsapiens",ordered_query = T)
      go <- go[order(go[,3]),]
      
      BP <- subset(go, domain == "BP")
      CC <- subset(go, domain == "CC")
      MF <- subset(go, domain == "MF")
      
      return(list(GO=go, BP=BP, CC=CC, MF=MF))
    } else {
      return("can not find the lncRNA")
    }
    
  } else {
    return("can not open files!")
  }
}

