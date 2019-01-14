source('functionalPredictio.R')

# Load information of all coding genes and lncRNAs
lnc.info <- read.table("lnc.info.txt", sep = "\t", head=T, row.names=1)
coding.info <- read.table("coding.info.txt", sep = "\t", head=T, row.names=1)

# Load gene regulatory network constructed by ARACNe-AP
network <- read.table("network.txt", sep = "\t", head=T)

# Input the name of your lncRNA of interest
lnc.name <- "LINC01614"

# Extract genes connected with the lncRNA and perform functional prediction
lnc.coding <- getnetwork(lnc.info, coding.info, network, lnc.name)
results <- makePrediction(lnc.name, lnc.coding)
gene.use <- results$gene.use
GO <- subset(results$GO, domain %in% c("BP", "CC", "MF"))
write.table(GO, paste(lnc.name,"_GO.txt",sep=""), sep="\t")

# Plotting significant GO terms associated with the networks of the lncRNA
library(ggplot2)
pdf( file=paste(lnc.name,"_GO.pdf",sep = "")) 
GO$term.name = factor(GO$term.name, levels = GO$term.name)
ggplot(aes(x = term.name, y = -log10(GO$p.value), fill = domain), data = GO) + 
  geom_bar(stat="identity",position="dodge") +
  labs(x = "", y = "-log10_Pvalue", title= lnc.name )  + coord_flip() +
  theme_classic() + scale_fill_manual(values = c('#8470FF','#87CEFA','#FFC125')) +
  theme(legend.title=element_blank(), axis.title.x = element_text(size=9) )
dev.off()
