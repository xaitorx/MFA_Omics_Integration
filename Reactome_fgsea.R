## Download Reactome annotation from Reactome website
## All levels of the pathway hierarchy
# ChEBI to All pathways
# ENSEMBL to All pathways
# miRBase to All pathways
# Paste them together rbind()
library(plyr)
library(fgsea)
library(ggplot2)
library(RCurl)

### data frame with ENSEMBL/CHEBI/miRNA identifiers annotated to reactome pathways
### split into 1 list of IDs/pathway
Reactome_all <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/MFA_Omics_Integration/v1/data/Reactome_all.csv"), stringsAsFactors=FALSE)
Reactome_all <- Reactome_all[order(Reactome_all$pathway),]

#count numer of elements/pathway
counts <- as.data.frame(table(Reactome_all$pathway))
magic_number <- length(unique(counts$Var1))

breaks <- rep(1:magic_number, counts$Freq)

# break IDs into pathways, rename each list accordingly
lista1 <- split(Reactome_all$ID, breaks)
names(lista1) <- unique(counts$Var1)

# save for later
anotacion <- as.data.frame(unique(Reactome_all[,2:3]))


### Perform GSEA with these lists on our preranked list of IDs
var_dim2_annot <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/MFA_Omics_Integration/v1/data/var_dim2_annot.csv"))

ranking <- var_dim2_annot[,c(7,4)]
ranking <- ranking[!is.na(ranking$REACTOME),]

# consolidate duplicated IDs, order ranking
ranking <- ddply(ranking,1,numcolwise(mean))
ranking <- ranking[order(ranking$contrib, decreasing = TRUE),]

ranking_ooo <- ranking$contrib
names(ranking_ooo) <- ranking$REACTOME

###Perform GSEA with list of reactome pathways
fgseaRes <- fgsea(pathways = lista1, 
                  stats = ranking_ooo,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)


fgseaRes_annot <- merge(anotacion, fgseaRes, by.x = 1, by.y = 1)
write.csv(fgseaRes_annot[,-9], "fgseaRes_annot.csv")

### plot
plotEnrichment(lista1[[(fgseaRes[130,])$pathway]],ranking_ooo) +
  ggtitle("Triglyceride catabolism", size) 
  # change for number of row


