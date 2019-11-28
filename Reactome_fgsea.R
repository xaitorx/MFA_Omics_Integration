## Download Reactome annotation from Reactome website
## All levels of the pathway hierarchy
# ChEBI to All pathways
# ENSEMBL to All pathways
# miRBase to All pathways
# Paste them together
library(plyr)
library(fgsea)

### data frame with ENSEMBL/CHEBI/miRNA identifiers annotated to reactome pathways
### split into 1 list of IDs/pathway
Reactome_all <- read.csv("C:/Users/Aitor/Desktop/Grecia/MFA/Reactome/Reactome_all.csv", stringsAsFactors=FALSE)
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
var_dim2_annot <- read.csv("C:/Users/Aitor/Desktop/Grecia/MFA/var_dim2_annot.csv")

ranking <- var_dim2_annot[,c(7,4)]
ranking2 <- ranking[!is.na(ranking$REACTOME),]

# consolidate duplicated IDs, order ranking
ranking3 <- ddply(ranking2,1,numcolwise(mean))
ranking3 <- ranking3[order(ranking3$contrib, decreasing = TRUE),]

rankings <- ranking3$contrib
names(rankings) <- ranking3$REACTOME

###Perform GSEA with list of reactome pathways
fgseaRes <- fgsea(pathways = lista1, 
                  stats = rankings,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)


fgseaRes_annot <- merge(anotacion, fgseaRes, by.x = 1, by.y = 1)

