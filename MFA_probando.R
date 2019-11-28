library(RCurl)
library(FactoMineR)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(factoextra)


log2_mfa_data <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/MFA_Omics_Integration/v1/data/log2_mfa_data.csv"), row.names=1)
log2_mfa_data <- as.data.frame(t(log2_mfa_data))

res <- MFA(log2_mfa_data, group=c(15445, 212, 798), type=c(rep("c",3)),
           ncp=5, name.group=c("mRNA","miRNA","metabolome"))

# A summary of the main results of the MFA function
summary(res)

### Barplot of variance explained by each dimension
eig <- as.data.frame(res$eig)
variances <- barplot(eig$`percentage of variance`, names.arg=1:nrow(eig),main = "Variance explained by dimension", xlab = "Dimensions",ylab = "Percentage of variance", col ="steelblue")
lines(x = variances, eig$`percentage of variance`, type="b", pch=19, col = "red")

### each "layer" contribution to the different components
contrib <- res$group$contrib
my_palette <- colorRampPalette(c("white", "steelblue", "navy", "black", "black"))(n = 500)
myBreaks <- c(seq(0, 100, length.out=ceiling(500) + 1))
pheatmap(contrib, cluster_cols=FALSE,show_rownames = TRUE,breaks = myBreaks, color = my_palette, legend = TRUE, border_color = 1, fontsize_col = 16, fontsize_row = 16)

### Graph of the groups
group_coord <- as.data.frame(res$group$coord)
ggplot(group_coord , aes(Dim.1, Dim.2, col= c("mRNA", "miRNA", "metabolome"))) + geom_point(aes(size=5)) + geom_text(aes(size=50,label = row.names(group_coord), vjust="inward",hjust="inward"), nudge_y = 0.01, parse = TRUE) + 
  labs (y = "Dim2 (26%)", x = "Dim1 (34%)") +
  theme_light()
fviz_mfa_axes(res)
#relationship between dimensions of the MFA and the dimensions from individual group PCAs

### Graph of the individuals
fviz_mfa_ind(res, col.ind = "cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
fviz_mfa_ind(res, partial = "all", repel = TRUE)
# cos2 = quality of representation
# R4-5 only metabolomics data. no mRNA/miRNA

# better looking individuals graph with ggplot2
ind_coord <- as.data.frame(res$ind$coord)
ggplot(ind_coord , aes(Dim.1, Dim.2, col= c(rep("DMSO", 5), rep("MKC", 5)))) + geom_point(aes(size=2)) + geom_text_repel(aes(label = row.names(ind_coord)), nudge_y = 0.25, fontface = "bold") + 
  labs (y = "Dim2 (26%)", x = "Dim1 (34%)") +
  theme_light()

# confidence ellipse (90%) + center
ggplot(ind_coord , aes(Dim.1, Dim.2, col= c(rep("DMSO", 5), rep("MKC", 5)))) + geom_point(aes(size=2)) + 
  stat_ellipse(geom = "path", alpha = 0.5, level = 0.9, linetype = 2, size=3) +
  stat_ellipse(geom = "polygon", type = "euclid", aes(fill = c(rep("DMSO", 5), rep("MKC", 5))), alpha = 0.5, level = 0.9) + 
  geom_text_repel(aes( label = row.names(ind_coord)), nudge_y = 0.25) +
  labs (y = "Dim2 (26%)", x = "Dim1 (34%)") +
  theme_light()

### Inspect variables relation with MFA dimensions
# Focus in dimension 2 (this dimension kind of separates individuals by biological group)
coord_var <- as.data.frame(res$quanti.var$coord)
contrib_var <- as.data.frame(res$quanti.var$contrib)
cor_var <- as.data.frame(res$quanti.var$cor)
var_Dim2 <- as.data.frame(cbind(coord_var$Dim.2, contrib_var$Dim.2 , cor_var$Dim.2))
colnames(var_Dim2) <- c("Dim2", "contrib", "correlation")
row.names(var_Dim2) <- row.names(coord_var)

# annotate (there is duplicate names for some variables)
#load annotation file
todo_annot <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/MFA_Omics_Integration/v1/data/todo_annot.csv"), row.names=1, header=T, na.strings=c("","NA"))

var_Dim2_annot <- cbind(todo_annot[match(row.names(var_Dim2), todo_annot$X),], var_Dim2)
var_Dim2_annot$type <- c(rep("mRNA", 15445), rep("miRNA", 212), rep("metabolome", 798))
write.csv(var_Dim2_annot, "variables_dimension2.csv", row.names = FALSE)

# contrib vs coordinates vs correlations
ppp <- ggplot(var_Dim2_annot, aes(log2(abs(var_Dim2_annot$Dim2)), log2(var_Dim2_annot$contrib))) 
ppp + geom_point(aes(colour = var_Dim2_annot$type, size = 0.1, alpha = 0.5)) +
  labs (y = "log2(Dim2 contribution)", x = "log2(abs(Dim2 coordinates(scores)))")

ppp <- ggplot(var_Dim2_annot, aes(log2(var_Dim2_annot$contrib), var_Dim2_annot$correlation)) 
ppp + geom_point(aes(colour = var_Dim2_annot$type, size = 0.1, alpha = 0.5)) +
  labs (y = "Dim2 correlation", x = "log2(var_Dim2$contrib)")

ppp <- ggplot(var_Dim2_annot, aes(var_Dim2_annot$Dim2, var_Dim2_annot$correlation))
ppp + geom_point(aes(colour = var_Dim2_annot$type, size = 0.1, alpha = 0.5)) +
  labs (y = "Dim2 correlation", x = "Dim2 coordinates(scores)") +
  coord_cartesian(xlim = c(-5, 5))


### inspect top weight variables from each group

var_Dim2_annot$Y <- as.character(var_Dim2_annot$Y)
var_Dim2_annot$Y[is.na(var_Dim2_annot$Y)] <- as.character(var_Dim2_annot$X[is.na(var_Dim2_annot$Y)])
# keep ENSEMBL IDs for the events without Gene Symbols

var_Dim2_annot <- var_Dim2_annot[order(var_Dim2_annot$contrib, decreasing = TRUE),]
var_Dim2_annot$Y <- factor(var_Dim2_annot$Y ,levels = unique(var_Dim2_annot$Y))
# order the rows of data.frame according to contribution (for example)

### Inspect top key variables 
grafico <- ggplot(var_Dim2_annot[1:20,], aes(y=contrib, x=Y, fill=var_Dim2_annot$type[1:20]))
grafico + geom_segment(aes(colour=var_Dim2_annot$type[1:20], x=Y, xend=Y, y=-5,size=1, yend=contrib-1)) +   
  geom_point( aes(colour=var_Dim2_annot$type[1:20]), size=6, alpha=0.8, shape=21, stroke=2) +   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top contributors to Dimension2",y = "Contribution (%)", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(ylim = c(0, 30))

# can do it with Factoextra aswell
fviz_contrib(res, choice = "quanti.var", axes = 2, top = 20)

### Inspect top variables from each group
mRNA_Dim2 <- subset(var_Dim2_annot, var_Dim2_annot$type== "mRNA")
miRNA_Dim2 <- subset(var_Dim2_annot, var_Dim2_annot$type== "miRNA")
metabolome_Dim2 <- subset(var_Dim2_annot, var_Dim2_annot$type== "metabolome")

grafico <- ggplot(mRNA_Dim2[1:20,], aes(y=contrib, x=Y))
grafico + geom_segment(colour="steelblue", aes(x=Y, xend=Y, y=-0.5,size=1, yend=contrib-0.01)) +   
  geom_point( size=6, color="steelblue", fill=alpha("steelblue"), alpha=0.8, shape=21, stroke=2)+   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top mRNA contributors",y = "Contribution (%)", x = "") + 
  coord_cartesian(ylim = c(0, 0.4)) + 
  theme(plot.title = element_text(hjust = 0.5))


grafico <- ggplot(miRNA_Dim2[1:20,], aes(y=contrib, x=Y))
grafico + geom_segment(colour="limegreen", aes(x=Y, xend=Y, y=-5,size=1, yend=contrib-1)) +   
  geom_point( size=6, color="limegreen", fill=alpha("limegreen"), alpha=0.8, shape=21, stroke=2)+   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top miRNA contributors",y = "Contribution (%)", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 30)) 


levels(metabolome_Dim2$Y) <- c(levels(metabolome_Dim2$Y), "D-Fructose-6-phosphate")
metabolome_Dim2[18,2] <- as.factor("D-Fructose-6-phosphate")
metabolome_Dim2$Y <- factor(metabolome_Dim2$Y ,levels = unique(metabolome_Dim2$Y))
grafico <- ggplot(metabolome_Dim2[1:20,], aes(y=contrib, x=Y))
grafico + geom_segment(colour="firebrick1", aes(x=Y, xend=Y, y=-5,size=1, yend=contrib-0.01)) +   
  geom_point( size=6, color="firebrick1", fill=alpha("firebrick1"), alpha=0.8, shape=21, stroke=2)+   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top metabolome contributors",y = "Contribution (%)", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 0.4))


  + scale_y_continuous(expand = c(0, 0)) # if you dont want to "expand" axis
