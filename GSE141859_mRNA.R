# Version info: R 4.3.2, Biobase 2.62.0, GEOquery 2.70.0, limma 3.58.1
#
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(affy)
library(umap)
library(frma)
library(oligo)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)



# ##### GSE141859 ######
# # load series and platform data from GEO
# # # Path for raw .tar files
# path_gse141859 <- "/Users/selcukkorkmaz/Documents/Studies/AMI_miRNA/GSE141859_RAW"
# 
# # 
# # # untar files
# untar(paste0(path_gse141859, ".tar"), exdir = "/Users/selcukkorkmaz/Documents/Studies/AMI_miRNA/GSE141859/")
# 
# # Reading CEL files
# celFiles <- list.files("/Users/selcukkorkmaz/Documents/Studies/AMI_miRNA/GSE141859/")
# full_cel_paths <- paste0("/Users/selcukkorkmaz/Documents/Studies/AMI_miRNA/GSE141859/", celFiles)
# celFiles_gse59867 <- ReadAffy (filenames = full_cel_paths)


# Perform RMA (Robust Multi-array Average) normalization
# This step includes background correction, normalization, and log2 transformation

# gset <- getGEO("GSE141859")
# if (length(gset) > 1) idx <- grep("GPL22760", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]


# gset <- frma(celFiles_gse59867)
# 
# # Load group information
# groups_df = read_excel("data/GSE141859_groups.xlsx")
# ids = as.data.frame(colnames(gset))
# colnames(ids) = "Accession"
# ids$Accession = gsub("\\..*","",ids$Accession)
# 
# # Merge group information with sample IDs
# groups_df2 <- left_join(ids, groups_df, by = "Accession")

counts = read.csv("GSE141859_control_vs_comparison.deseq.counts.csv")
colnames(counts)[1] = "ID"

# Load necessary libraries
library(limma)
library(edgeR)

# Step 2: Load your data
# Step 3: Create a DGEList object
y <- DGEList(counts = counts)
y$samples$group = c(1,1,1,1,2,2,2,2)

# Step 4: Filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,]

# Step 5: Normalization
y <- calcNormFactors(y)

# Step 6: Design matrix
group <- factor(c(rep("Control", 4), rep("MI", 4)))
design <- model.matrix(~0 + group)

# Step 7: Differential expression analysis using voom and limma
y <- voom(y, design, plot = TRUE)
fit <- lmFit(y, design)
contrast.matrix <- makeContrasts(MIvsControl = groupMI - groupControl, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Step 8: Results
# results <- topTable(fit2, coef = "MIvsControl", n = Inf)
# write.csv(results, file = "/path/to/your/results.csv")

# View the top results
# head(results)


groups = as.factor(c("Control","Control","Control","Control", "MI","MI","MI","MI"))
group_names <- make.names(c("Control","MI"))


# # Log2 transformation check and application
# ex <- gset[,-1]
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0)
# 
# if (LogC) { #ex[ex <= 0] <- NaN
# gset_exp <- log2(ex) }else{
#   gset_exp <- ex
# }
# 
# 
# # Preprocessing of transformed data
# gset_exp[gset_exp <= 0] <- NaN
# rownames(gset_exp) = gset$X
# gset_exp <- gset_exp[complete.cases(gset_exp),]
# gset_t = gset_exp  %>% t() %>% as.data.frame()
# gset_t$group <- groups
# 
# 
# # Box plots before normalizing
# ex <- as.matrix(gset_t[,-ncol(gset_t)] %>% t())
# ex_df <- as.data.frame(ex)
# ex_df$SampleID <- rownames(ex_df)
# long_df <- melt(ex_df, id.vars = "SampleID", variable.name = "miRNA", value.name = "Expression")
# 
# long_df <- long_df %>%
#   mutate(Group = rep(groups, each = nrow(long_df) / length(groups)))
# 
# long_df$miRNA_Group <- interaction(long_df$miRNA, long_df$Group, sep = " - ")
# head(long_df)
# 
# before=ggplot(long_df, aes(x = miRNA_Group, y = Expression, fill = Group)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Rotate x labels for readability
#   labs(x = "", y = "") +
#   scale_fill_manual(values = c("Control" = "#1B9E77", "MI" = "#7570B3"))+
#   theme(legend.text=element_text(size=10), legend.title=element_blank())+
#   scale_x_discrete(breaks=long_df$miRNA_Group,
#                    labels=long_df$miRNA)+ theme(legend.text=element_text(size=28))
# 
# # Normalize expression data
# gset_exp_norm <- normalizeBetweenArrays(gset_exp, method = "quantile") # normalize data
# 
# 
# # Prepare for linear modeling
# gset_exp2 = gset_exp_norm  %>% t() %>% as.data.frame()
# colnames(gset_exp2) = gset$ID_REF
# gset_exp2$group <- groups
# write.table(gset_exp2, "data/GSE141859_expression.txt", quote = F, sep = "\t")
# design <- model.matrix(~group + 0, gset_exp2)
# colnames(design) <- levels(groups)
# 
# # Vooma transformation and linear modeling
# transposed_gset_exp2 <- as.matrix(gset_exp2[,-ncol(gset_exp2)]) %>% t()
# v <- voomaByGroup(transposed_gset_exp2, group=group_names, design, plot=T, cex=0.1, pch=".", col=1:nlevels(groups))
# v$genes <- rownames(gset_exp)
# 
# # Fit linear model to the voom-transformed data
# 
# fit  <- lmFit(v)
# 
# # Set up contrasts of interest and recalculate model coefficients
# cts <- c(paste(group_names[1],"-",group_names[2],sep=""))
# cont.matrix <- makeContrasts(contrasts=cts, levels=design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# 
# # Compute statistics and table of top significant genes
# fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","AveExpr"))
# gpl8179 = getGEO("GPL6101")
# annotation = gpl8179@dataTable@table
# topGenes = annotation[annotation$ID %in% tT$ID, c("ID", "Symbol")]
# topGenes = left_join(tT, topGenes, by="ID")
# topGenes = tT %>% dplyr::select(c(1,8, 2:7))
write.table(topGenes, file="data/TopGenes_GSE141859.txt", row.names=F, sep="\t", quote = F)


# Summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0.5)

# Filter and save differentially expressed genes
degs = dT@.Data[dT@.Data[,1]!=0,]
length(degs)
degs_data = tT[abs(tT$logFC)>=0.5 & tT$adj.P.Val<0.05,]
# degs_data = degs_data[degs_data$ID %in% names(degs),]
degs_data <- subset(degs_data, select=c("ID","adj.P.Val","P.Value","t","B","logFC"))
probe_ids <- degs_data$ID
miRNAs = annotation[annotation$ID %in% degs_data$ID, c("ID", "Symbol")]
colnames(miRNAs)[1] = "ID"
degs_data = left_join(degs_data, miRNAs, by="ID")
degs_data = degs_data %>% dplyr::select(c(1,7, 2:6))
write.table(degs_data, file="data/DEGs_GSE141859.txt", row.names=F, sep="\t", quote = F)

# Creating plots using layout: voom:mean-variance trend plot, moderated t statistic plot, 
# p-adj value distribuution plot, expression value distribuution plot

# General expression data analysis
ex <- as.matrix(gset_exp2[,-ncol(gset_exp2)] %>% t())
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
t.good <- which(!is.na(fit2$F)) 
layout <- layout(matrix(c(1,2,3,4), 2, 2, byrow = F)) 
layout <- voomaByGroup(transposed_gset_exp2, group=group_names, design, plot=T, cex=0.1, pch=".", col=1:nlevels(groups))
layout <- hist(tT$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
               ylab = "Number of genes", main = "P-adj value distribution")
layout <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

layout <- plotDensities(ex, group=groups, main="GSE141859: Expression value distribution", legend ="topright")

# Box plot after notmalization
ex_df <- as.data.frame(ex)
ex_df$SampleID <- rownames(ex_df)

# Melt the data frame to long format
long_df <- melt(ex_df, id.vars = "SampleID", variable.name = "miRNA", value.name = "Expression")

# Add group information to the long data frame
long_df <- long_df %>%
  mutate(Group = rep(groups, each = nrow(long_df) / length(groups)))

# Create a new variable for interaction between miRNA and Group to plot side by side
long_df$miRNA_Group <- interaction(long_df$miRNA, long_df$Group, sep = " - ")
head(long_df)

# Now, use ggplot2 to create the box plot with Control and HCM side by side
after =ggplot(long_df, aes(x = miRNA_Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Rotate x labels for readability
  labs(x = "", y = "") +
  scale_fill_manual(values = c("Control" = "#1B9E77", "HCM" = "#7570B3"))+
  theme(legend.text=element_text(size=10), legend.title=element_blank())+
  scale_x_discrete(breaks=long_df$miRNA_Group,
                   labels=long_df$miRNA)+ theme(legend.text=element_text(size=28))


p2 = before / after

p = p+ plot_annotation(tag_levels = 'A')

ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/GSE141859_RLE.png', p, width = 16, height = 8)


# Volcano plot

adjPVal_cutoff = 0.05
logFC_cutoff = 2.5

down <- rownames(tT)[tT$adj.P.Val < adjPVal_cutoff & tT$logFC < -logFC_cutoff]
up <- rownames(tT)[tT$adj.P.Val < adjPVal_cutoff & tT$logFC > logFC_cutoff]

## Print volcano plots with up genes in purple and down genes in green
keyvals.colour1 <- ifelse(
  rownames(tT) %in% up, 'red',
  ifelse(rownames(tT) %in% down, 'blue',
         'black'))

keyvals.colour1[is.na(keyvals.colour1)] <- 'black'
names(keyvals.colour1)[keyvals.colour1 == 'black'] <- 'Not significant'
names(keyvals.colour1)[keyvals.colour1 == 'red'] <- 'Up-regulated'
names(keyvals.colour1)[keyvals.colour1 == 'blue'] <- 'Down-regulated'


p2 <-  EnhancedVolcano(tT,
                       lab = NA,
                       x = 'logFC',
                       y = 'adj.P.Val',
                       xlim = c(-4,4),
                       title = NULL,  
                       subtitle = NULL, 
                       caption = NULL,
                       pCutoff = adjPVal_cutoff,
                       FCcutoff = logFC_cutoff,
                       pointSize = c(ifelse(tT$adj.P.Val < adjPVal_cutoff & abs(tT$logFC) > logFC_cutoff, 4, 2)),
                       colCustom = keyvals.colour1,
                       legendLabSize = 12,
                       legendIconSize = 8,
                       axisLabSize = 14,
                       legendPosition = 'bottom')




# Mean-Difference Plot of Expression Data
df <- data.frame(
  logFC = tT$logFC,
  AveExpr = tT$AveExpr,
  Significant = names(keyvals.colour1)   # Assuming 'dT' indicates significant genes as non-zero
)

# Create the Mean-Difference Plot using ggplot2
p2 <- ggplot(df, aes(x = AveExpr, y = logFC, color = Significant)) +
  geom_point(alpha = 0.6) +  # Adjust point transparency with 'alpha'
  scale_color_manual(values = c("Not significant" = "grey7", "Up-regulated" = "red", "Down-regulated" = "blue")) +  # Change colors as needed
  labs(x = "Average expression", y = bquote(~Log[2] ~ "fold change")) +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if not needed


# UAMP plot
df <- data.frame(UMAP1 = ump$layout[,1], UMAP2 = ump$layout[,2], Group = groups)

# Creating the UMAP plot with ggplot2 and ggrepel for label repulsion
p3 <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(aes(shape = Group), size = 3) +  # Adjust size as needed
  scale_shape_manual(values=c(19,19)) +  # Ensure unique shapes for groups if desired
  # geom_text_repel(aes(label = rownames(df)), size = 2, max.overlaps = Inf) +  # Adjust size for visibility
  scale_color_manual(values = c("#179f76", "#7570a7")) +  # Customize colors if needed
  theme_minimal() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right") + 
  theme(legend.text=element_text(size=10), legend.title=element_blank())
# guides(color = guide_legend(title = "Group"), shape = guide_legend(title = "Group"))  # Customize legend


p = p1/p2|p3

p = p+ plot_annotation(tag_levels = 'A')

ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/GSE141859_volcano_umap.png', p, width = 12, height = 8)




