###### Part 1: Venn diagrams show common genes between DEGs in GSE35088 and and target genes of DE-miRNAs from GSE36946 #####

# Load necessary libraries for data manipulation, miRNA-target interaction analysis, and plotting
library(GEOquery)   # For downloading and accessing GEO datasets
library(multiMiR)   # For querying miRNA-target interactions databases
library(ggvenn)     # For creating Venn diagrams to visualize overlaps
library(patchwork)

# Read differential expression genes (DEGs) file for GSE36946 with headers, separated by tabs
demirna = read.table("data/DEGs_GSE95855.txt", header=TRUE, sep="\t")

# Define a vector of miRNAs involved in the study
mirnas = demirna$ID
mirnas

# Define a corresponding vector of precursors for the miRNAs
precursor = demirna$ID
precursor

# Define a vector indicating whether the miRNAs are up or down-regulated
upDown = ifelse(demirna$logFC<0, "down", "up")
upDown

# Combine miRNA, precursor, and expression vectors into a data frame
demirnas = cbind.data.frame(miRNA = mirnas, precursor = precursor, expr = upDown)
demirnas

# Retrieve target genes for the specified miRNAs from the multiMiR database
getTargets = get_multimir(org = "rno", mirna = demirnas$precursor, summary = T, table = "predicted")
getTargets

# Read DEGs data for GSE35088 and store it in a variable
degs_gse35088_data = read.table("data/DEGs_GSE35088.txt", header = TRUE, sep = "\t")
# Display dimensions and the first few rows of the dataset to understand its structure
dim(degs_gse35088_data)
head(degs_gse35088_data)

# Find up-regulated and down-regulated common DEGs based on log fold change (logFC)
upCommonDEGs = degs_gse35088_data[degs_gse35088_data$logFC > 0, "Symbol"]
upCommonDEGs

downCommonDEGs = degs_gse35088_data[degs_gse35088_data$logFC < 0, "Symbol"]
downCommonDEGs

# Filter miRNA targets from the retrieved target list based on their expression direction
upMiRNAtargets = getTargets@data[getTargets@data$mature_mirna_id %in% demirnas[demirnas$expr == "up", "precursor"], "target_symbol"]
length(upMiRNAtargets)
downMiRNAtargets = getTargets@data[getTargets@data$mature_mirna_id %in% demirnas[demirnas$expr == "down", "precursor"], "target_symbol"]
length(downMiRNAtargets)


# Identify and save up-regulated DEGs that are also targeted by up-regulated miRNAs
updegs = degs_gse35088_data[degs_gse35088_data$logFC > 0,]
upDEGs = updegs[updegs$Symbol %in% upMiRNAtargets,]
dim(upDEGs)
write.table(upDEGs, "data/common_upDEGs.txt", quote = FALSE, row.names = FALSE, sep = "\t")

upDEGs = upDEGs[abs(upDEGs$logFC) > 0.9,]
upDEGs

# Read top genes data from GSE32453 for further analysis
# topgenes_gse32453 = read.table("data/TopGenes_GSE32453.txt", 
#                                header=TRUE, sep="\t")
# 
# # Retrieve platform annotation data for GPL6104 to map gene symbols
# gpl6104 = getGEO("GPL6104")
# annotation = gpl6104@dataTable@table %>% dplyr::select(ID, ILMN_Gene)
# 
# # Merge top genes data with annotation data based on ID
# topgenes_gse32453 = left_join(topgenes_gse32453, annotation, by="ID")
# 
# # Find up-regulated DEGs in top genes list and annotate them
# upDEGs_test = topgenes_gse32453[topgenes_gse32453$ID %in% annotation[annotation$ILMN_Gene %in% updegs[updegs$ID %in% upMiRNAtargets,"ID"],"ID"],]

# Generate a list for Venn diagram comparison of DEGs and DE-miRNA Targets
x = list(
  DEGs = upCommonDEGs,
  "DE-miRNA Targets" = upMiRNAtargets
)

# Create a Venn diagram for up-regulated DEGs and miRNA targets
p1 = ggvenn(
  x, 
  fill_color = c("#FF2C21", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 6, text_size = 6
)

# Exclude a specific gene from down-regulated DEGs for further analysis
downdegs = degs_gse35088_data[degs_gse35088_data$logFC < 0,]
dim(downdegs)
downdegs
# Save down-regulated DEGs that are targeted by down-regulated miRNAs
write.table(downdegs[downdegs$Symbol %in% downMiRNAtargets,], "data/common_downDEGs.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# # Find down-regulated DEGs in top genes list and annotate them
# downDEGs_test = topgenes_gse32453[topgenes_gse32453$ID %in% annotation[annotation$ILMN_Gene %in% downdegs[downdegs$ID %in% downMiRNAtargets,"ID"],"ID"],]
# downDEGs_test$ILMN_Gene

# Prepare a list for Venn diagram of down-regulated DEGs excluding a specific gene
x = list(
  DEGs = downCommonDEGs,
  "DE-miRNA Targets" = downMiRNAtargets
)

# Identify overlaps between down-regulated DEGs and miRNA targets, excluding a specific gene
downDEGs = downdegs[downdegs$Symbol %in% downCommonDEGs[downCommonDEGs %in% downMiRNAtargets], ]
dim(downDEGs)
downDEGs = downDEGs[abs(downDEGs$logFC) > 0.9,]
downDEGs

# Create a Venn diagram for down-regulated DEGs and miRNA targets
p2 = ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 6, text_size = 6
)

# Combine both Venn diagrams for up and down-regulated genes and targets
p = p1 + p2
p = p + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = '18'))

# Save the combined Venn diagram plot to a file
ggsave('~/Documents/Studies/AMI_miRNA/Figures/Figure4.png', p, width = 12, height = 8)


# ###### Part 2: Expression patterns in GSE32453 of common miRNAs with targets from the DEGs #####
# 
# # Load necessary libraries for statistical analysis and enhanced plotting
# library(rstatix)    # Helper functions for statistical tests and summaries
# library(ggpubr)     # For enhancing 'ggplot2' plots for publication quality
# 
# # Read expression data for the study GSE32453
# exprs_gse32453 = read.table("data/exprs_GSE32453.txt", header = TRUE, sep = "\t")
# # Display the dimensions of the expression data
# dim(exprs_gse32453)
# # Preview a specific subset of the data
# exprs_gse32453[1:5,13131:13132]
# 
# # Prepare data for upregulated DEGs for analysis
# # Select relevant columns based on upregulated DEGs and group them
# upDEGs_test_data = exprs_gse32453 %>% dplyr::select(upDEGs_test$ID, group)
# # Rename the first two columns based on gene symbols for clarity
# colnames(upDEGs_test_data)[1:2] = upDEGs_test$ILMN_Gene
# 
# # Convert the wide data to a long format for easier analysis
# upDEGs_test_data2 <- upDEGs_test_data %>%
#   gather(key = "Expression", value = "Value", upDEGs_test$ILMN_Gene) %>%
#   convert_as_factor(Expression) %>% na.omit()
# 
# # Perform t-tests on the long-format data, adjust p-values for multiple testing using Bonferroni method, and add significance labels
# stat.test <- upDEGs_test_data2 %>%
#   group_by(Expression) %>%
#   t_test(Value ~ group) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance()
# 
# # Plot the results of the upregulated DEGs analysis
# myplot <- upDEGs_test_data2 %>%
#   drop_na() %>%
#   ggplot(aes(x=Expression, y=Value, color=group)) +
#   geom_boxplot() +
#   geom_point(position=position_jitterdodge(), alpha = 0.3) +
#   scale_color_manual(values=c("#1B9E77", "#7570B3")) +
#   ylab("Target Gene Expression") +
#   xlab("") +
#   theme_bw() +
#   theme(text = element_text(size = 20)) +
#   theme(legend.position="top", legend.title = element_blank())
# 
# # Add statistical annotations to the plot
# stat.test <- stat.test %>% add_xy_position(x = "Expression")
# p1 = myplot + stat_pvalue_manual(stat.test, label = "p")
# 
# # Repeat the process for downregulated DEGs
# downDEGs_test_data = exprs_gse32453 %>% dplyr::select(downDEGs_test$ID, group)
# colnames(downDEGs_test_data)[1:6] = downDEGs_test$ILMN_Gene
# 
# downDEGs_test_data2 <- downDEGs_test_data %>%
#   gather(key = "Expression", value = "Value", downDEGs_test$ILMN_Gene) %>%
#   convert_as_factor(Expression) %>% na.omit()
# 
# stat.test <- downDEGs_test_data2 %>%
#   group_by(Expression) %>%
#   t_test(Value ~ group) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance()
# 
# # Plot the results of the downregulated DEGs analysis
# myplot <- downDEGs_test_data2 %>%
#   drop_na() %>%
#   ggplot(aes(x=Expression, y=Value, color=group)) +
#   geom_boxplot() +
#   geom_point(position=position_jitterdodge(), alpha = 0.3) +
#   scale_color_manual(values=c("#1B9E77", "#7570B3")) +
#   ylab("Target Gene Expression") +
#   xlab("") +
#   theme_bw() +
#   theme(text = element_text(size = 20)) +
#   theme(legend.position="top", legend.title = element_blank())
# 
# stat.test <- stat.test %>% add_xy_position(x = "Expression")
# p2 = myplot + stat_pvalue_manual(stat.test, label = "p")
# 
# # Combine plots for both upregulated and downregulated DEGs
# p = p1 + p2
# p = p + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = '18'))
# 
# # Save the combined plot to a file
# ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/clinical_validation.png', p, width = 14, height = 8)
# 
