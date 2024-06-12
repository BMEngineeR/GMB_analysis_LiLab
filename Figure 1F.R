library(fgsea)
library(pheatmap)

load("TCGA_cell_2013_meta.Rdata")
load("TCGA_cell_2013_exp.Rdata")
angiogenesis_geneset <- read.table("../pathway_GSEA/ANGIOGENESIS.v2022.1.Hs.grp",header = T)
angiogenesis_geneset <- angiogenesis_geneset$ANGIOGENESIS
myeloid_compartment_geneset <- read.table("../pathway_GSEA/BROWN_MYELOID_CELL_DEVELOPMENT_UP.v2022.1.Hs.grp",header = T)
myeloid_compartment_geneset <- myeloid_compartment_geneset$BROWN_MYELOID_CELL_DEVELOPMENT_UP
Tcell_signal_geneset <- read.table("../pathway_GSEA/KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY.v2022.1.Hs.grp",header = T)
Tcell_signal_geneset <- Tcell_signal_geneset$KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY
Tcell_exhaustion_geneset <- read.table("../pathway_GSEA/GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN.v2022.1.Hs.grp",header = T)
Tcell_exhaustion_geneset <- Tcell_exhaustion_geneset$GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN

GARP_list <- c("LRRC32", "ITGB6", "ITGB8", "ITGAV",
               "ITGA2B", "SELP", "GP9", "GP1BA","F2")
pathway_list <- list(angiogenesis_geneset = angiogenesis_geneset,
                     myeloid_compartment_geneset = myeloid_compartment_geneset,
                     Tcell_signal_geneset = Tcell_signal_geneset,
                     Tcell_exhaustion_geneset = Tcell_exhaustion_geneset,
                     GARP_geneset = GARP_list)

# 

TCGA_cell_2013_meta <- TCGA_cell_2013_meta[TCGA_cell_2013_meta$tutmor_type != "G-CIMP" & !is.na(TCGA_cell_2013_meta$tutmor_type),]
TCGA_cell_2013_meta$tutmor_type_binary <- ifelse(TCGA_cell_2013_meta$tutmor_type == "Mesenchymal","Mesenchymal","non_Mesenchymal")

table(TCGA_cell_2013_meta$tutmor_type)
TCGA_cell_2013_exp <- TCGA_cell_2013_exp[,TCGA_cell_2013_meta$ID]
TCGA_cell_2013_exp <- TCGA_cell_2013_exp[-which(rownames(TCGA_cell_2013_exp)==""),]
intersect(rownames(TCGA_cell_2013_exp),angiogenesis_geneset)
# DESEQ2
dim(TCGA_cell_2013_exp)
gene_id <- cbind('Probe Set ID' = rownames(TCGA_cell_2013_exp), "Gene Symbol" = rownames(TCGA_cell_2013_exp), "Gene Title" = "na")


# heatmap
data_set <- myeloid_compartment_geneset

TCGA_cell_2013_meta_df <-TCGA_cell_2013_meta[order(TCGA_cell_2013_meta$tutmor_type_binary),]
TCGA_cell_2013_meta_df_angiogenesis <- TCGA_cell_2013_exp[data_set,]
TCGA_cell_2013_meta_df_angiogenesis <- TCGA_cell_2013_meta_df_angiogenesis[rowSums(TCGA_cell_2013_meta_df_angiogenesis) > 0,]
TCGA_cell_2013_meta_df_angiogenesis <- TCGA_cell_2013_meta_df_angiogenesis[!is.na(TCGA_cell_2013_meta_df_angiogenesis$`TCGA-02-0047`),]
TCGA_cell_2013_meta_df_angiogenesis <- TCGA_cell_2013_meta_df_angiogenesis[,rownames(TCGA_cell_2013_meta_df)]
TCGA_cell_2013_meta_df_angiogenesis_log <- log1p(TCGA_cell_2013_meta_df_angiogenesis)
Myeloid_rank <- read.csv("BROWN_MYELOID_CELL_DEVELOPMENT_UP.csv")

my_sample_col <- data.frame(sample = TCGA_cell_2013_meta_df$tutmor_type_binary)
rownames(my_sample_col) <- TCGA_cell_2013_meta_df$ID
pheatmap(mat = TCGA_cell_2013_meta_df_angiogenesis_log[Myeloid_rank$GENE.SYMBOL[c(1:20,(nrow(Myeloid_rank)-10):nrow(Myeloid_rank))],-24],
         scale = "row",cluster_cols = F,cluster_rows = F,annotation_col = my_sample_col,
         breaks = c(seq(-2,2,length.out = 200)),color = colorRampPalette(c("blue","white","red"))(n=200),labels_col = NA,fontsize = 18)

my.comparison <- list(c("Mesenchymal","non_Mesenchymal"))

# geneset
TCGA_cell_2013_meta_GARP <-TCGA_cell_2013_meta[,c("tutmor_type_binary","GARP_geneset")]
TCGA_cell_2013_meta_GARP <- TCGA_cell_2013_meta_GARP[!is.na(TCGA_cell_2013_meta_GARP$GARP_geneset),]
my.comparison <- list(c("Mesenchymal","non_Mesenchymal"))


TCGA_cell_2013_meta_GARP$tutmor_type_binary <- factor(TCGA_cell_2013_meta_GARP$tutmor_type_binary,levels = c("Mesenchymal","non_Mesenchymal"))
ggboxplot(TCGA_cell_2013_meta_GARP, x = "tutmor_type_binary", y = "GARP_geneset",color = "tutmor_type_binary",palette =c("#fe7e83", "#8b8b95"),add = "jitter")+
  stat_compare_means(comparisons = my.comparison)  



