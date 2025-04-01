##############################################################################################################
# Fig.1H Donut Plot for TCRβ CDR3 Distribution in TIL Products
##############################################################################################################
#library packages
library(ggplot2)
library(ggrepel)
library(ggnewscale)

TIL_product <- as.matrix(read.csv2(file.choose(), header = TRUE, sep = ","))
corresponding_tissue <- as.matrix(read.csv2(file.choose(), header = TRUE, sep = ","))

matched_CDR3 <- intersect(TIL_product[,1], corresponding_tissue[,1])

TIL_product <- as.data.frame(TIL_product)
TIL_product$count <- as.numeric(TIL_product$count)

# Calculate fraction and cumulative sum
TIL_product$fraction <- TIL_product$count / sum(TIL_product$count) * 100
TIL_product$ymax <- cumsum(TIL_product$fraction)
TIL_product$ymin <- c(0, head(TIL_product$ymax, n = -1))

# Filter matched CDR3s for annotation
label_df <- TIL_product[TIL_product$cdr3 %in% matched_CDR3, ]
label_df$ypos <- (label_df$ymin + label_df$ymax) / 2
label_df$x_label <- 4

TIL_product$matched.CDR3 <- ifelse(TIL_product$cdr3 %in% matched_CDR3, "Matched", "Non-Matched")

ggplot(TIL_product) +
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = 4.12, xmin = 4.05, fill = matched.CDR3), color = "white") +
  scale_fill_manual(values = c("Matched" = "black", "Non-Matched" = "white"))+
  new_scale_fill() + 
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = fraction),color = "white", size = 0.5) +
  scale_fill_gradient(low = "#FFF8DC", high = "#dd2b19", name = "Persentage", limits = c(0, 100))+
  coord_polar(theta = "y") +
  xlim(c(2, 5)) +  
  theme_void()


##############################################################################################################
# Fig.1I Dot Plot for TCRβ Clonity Metrics
##############################################################################################################
#library packages
library(ggplot2)
library(patchwork)

stat_dat <- read.csv2(file.choose(), header = TRUE, sep = ",")
stat_dat <- as.data.frame(stat_dat)
stat_dat$TIL.ID.. <- factor(stat_dat$TIL.ID.., levels = c("GBT23", "GBT28", "GBT29"))
stat_dat$D50 <- as.numeric(stat_dat$D50)
stat_dat$Entropy <- as.numeric(stat_dat$Entropy)
stat_dat$Unique.CDR3 <- as.numeric(stat_dat$Unique.CDR3)
stat_dat$Diversity.Index <- as.numeric(stat_dat$Diversity.Index)

p1 <- ggplot(stat_dat, aes(x = Species , y = CDR3, colour = TIL.ID.. )) +
  geom_point(size = 3)+
  scale_color_manual(values = c("GBT23" = "#ffad00", "GBT28" = "#F4651A", "GBT29" = "#DD0000")) +
  labs(x = "",
       y = "TCR count"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))+
  theme(legend.position = "none")

p2 <- ggplot(stat_dat, aes(x = Species, y = Unique.CDR3, colour = TIL.ID.. )) +
  geom_point(size = 3)+
  scale_color_manual(values = c("GBT23" = "#ffad00", "GBT28" = "#F4651A", "GBT29" = "#DD0000")) +
  labs(x = "",
       y = "Unique TCR count"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))+
  theme(legend.position = "none")

p3 <- ggplot(stat_dat, aes(x = Species, y = D50, colour = TIL.ID.. )) +
  geom_point(size = 3)+
  scale_color_manual(values = c("GBT23" = "#ffad00", "GBT28" = "#F4651A", "GBT29" = "#DD0000")) +
  labs(x = "",
       y = "D50 Index"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))+
  theme(legend.position = "none")

p4 <- ggplot(stat_dat, aes(x = Species, y = Diversity.Index, colour = TIL.ID.. )) +
  geom_point(size = 3)+
  scale_color_manual(values = c("GBT23" = "#ffad00", "GBT28" = "#F4651A", "GBT29" = "#DD0000")) +
  labs(x = "",
       y = "Diversity Index"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))+
  theme(legend.position = "none")

p5 <- ggplot(stat_dat, aes(x = Species, y = Entropy, colour = TIL.ID.. )) +
  geom_point(size = 3)+
  scale_color_manual(values = c("GBT23" = "#ffad00", "GBT28" = "#F4651A", "GBT29" = "#DD0000")) +
  labs(x = "",
       y = "Entropy"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))

p1+p2+p3+p4+p5 + plot_layout(ncol = 5)



##############################################################################################################
# Fig.2A-B Overlay Histogram for TCRβ CDR3 Distribution in Bulk Tumor Tissue
##############################################################################################################
#library packages
library(ggplot2)
library(ggrepel)
library(dplyr)

tumor_tissue <- as.matrix(read.csv2(file.choose(), header = TRUE, sep = ","))
corres_TIL <- as.matrix(read.csv2(file.choose(), header = TRUE, sep = ","))

tumorinput <- as.matrix(as.numeric(tumor_tissue[,2]))
row.names(tumorinput) <- tumor_tissue[,1]

tilinput <- as.matrix(as.numeric(corres_TIL[,2]))
row.names(tilinput) <- corres_TIL[,1]

matched_CDR3 <- intersect(tumor_tissue[,1], corres_TIL[,1])
cdr3mark <- which(rownames(tumorinput) %in% matched_CDR3)
labs <- rownames(tumorinput)[cdr3mark]

tumor_tissue <- as.data.frame(tumor_tissue)
tumor_tissue$count <- as.numeric(tumor_tissue$count)
tumor_tissue$rank <- c(1:nrow(tumor_tissue))

tumor_tissue <- tumor_tissue %>%
  mutate(
    percentage = count / sum(count) * 100,
    cumulative_count = cumsum(count), 
    cumulative_percentage = cumulative_count / sum(count) * 100
  )
label_tumor_tissue <- data.frame(x = cdr3mark, y = rep(0.5, length(matched_CDR3)), text = matched_CDR3)

ggplot(tumor_tissue, aes(x = rank, y = percentage)) +
  geom_vline(xintercept = cdr3mark, color="red", linetype="dashed", alpha=0.7) +
  geom_col(width = 1) +
  scale_x_continuous(breaks=c(1, nrow(tumor_tissue)))+
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Clone Ranked by Frequency",
    y = "Frequency (%)"
  ) +
  theme_classic() +
  theme(text = element_text(size = 14))



##############################################################################################################
# Fig.2C GSVA Heatmap for Immunophynotyping
##############################################################################################################
#library packages
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(data.table)
library(BiocParallel)
library(circlize)

# Create annotation matrix
Th1 <- c("IFNG", "TBX21", "CCL5", "CXCR3", "TBX21", "STAT11")
Th2 <- c("IL4", "IL5", "IL13", "GATA3", "STAT6", "IL1RL1")
Th9 <- c("IL9", "IL10", "SPI1", "STAT5A")
Th17 <- c("IL23R", "RORC", "STAT3")
Tfh <- c("BCL6", "ICOS", "CXCR5", "IL21", "CD38")
Treg <- c("FOXP3", "IKZF2", "IL10", "LAG3", "TGFB1", "IL2RA", "STAT5A")
CD8Tn <- c("CD8A", "CD8B", "LAG3", "TCF7", "SELL")
TcmTn <- c("CCR7", "CD8A", "SELL")
Exhaustion <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1", "PVR")
Bcell <- c("TNFRSF17", "BTK", "CD19", "CD22", "CD40", "CD79A", "CD79B", "CXCR5", "FCGR2B", "BCL6", "CD38", "CD27", "CD5")
NKcell <- c("KLRD1", "KLRB1","KLRC1", "NCR2", "NCR3", "KLRC1", "SLAMF7", "FCGR3A")
Myeloid <- c("MPO", "CD14", "CD68", "ITGAM")
Monocyte <- c("ITGAM", "CD14", "FCGR3A", "FCGR2A", "FCGR1A", "CD68", "MRC1", "FCGR2B", "MPO", "S100A8","PADI4")
DC <- c("ITGAX", "CD40", "CD83", "IL3RA", "HLA-DRA", "ITGAE", "S100A9")
Eosinophils <- c("EPX", "IL5", "IL4", "IL13")
Neutrophils <- c("MPO", "PADI4", "CD14", "FCGR3A")
Proinflammatory <- c("CCL5", "IFNG", "IL9", "IL12B", "IL13", "IL21", "TNF")
Antiinflammatory <- c("IL10", "IL4", "IL1R1", "IL5")
Pleiotropic <- c("IL5", "IL21", "IL23A", "IL27")
Activation <- c("TNFRSF9", "IL2RA", "CD38", "CD63", "CD69", "CD226", "CXCR3", "ICOS", "KLRB1", "KLRD1", "TFRC")
Cytotoxic <- c("GZMB", "PRF1", "FAS", "IFNG")

anno <- list(Th1 = Th1, Th2 = Th2, Th9 = Th9, Th17 = Th17, Tfh = Tfh, Treg = Treg, CD8Tn = CD8Tn, TcmTn = TcmTn, Exhaustion = Exhaustion, Bcell = Bcell, NKcell = NKcell, Myeloid = Myeloid, Monocyte = Monocyte, DC = DC, Eosinophils = Eosinophils, Neutrophils = Neutrophils, Proinflammatory = Proinflammatory, Antiinflammatory = Antiinflammatory, Pleiotropic = Pleiotropic, Activation = Activation, Cytotoxic = Cytotoxic)


#Load data
data <- read.csv(file = file.choose(), header = T, row.names = 1)
data[is.na(data)] <- 0
data <- data[,5:12]

# Run GSVA
params <- gsvaParam(as.matrix(data), anno, minSize = 1, maxSize = Inf, kcdf = "Gaussian", tau = 1, maxDiff = TRUE, absRanking = FALSE)
test <- gsva(params, verbose = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))

#Plot
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#424da7","#ffffff","#dd2b19"))
column_groups <- c("Non_TILs", "TILs", "Non_TILs", "TILs", "TILs", "TILs", "Non_TILs", "Non_TILs")
column_colors <- c("Non_TILs" = "#800000", "TILs" = "#00008B")
col_annotation <- rowAnnotation(Group = column_groups, col = list(Group = column_colors))

Heatmap(scale(t(test)), left_annotation = col_annotation, 
        name = "GSVA Score",rect_gp = gpar(col = "white", lwd = 2), col = col_fun, 
        clustering_method_rows = "complete", clustering_method_columns = "ward.D2", column_names_rot = 90)



##############################################################################################################
# Fig.S2 Boxplot for TCRβ Clonity Metrics
##############################################################################################################
#library packages
library(ggplot2)
library(ggsci)
library(patchwork)

stat_dat <- read.csv2(file.choose(), header = TRUE, sep = ",")
stat_dat <- as.data.frame(stat_dat)
stat_dat$D50 <- as.numeric(stat_dat$D50)
stat_dat$Entropy <- as.numeric(stat_dat$Entropy)
stat_dat$Unique.CDR3 <- as.numeric(stat_dat$Unique.CDR3)
stat_dat$CDR3 <- as.numeric(stat_dat$CDR3)
stat_dat$Diversity.Index <- as.numeric(stat_dat$Diversity.Index)


p1 <- ggplot(stat_dat, aes(x = Mark, y = CDR3)) +
  geom_boxplot(aes(color=Mark),linewidth=0.8,fill=NA)+
  scale_color_manual(values = c("TIL-" = "#800000", "TIL+" = "#00008B")) +
  labs(x = "TIL ID",
       y = "CDR3") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none")

p2 <- ggplot(stat_dat, aes(x = Mark, y = Unique.CDR3)) +
  geom_boxplot(aes(color=Mark),linewidth=0.8,fill=NA)+
  scale_color_manual(values = c("TIL-" = "#800000", "TIL+" = "#00008B")) +
  labs(x = "TIL ID",
       y = "Unique.CDR3") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none")

p3 <- ggplot(stat_dat, aes(x = Mark, y = D50)) +
  geom_boxplot(aes(color=Mark),linewidth=0.8,fill=NA)+
  scale_color_manual(values = c("TIL-" = "#800000", "TIL+" = "#00008B")) +
  labs(x = "TIL ID",
       y = "D50") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none")

p4 <- ggplot(stat_dat, aes(x = Mark, y = Diversity.Index)) +
  geom_boxplot(aes(color=Mark),linewidth=0.8,fill=NA)+
  scale_color_manual(values = c("TIL-" = "#800000", "TIL+" = "#00008B")) +
  labs(x = "TIL ID",
       y = "Diversity.Index") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none")

p5 <- ggplot(stat_dat, aes(x = Mark, y = Entropy)) +
  geom_boxplot(aes(color=Mark),linewidth=0.8,fill=NA)+
  scale_color_manual(values = c("TIL-" = "#800000", "TIL+" = "#00008B")) +
  labs(x = "TIL ID",
       y = "Entropy") +
  theme_classic() +
  theme(axis.title.x = element_blank()) 

p1+p2+p3+p4+p5 + plot_layout(ncol = 5)