rm(list=ls())
library(Seurat)
library(readxl)
library(tidydr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(parallel)
library(parallelMap)
parallelStartSocket(cpus = 12)
load("tumor_merge_manual_copykat.Rda")

pre_var <- read.delim("D:/nb/NBL2024.3.4/bulk/GSE49710/2.model construction/ML/pre_var.txt")
feature=list(feature.score = pre_var$x)
sco <- AddModuleScore(sco,feature)
names(sco@meta.data)[15] <- "feature.score"


ggplot(sco@meta.data,aes(Major.manual.subtype.copykat,feature.score))+
  geom_boxplot(aes(fill=Major.manual.subtype.copykat),width=.5,
               outlier.size = .2)+
  scale_fill_manual(values =  c("#7fb961","#356d67","#5066a1",
                                "#76afda","#dcf2ff","#ffe788",
                                "#ffc556",
                                "#e8743c",
                                "#b20000",
                                "#a14462",
                                "#cca69c",
                                "#7d4444"))+
  theme_classic()+ylab("Signature score")+
  coord_flip()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold",color = "black",size = 12),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,
                                   vjust = 0,hjust = 0.5))

ggsave(filename ="1.annotation/Major.manual.subtype.copykat boxplot.pdf",width = 3.5,height =5.5)


data <- as.data.frame(sco[["tsne"]]@cell.embeddings)
data$feature.score <- sco$feature.score

ggplot(data, aes(x = tSNE_1, 
                 y = tSNE_2, 
                 color = feature.score))+
  geom_point(size = .5,shape=16) +
  scale_color_viridis_c(begin = .2)+
  theme_dr()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust=.5,size=12),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(face = "bold.italic",size=10))+
  xlab("tSNE1")+ylab("tSNE1")
ggsave(filename ="1.annotation/Major.manual.subtype.copykat score tsne.pdf",width = 3.8,height = 3)



data_pro <- sco@meta.data[,c("Major.manual.subtype","copykat.pred")]

data_pro$copykat.pred <- ifelse(data_pro$copykat.pred=="aneuploid","malignant","non-malignant")

data_plot=table(data_pro$Major.manual.subtype,data_pro$copykat.pred) %>%
  as.data.frame()


ggplot(data_plot,aes(Var2,Freq,fill=Var1))+
  geom_col(position = "fill",width = .5)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values =  c("#7fb961","#356d67","#5066a1",
                                "#76afda","#dcf2ff","#ffe788",
                                "#ffc556",
                                "#e8743c",
                                "#b20000",
                                "#a14462",
                                "#cca69c",
                                "#7d4444"))+
  xlab("")+ylab("Proportion")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12,color = "black"),
        axis.text = element_text(colour = "black" ),
        axis.title.y = element_text(colour = "black",face = "bold",size = 12 ),
        legend.position = "none")


ggsave(filename ="1.annotation/Major.manual.subtype.copykat score Proportion.pdf",width = 2.5,height = 3)





#Malignant----

Mal_sce <- subset(sco,idents ="Malignant" )
Mal_sce <- NormalizeData(Mal_sce) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)


Mal_sce <- FindClusters(Mal_sce, resolution = .3)
Mal_sce$Subcluster <- paste0("Mal",Mal_sce$RNA_snn_res.0.3)

DimPlot(Mal_sce, reduction = "tsne" ,group.by = "Subcluster")+
  scale_color_manual(values =c("#ffc556","#76afda","#4c9568") )+
  theme_dr()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(face = "bold.italic",size=10))+
  xlab("tSNE1")+ylab("tSNE1")

ggsave(filename ="1.annotation/Mal cluster tsne.pdf",width = 3.8,height = 3)


FeaturePlot(Mal_sce, features = "PLK1",reduction = "tsne" )+
  scale_color_gradientn(colours =viridis(32,begin =0.3,end=1),name = "PLK1")+
  theme_dr()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14),
        axis.line = element_line(linewidth = 1.5),
        legend.text = element_text(size = 10),
        axis.title = element_text(face = "bold.italic",size=15))+
  xlab("tSNE1")+ylab("tSNE1")

ggsave(filename ="1.annotation/Mal cluster PLK1 tsne.pdf",width =4,height = 4)






library(viridis)
genes <- rev(c("RFC4","NCAPG","PLK1","MAD2L1","MCM2","PRC1","SPC24","RAD51","POLE2","SMC4"))
DotPlot(sco, features = genes,group.by = "Subcelltype")+
  scale_color_gradientn(colours =viridis(32,begin =0.3,end=1),name = "PLK1")+
  # theme_dr()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=14),
        axis.title = element_blank())

ggsave(filename ="1.annotation/Mal cluster genes DotPlot.pdf",width =4.5,height =3.5)


PLK1_signaling <- read_xlsx("PLK1 signaling.xlsx")
feature=list(PLK1_signaling = PLK1_signaling$Gene)
Mal_sce <- AddModuleScore(Mal_sce,feature)
names(Mal_sce@meta.data)[26] <- "PLK1 signaling events"


data <- as.data.frame(Mal_sce[["tsne"]]@cell.embeddings)
data$feature.score <- Mal_sce$`PLK1 signaling events`

ggplot(data, aes(x = tSNE_1, 
                 y = tSNE_2, 
                 color = feature.score))+
  geom_point(size = 1.2,shape=16) +
  scale_color_viridis_c(begin = .2,name = "PLK1 signaling events")+
  theme_dr()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust=.5,size=12),
        legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        axis.title = element_text(face = "bold.italic",size=10))+
  xlab("tSNE1")+ylab("tSNE1")
ggsave(filename ="1.annotation/Mal cluster PLK1 signaling events tsne.pdf",width = 4.7,height = 4)


Mal_sig <- read.delim("D:/nb/NBL2024.3.4/GSE192906/Mal_genes.txt")
Mal1_sig <- list(Mal1_sig=Mal_sig[Mal_sig$cluster=="1","gene"])
Mal_sce <- AddModuleScore(Mal_sce,Mal1_sig)
names(Mal_sce@meta.data)[26] <- "Mal1_sig"

PLK1_gene <- GetAssayData(Mal_sce,slot = "data") %>% 
  as.matrix() %>% .["PLK1",,drop=F] %>%
  t() %>% as.data.frame()

Mal_sce <- AddMetaData(Mal_sce,PLK1_gene)

data <- Mal_sce@meta.data[,c(16,26:27)] %>% .[.$Subcluster == "Mal1",]


ggscatter(data,x="PLK1",y="Mal1_sig",
          size = 1.5,color ="#76afda",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#5066a1", fill = "gray"), # Customize reg. line
          conf.int = TRUE)+
  stat_cor(method = "spearman",size=5)+
  geom_rug(color ="#76afda")+
  labs(y='Mal1 score',x='PLK1 expression level')+
  theme_classic(base_line_size = .7)+
  theme(legend.position = c(0.12,0.3),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'black',face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'black',face='bold'))


ggsave(filename = paste0("D:/nb/NBL2024.3.4/GSE192906/1.annotation/tumor Cell PLK1",names(ee)[i],".pdf"),width = 4,height =4)



#----
load("Mal_markers.rda")
library(IOBR)
library(org.Hs.eg.db)
library(msigdbr)
m_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
m_C2_kegg <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)
m_C2_Reactome <- msigdbr(species = "Homo sapiens", category = "C2",subcategory ="REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)
m_C5_BP <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- rbind(m_h,m_C2_kegg) %>% rbind(.,m_C2_Reactome) %>% rbind(.,m_C5_BP)
head(m_t2g)


gsea_res_list <- list()
result_data <- data.frame()
for (i in unique(Mal_markers$cluster)) {
  
  data <- subset(Mal_markers,cluster==i)
  
  entrezid <- bitr(data$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  
  list <- merge(entrezid,data,by.x=1,by.y=7)
  
  genelist = list$avg_log2FC
  names(genelist) = list$ENTREZID
  genelist <- sort(genelist,decreasing = T)
  
  
  gsea_res <- GSEA(genelist, 
                   TERM2GENE = m_t2g,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH"
  )
  
  gsea_res_list[[i]] <- gsea_res
  
  result = gsea_res@result %>% mutate(cluster = i)
  result_data <- rbind(result_data,result)
}


result_data1 <- result_data %>% filter(p.adjust < 0.05 & NES >0)

Cluster0_pathways <- c(#Cluster 0
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN",
  "GOBP_LEUKOCYTE_DEGRANULATION",
  "HALLMARK_KRAS_SIGNALING_DN",
  "GOBP_MAST_CELL_ACTIVATION",
  "KEGG_AUTOIMMUNE_THYROID_DISEASE",
  "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION")

#Cluster 1
Cluster1_pathways <- c("HALLMARK_MYC_TARGETS_V1",
                       "HALLMARK_E2F_TARGETS",
                       "HALLMARK_MTORC1_SIGNALING",
                       "REACTOME_DNA_REPLICATION",
                       "REACTOME_G2_M_CHECKPOINTS",
                       "REACTOME_CELL_CYCLE_CHECKPOINTS",
                       "HALLMARK_MYC_TARGETS_V2",
                       "REACTOME_S_PHASE",
                       "REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION")

#Cluster 2
Cluster2_pathways <- c("GOBP_NEURON_PROJECTION_EXTENSION",
                       "GOBP_NEUROTRANSMITTER_SECRETION",
                       "GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT",
                       "GOBP_AXON_DEVELOPMENT",
                       "GOBP_SYNAPTIC_SIGNALING",
                       "GOBP_AXONAL_TRANSPORT",
                       "GOBP_NEURON_DEVELOPMENT",
                       "REACTOME_NCAM_SIGNALING_FOR_NEURITE_OUT_GROWTH",
                       "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY",
                       "REACTOME_NEURONAL_SYSTEM")

df <- result_data1[result_data1$Description %in% Cluster0_pathways & result_data1$cluster=="0",] %>%
  rbind(.,result_data1[result_data1$Description %in% Cluster1_pathways & result_data1$cluster=="1",] )%>%
  rbind(.,result_data1[result_data1$Description %in% Cluster2_pathways & result_data1$cluster=="2",])
df$Description <- str_remove_all(df$Description,"HALLMARK_|GOBP_|REACTOME_")
df$Description <- str_replace_all(df$Description,"_"," ")
df$Description <- factor(df$Description, levels = rev(df$Description))


ggplot() +
  geom_bar(data = df,
           aes(x = NES, y = Description,fill=cluster),
           width = .8, 
           stat = 'identity') +
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  geom_text(data = df,
            aes(x = 0.1, 
                y = Description,
                label = Description),
            size = 3,fontface = 'italic', 
            hjust = 0) +#
  labs(x = 'NES',
       y = 'Mal2                            Mal1                               Mal0', 
       title = 'Functional analysis of malignant cluster')+
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 12, face = 'bold'),
    axis.title = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = 'bold' ),
    axis.text = element_text(size = 10,color="black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, face = 'bold.italic'),
    axis.ticks.y = element_blank())+
  scale_fill_manual(values = c("#ffc556","#76afda","#4c9568"))

ggsave(filename ="1.annotation/Mal cluster Functional analysis.pdf",width = 3.6,height = 6)


#NB signature----
sign <- read_xlsx("NB signature.xlsx",sheet =3)
sign <- split(as.matrix(sign[,1]),sign[,2])
Mal_sce <- AddModuleScore(Mal_sce,features = sign,name = names(sign))
names(Mal_sce@meta.data)

library(gghalves)
ggplot(Mal_sce@meta.data,aes(Subcluster,feature.score,fill=Subcluster))+
  geom_half_violin(side = "l",color="white",width=0.5)+
  geom_half_point(side = "r",size=1,aes(color=Subcluster),range_scale = 0.5
  )+
  scale_fill_manual(values = c("#ffc556","#76afda","#4c9568"))+
  scale_color_manual(values = c("#ffc556","#76afda","#4c9568"))+
  stat_compare_means(size=5,label.x.npc = .05,label.y.npc = .95)+
  labs(y="Signature score")+
  theme_classic()+
  theme(
    legend.position = 'none',
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = 'bold' ),
    axis.text.x = element_text(size = 16,color="black", face = 'bold'),
    axis.text.y = element_text(size = 10,color="black")
  )


ggsave(filename ="1.annotation/Mal cluster Signature.pdf",width = 4,height = 3.5)




Mal_sce.marker.sig <- Mal_markers %>% 
  mutate(Ratio = round(pct.1/pct.2,3),
         pct.fc = pct.1-pct.2) %>%
  filter(p_val_adj < 0.05, 
         pct.1 > 0.7,
         pct.2 < 0.2,
         pct.fc > 0.5,
         avg_log2FC > .5)
write.table(Mal_sce.marker.sig,file = "Mal_genes.txt",row.names = T,sep = "\t",quote = F)


VlnPlot(Mal_sce,features = "PLK1")
DotPlot(Mal_sce,features = "PLK1")


# CytoTRACE2-----
library (CytoTRACE2)  #loading
cytotrace2_result <- cytotrace2(Mal_sce, is_seurat = TRUE, 
                                slot_type = "counts", 
                                species = 'human')
annotation <- data.frame(phenotype = Mal_sce@active.ident) %>% 
  set_rownames(., colnames(Mal_sce))
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)
plots$CytoTRACE2_Boxplot_byPheno
plots$CytoTRACE2_Potency_UMAP

Mal_sce <- AddMetaData(Mal_sce, cytotrace2_result@meta.data[,21:25])


DimPlot(Mal_sce, reduction = "tsne" ,group.by = "CytoTRACE2_Potency")+
  scale_color_manual(values = rev(c("#9E0142", "#F46D43", "#FEE08B", 
                                    "#E6F598", "#66C2A5", "#5E4FA2")))+
  theme_dr()+
  labs(title="CytoTRACE")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust=.5,size=12),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(face = "bold.italic",size=10))+
  xlab("tSNE1")+ylab("tSNE1")
ggsave(filename ="1.annotation/CytoTRACE DimPlot.pdf",width = 4,height = 3)



FeaturePlot(Mal_sce, reduction = "tsne" ,#pt.size = .5,
            features = "CytoTRACE2_Relative")+
  scale_color_gradientn(colors =  rev(c("#9E0142", "#F46D43", "#FEE08B", 
                                        "#E6F598", "#66C2A5", "#5E4FA2"))
  )+
  theme_dr()+
  labs(title="CytoTRACE")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust=.5,size=14),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.width = unit(4,"mm"),
        axis.title = element_text(face = "bold.italic",size=11))+
  xlab("tSNE1")+ylab("tSNE1")
ggsave(filename ="1.annotation/CytoTRACE FeaturePlot.pdf",width = 4,height = 4)




# SingleR----
library(SingleR)
library(BiocParallel)

count_matrix <- GetAssayData(Mal_sce,slot = "count")

ref <- readRDS("neuroblastoma-v1.2.0/data_raw/adrmed/adrenal_medulla_Seurat.RDS")
ref <- updateObject(ref)

results <- SingleR(
  test = count_matrix,
  ref = ref$RNA@data,
  labels = Idents(ref),
  de.method = "wilcox",
  BPPARAM = MulticoreParam(workers = 72, progressbar = TRUE, RNGseed = 42)
)

df <- as_tibble(results, rownames = "cell")

meta <- Mal_sce@meta.data[,"Subcluster",drop=F]
data <- merge(meta,df,by.x=0,by.y=1)
table(data$labels,data$Subcluster)
data$Subcluster <- as.character(data$Subcluster)


vis_data <- 
  data %>% 
  group_by(Subcluster) %>% 
  dplyr::count(labels) %>%
  mutate(n_rel = n / sum(n) * 100) %>%
  ungroup()

colors <- c("#921813","#be202e","#be6867","#ef845f","#33ac75","#006e3b","#686c58","#aacedc","#244a94","#303d63","gray50")
names(colors) <- c("late SCPs","SCPs","cycling SCPs", "Bridge",                     
                   "connecting Chromaffin cells", "Chromaffin cells",           
                   "late Chromaffin cells","cycling Neuroblasts",        
                   "Neuroblasts","late Neuroblasts","other")
# make plot
ggplot(vis_data, aes(labels, n_rel)) +
  geom_hline(yintercept = 0, size = .4) +
  geom_col(aes(fill = labels),width = .8, show.legend = T) +
  scale_x_discrete("adrenal medullary cell type") +
  labs(fill ="",y="Percentage of cells")+
  scale_fill_manual(values = colors) +
  coord_flip() +
  facet_wrap(vars(Subcluster), nrow = 1) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold",size = 15),
    axis.text.x =  element_text(size = 12,color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold",size = 14),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 14),
    legend.position = "none"
  )

ggsave(filename ="1.annotation/Percentage of cells.pdf",width = 6.5,height = 4)


save(Mal_sce,file = "Mal_sce.Rda")


#----
library(SCENIC)

setwd("SCENIC")
library(tidyverse)
library(Seurat) 
library(SCENIC)
library(parallel)
library(parallelMap)
parallelStartSocket(cpus = detectCores())
load("../Mal_sce.Rda")

exprMat  <-  as.matrix(Mal_sce@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4]
cellInfo <-  Mal_sce@meta.data[,c("Subcluster","nCount_RNA","nFeature_RNA")]
colnames(cellInfo)=c('CellType','nUMI', 'nGene' )

saveRDS(cellInfo,file = "cellInfo.RDS")


hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather',  '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather')

db_path <- '.'

scenicOptions <- initializeScenic(org = 'hgnc', dbDir = db_path, 
                                  dbs = hg38_dbs, 
                                  dbIndexCol = "motifs",
                                  nCores=1)

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)


runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

save(exprMat,scenicOptions,file = "runGenie3.Rda")


### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"),
                                            dbIndexCol = "motifs") # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings


export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 




#Slingshot-----

library(slingshot)
library(SingleCellExperiment)
library(qs)
library(Seurat)
library(tidyverse)
library(RColorBrewer)


counts <- Mal_sce@assays$RNA@counts
sim <- SingleCellExperiment(assays = List(counts = counts)) 

tsne = Mal_sce@reductions$tsne@cell.embeddings
colnames(tsne) = c('TSNE-1', 'TSNE-2')
reducedDims(sim) = SimpleList(TSNE = tsne)

# metadata
meta = Mal_sce@meta.data

colData(sim)$sampleId = meta$orig.ident
colData(sim)$celltype = meta$Subcluster


sim <- slingshot(sim, 
                 clusterLabels = 'celltype', 
                 reducedDim = 'TSNE',  
                 start.clus= "Mal1",  
                 end.clus = NULL    
)     

sling <- sim@colData[,"slingPseudotime_1",drop=F] %>% as.data.frame() 
Mal_sce <- AddMetaData(Mal_sce,sling)


#monocle-----
library(monocle)

Mal_mon <- Mal_sce

expr_matrix <- as(as.matrix(Mal_mon@assays$RNA@counts), 'sparseMatrix')

p_data <- Mal_mon@meta.data 
p_data$celltype <- Mal_mon$Subcluster  
f_data <- data.frame(gene_short_name = row.names(Mal_mon),row.names = row.names(Mal_mon))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,  phenoData = pd,featureData = fd,
                      lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) 
expressed_genes<- VariableFeatures(Mal_mon)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype",cores=12)
head(diff)

deg <- subset(diff, qval < 0.01)           
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)


ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene) 
plot_ordering_genes(cds)
ordergene <-row.names(deg)[order(deg$qval)]

cds <- reduceDimension(cds,max_components = 2,method = 'DDRTree')
cds <- orderCells(cds,reverse = T)


library(monocle)
library(viridis)
plot_cell_trajectory(cds,color_by="Pseudotime",show_backbone=TRUE,
                     cell_size = 1, cell_link_size = 1)+
  theme_dr()+
  scale_color_gradientn(colours =viridis(32,begin =0.3,end=1))+ #
  theme(legend.position=c(.5,.9),
        legend.direction = 'horizontal',
        legend.title = element_text(hjust = 0.5 ,size = 11,color = 'black',face = 'bold'),
        legend.text = element_text(hjust = 0.5 ,size = 10,color = 'black'),
        axis.title = element_text(size = 12,face = 'bold',colour = 'black'), #darkred
        axis.text = element_blank(),#element_text(size = 12,colour = 'black'),
        panel.grid = element_blank()
  )
ggsave(filename = "1.annotation/tumor Cell colored by Pseudotime.pdf",width = 4,height =4)
# dev.off()

plot_cell_trajectory(cds,color_by="slingPseudotime_1",show_backbone=TRUE,
                     cell_size = 1, cell_link_size = 1)+
  theme_dr()+
  scale_color_gradientn(colours =viridis(32,begin =0.3,end=1),name = "slingPseudotime")+ #
  theme(legend.position=c(.5,.9),
        legend.direction = 'horizontal',
        legend.title = element_text(hjust = 0.5 ,size = 11,color = 'black',face = 'bold'),
        legend.text = element_text(hjust = 0.5 ,size = 10,color = 'black'),
        axis.title = element_text(size = 12,face = 'bold',colour = 'black'), #darkred
        axis.text = element_blank(),#element_text(size = 12,colour = 'black'),
        panel.grid = element_blank()
  )

ggsave(filename = "1.annotation/tumor Cell colored by slingPseudotime.pdf",width = 4,height = 4)
# dev.off()

gene_key <- "PHOX2B"
plot_genes_in_pseudotime(cds[gene_key,], color_by = "celltype",cell_size = 1.5)+
  scale_color_manual(values = c("#ffc556","#76afda","#4c9568"))+
  theme_classic()+
  theme(legend.position="top",
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(hjust = 0.5 ,size = 10,face = 'bold',color = 'black'),
        axis.title = element_text(size = 12,face = 'bold',colour = 'black'), #darkred
        axis.text = element_text(size = 12,colour = 'black'),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = 'bold',colour = 'black')
  )


ggsave(filename = "1.annotation/tumor Cell colored by PHOX2B.pdf",width = 3.5,height = 3)


gene_key <- "PLK1"
plot_genes_in_pseudotime(cds[gene_key,], color_by = "celltype",cell_size = 1)+
  scale_color_manual(values = c("#ffc556","#76afda","#4c9568"))+
  theme_classic()+
  theme(legend.position="right",
        # legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(hjust = 0.5 ,size = 10,color = 'black'),
        axis.title = element_text(size = 13,face = 'bold',colour = 'black'), #darkred
        axis.text = element_text(size = 13,colour = 'black'),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = 'bold',colour = 'black')
  )


ggsave(filename = "1.annotation/tumor Cell colored by PLK1.pdf",width = 4,height = 3)

