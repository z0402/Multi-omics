library(ggsci)
library(tidyverse)
load("D:/nb/NBL2024.3.4/bulk/GSE49710/0.data/NB_exprs_list.Rda")
load("GSE49710/2.model construction/ML/cl_risk_list.Rda")
#E_MTAB_8248
E_MTAB_8248_group <- cl_list[[1]] %>% rownames_to_column("Sample") %>%
  mutate(group = ifelse(.$RiskScore > median(.$RiskScore),"high","low"))
rownames(E_MTAB_8248_group) <- E_MTAB_8248_group$Sample
expr <- ee[[1]]


write.table(expr,file = "GSE49710/0.data/E_MTAB_8248_expr.txt",row.names = T,sep = "\t",quote = F)


identical(colnames(expr),E_MTAB_8248_group$Sample)

library(limma)

group = E_MTAB_8248_group[,"group",drop=F] %>% as.data.frame()

design <- model.matrix(~0+factor(group$group))
colnames(design) = levels(factor(group$group))
rownames(design) = colnames(expr)

cts <- paste("high","low",sep="-")
contrast_matrix <- makeContrasts(contrasts=cts,levels = design)
fit <- lmFit(expr,design)
fit2 <- contrasts.fit(fit,contrast_matrix)
efit <- eBayes(fit2)

tempOutput <- topTable(efit, coef=paste0("high",'-',"low"), n=Inf) %>%na.omit()
tempOutput$gene <- rownames(tempOutput)

write.csv(tempOutput,file="GSE49710/3.functional analysis/Enrichment/deg.csv",quote = F)

deg_up <- subset(tempOutput,logFC > 1 & adj.P.Val < 0.05)
deg_down <- subset(tempOutput,logFC < -1 & adj.P.Val < 0.05)


library(clusterProfiler)
library(org.Hs.eg.db)


library(GseaVis)
library(ggrepel)
deg <- bitr(tempOutput$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb ="org.Hs.eg.db" ) %>%
  merge(.,tempOutput,by.x="SYMBOL",by.y="gene")
genelist <- deg$logFC
names(genelist) <- deg$ENTREZID
genelist <- sort(genelist,decreasing = T)

deg_up_gsea <- gseKEGG(genelist,organism = "hsa",keyType = "kegg")


df <- deg_up_gsea@result %>% arrange(desc(NES))
df$type <- ifelse(df$NES>0,"Activated","Repressed")
df_selected <- df[c(1:10,159,157,151,149,144,141,136,127:129),]


ggplot(df, aes(y = -log10(p.adjust), x = NES)) + 
  geom_point(aes(color = type),shape=16, alpha = 1, size = 2) + 
  geom_vline(xintercept = c(-1.5, 1.5), size = .75, lty = "dashed",
             color = "grey40") +
  geom_text_repel(data =df_selected ,force = 10,size=5,
                  # nudge_y = c(-0.9,0.9),
                  aes(y = -log10(p.adjust), x = NES,label = Description))+
  # geom_hline(yintercept = c(-1.5, 1.5), 
  #            size = 1, lty = "dashed", color = "grey") + 
  scale_colour_manual(name = "", values = c(`Activated` = "#ffc556",
                                            # `none sig` = "gray",
                                            `Repressed` = "#7fb961")) +
  xlab("NES") + ylab(paste("-log10(p.adjust)"))+
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.title = element_text(colour = "black", size = 16,face = "bold"), 
        legend.position = "top",
        legend.text = element_text(size = 14)) + 
  guides(color = guide_legend(override.aes = list(size = 5))) 

ggsave("GSE49710/3.functional analysis/Enrichment/KEGG gsea.pdf",width = 6,height =5.5)



#CIBERSORT----

cibersort <- read.delim("GSE49710/3.functional analysis/Enrichment/E_MTAB_8248_expr_cibersort_result.txt",check.names = F)

cibersort_group <- merge(cibersort,E_MTAB_8248_group[,c("Sample","group"),drop=F],by=1) %>%
  pivot_longer(cols = 2:23,names_to = "Type",values_to = "Fraction")
colnames(cibersort_group)

ggplot(cibersort_group,aes(Type,Fraction,fill=group))+
  geom_boxplot(outlier.shape = NA,width=.8,position = position_dodge(1.1))+
  scale_fill_manual(values = c("#ffc556","#7fb961"),name="")+
  stat_compare_means(label = "p.signif",hide.ns = T,size=6)+
  ylab("Cell fraction")+
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(angle = -90,hjust = 0,vjust = .5, size = 12),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = "black", size = 16,face = "bold"), 
        legend.position = c(.1,.85)) 

ggsave("GSE49710/3.functional analysis/Enrichment/cibersort.pdf",width =9,height =5.5)

# TCellSI-----
library(TCellSI)
library(ggpubr)
ResultScores_TCSS <- TCSS_Calculate(expr)
identical(colnames(ResultScores_TCSS),E_MTAB_8248_group$Sample)

TCSS_group <- merge(t(ResultScores_TCSS),E_MTAB_8248_group[,c("group"),drop=F],by=0) %>%
  pivot_longer(cols = 2:9,names_to = "Type",values_to = "Score")
colnames(TCSS_group)

ggplot(TCSS_group,aes(Type,Score,fill=group))+
  geom_boxplot(outlier.shape = NA,width=.4,position = position_dodge(.6),color="gray40")+
  scale_fill_manual(values = c("#ffc556","#7fb961"),name="")+
  stat_compare_means(label = "p.signif",hide.ns = T,size=6)+
  ylab("TCSS")+
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(angle = -90,hjust = 0,vjust = .5, size = 10),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = "black", size = 14,face = "bold"), 
        legend.position = "top") 

ggsave("GSE49710/3.functional analysis/Enrichment/TCSS.pdf",width =5,height =4)


#TISID----
library(GSVA)
immune_cell <- read.delim("GSE49710/0.data/TISID_immune_marker_genes.txt") 

library(reshape2)
geneSets = split(1:dim(immune_cell)[1], immune_cell$Cell.type) %>%
  lapply(function(x) {
    immune_cell[unlist(x), 1]
  })
head(geneSets, 2)

immune_cell_score <- gsva(as.matrix(expr),geneSets,
                          mx.diff=F,verbose=F,method="ssgsea",kcdf="Gaussian",
                          parallel.sz=4)

immune_cell_score <- t(immune_cell_score) %>% as.data.frame()
immune_cell_score_group <- merge(immune_cell_score,E_MTAB_8248_group[,c("Sample","group"),drop=F],by=0) %>%
  pivot_longer(cols = 2:29,names_to = "Type",values_to = "Score")


ggplot(immune_cell_score_group,aes(Type,Score,fill=group))+
  geom_boxplot(outlier.shape = NA,width=.7,position = position_dodge(1),color="gray40")+
  scale_fill_manual(values = c("#ffc556","#7fb961"),name="")+
  stat_compare_means(label = "p.signif",hide.ns = T,size=5)+
  ylab("Immune score")+
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", size = 10), 
        axis.text.x = element_text(angle = -90,hjust = 0,vjust = .5, size = 11),
        axis.title.x = element_blank(),
        axis.title = element_text(colour = "black", size = 14,face = "bold"), 
        legend.position = "top",
        legend.text = element_text(size = 11))

ggsave("GSE49710/3.functional analysis/Enrichment/TISID immune.pdf",width = 10.5,height =5.5)


#immune immunomodulator----
immune_molecular <- read.delim("GSE49710/0.data/immunomodulator.txt") %>%
  .[.$Gene %in% rownames(expr),]

immune_molecular_group <- expr[immune_molecular$Gene,] %>% t() %>% as.data.frame() %>%
  merge(.,E_MTAB_8248_group[,"group",drop=F],by=0) %>% .[,-1] %>%
  aggregate(.~group,.,mean)%>%
  column_to_rownames("group")%>% t() %>% as.data.frame()


identical(immune_molecular$Gene,rownames(immune_molecular_group))
# 创建列注释
annCol <- data.frame(group = c("high","low"),
                     row.names = c("high","low"))
library(ComplexHeatmap)
top_anno <- HeatmapAnnotation(df                   = annCol,
                              col                  = list(group= c(high="#ffc556",low="#7fb961")),
                              gp                   = gpar(col = "grey80"), 
                              simple_anno_size     = unit(3.5, "mm"), 
                              show_legend          = T, 
                              show_annotation_name = F,
                              border               = FALSE) 


annRow <- immune_molecular#Antigen presentation
annRow$Category <- factor(annRow$Category, 
                          levels = c("Co-inhibitor","Ligand","Receptor","Cell adhesion","Antigen presentation","Other")) # 由于行需要按照类分割，所以需要定义因子顺序，否则按照字母表
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "gray40","N/A" = "#888888","Stimulatory" = "#E59E02"))

left_anno <- HeatmapAnnotation(df                   = annRow[,"ICI",drop = F],
                               which                = "row", 
                               gp                   = gpar(col = "grey80"), 
                               col                  = annRowColors,
                               simple_anno_size     = unit(3.5, "mm"),
                               show_annotation_name = F,
                               border               = F)


library(circlize)
col_expr <- colorRamp2(c(-1,-.5,0,.5,1), 
                       c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")) 

col_expr <- colorRamp2(c(-1,-.5,0,.5,1), c("#440154","#404387","#345F8D","#44BF70","#7AD150","#FDE724"))


pdf("GSE49710/3.functional analysis/Enrichment/immune_molecular.pdf",height = 8,width = 4)
Heatmap(matrix             = as.matrix(t(scale(t(immune_molecular_group)))),#
        col                = viridis(10,begin = .5,alpha = .8),#col_expr,
        border             = NA, 
        rect_gp = gpar(col = "grey80"),
        cluster_rows       = F, 
        cluster_columns    = F, 
        show_row_names     = T, 
        row_names_side     = "left",
        row_names_gp       = gpar(fontsize = 8), 
        show_column_names  = F, 
        column_names_side  = "top", 
        row_split          = annRow$Category, 
        top_annotation     = top_anno, 
        left_annotation    = left_anno, 
        name               = "mRNA\nExpression",
        width              = ncol(immune_molecular_group) * unit(4, "mm"), 
        height             = nrow(immune_molecular_group) * unit(3.5, "mm")) 


dev.off()







#TIDE----

tide_group <- read.csv("GSE49710/3.functional analysis/Enrichment/E_MTAB_8248_tide_result.csv") %>%
  merge(.,E_MTAB_8248_group[,c("Sample","group"),drop=F],by=1)


ggplot(tide_group,aes(group,TIDE,fill=group))+
  geom_half_violin(side = "r",trim = T,color="white",width=.7,
                   position = position_nudge(x = .2))+
  scale_fill_manual(values = c("#ffc556","#7fb961"))+
  geom_half_boxplot(aes(color=group),side = "l",width=0.4,center = T,
                    outlier.shape = NA,errorbar.draw = F,
                    position = position_nudge(x = .1))+
  scale_color_manual(values = c("#ffc556","#7fb961"))+
  stat_compare_means(label="p.signif",hide.ns = F,size=4)+
  coord_flip()+
  ylab("TIDE")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y= element_text(size=12,color = "black"),
        axis.text.x = element_text(size=9,color = "black"),
        axis.title.x = element_text(face = "bold",size=12),
        legend.text = element_text(size = 10),
        legend.title = element_blank())


ggsave("GSE49710/3.functional analysis/Enrichment/TIDE.pdf", width =4, height = 2.5)










p.value = chisq.test(table(tide_group$Responder,tide_group$group))$p.value
p.value

p_data <- table(tide_group$Responder,tide_group$group) %>%
  as.data.frame() %>% set_names(.,c("Responder","group","n"))
ggplot(p_data,aes(group, n))+
  geom_col(aes(fill = Responder),position ="fill",width = .5)+
  ylab("Proportion")+
  theme_classic(base_rect_size = 1.5)+
  annotate(
    "text", label ="****",x = 1.5, y = 1.05, size =4
  )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face="bold",color="black"),
        panel.grid.major = element_line(linetype=2),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, size = 12,face="bold"),
        legend.key.size = unit(4,"mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "top")+
  scale_fill_manual(values = c("#ffc556","#7fb961"),label= c("NR","R"),name="")

ggsave("GSE49710/3.functional analysis/Enrichment/E_MTAB_8248_tide Responder.pdf",width = 3,height = 3.5)


ggboxplot(tide_group,x="group",y="TAM.M2")+ stat_compare_means()
#雷达图
library(fmsb)
radar <- tide_group[,11:16] %>% aggregate(.~group,.,median) %>%
  column_to_rownames("group")

radar <- rbind(data.frame(Dysfunction = c(0.5,-0.6),
                          Exclusion   = c(0.5,-0.6),   
                          MDSC  = c(0.5,-0.6),
                          CAF   = c(0.5,-0.6),
                          TAM.M2= c(0.5,-0.6)),
               radar
)

pdf("GSE49710/3.functional analysis/Enrichment/E_MTAB_8248_tide radarplot.pdf", width = 4, height = 4)

radarchart(radar, axistype =1,
           pcol = c("#ffc556","#7fb961"), 
           # pfcol = colors_in, 
           plwd = 2,
           cglcol = "grey40",cglty = 2,calcex = .6,
           axislabcol = "grey50", 
           caxislabels = format(round(seq(0, max(as.numeric(as.character(unlist(radar)))), length.out = 5), 1), nsmall = 0), 
           cglwd = 0.8, vlcex = .8)

legend(
  x = "bottom", legend = rownames(radar[-c(1,2),]), horiz =T,
  bty = "n", pch = 20 , col = c("#ffc556","#7fb961"),
  text.col = "black", cex = 1, pt.cex = 1
)
dev.off()


#-----

genesets <- readxl::read_xlsx("D:/bio/Immunotherapy Response Signature/immunotherapy-related pathways.xlsx",col_names = F) %>%
  column_to_rownames("...1") %>%
  t() %>% as.data.frame() %>% as.list()

genesets <- lapply(genesets, na.omit)

library(GSVA)
library(gghalves)
library(ggpubr)
immune_resistance_score <- gsva(as.matrix(expr),genesets,
                                mx.diff=F,verbose=F,method="ssgsea",kcdf="Gaussian",
                                parallel.sz=4)

immune_resistance_score <- t(immune_resistance_score) %>% as.data.frame()
immune_resistance_score_group <- merge(immune_resistance_score,E_MTAB_8248_group[,c("Sample","group"),drop=F],by=0) %>%
  pivot_longer(cols = 2:15,names_to = "Type",values_to = "Score")


ggplot(immune_resistance_score_group,aes(Type,Score,fill=group))+
  geom_half_violin(side = "r",trim = T,color="white",width=1.5,
                   position = position_nudge(x = .2))+
  scale_fill_manual(values = c("#ffc556","#7fb961"))+
  geom_half_point(aes(color=group),side = "l",size=.2,width=0.5,
                  position = position_nudge(x = .2))+
  scale_color_manual(values = c("#ffc556","#7fb961"))+
  stat_compare_means(label="p.signif",hide.ns = F,size=4)+
  coord_flip()+
  ylab("Pathway score")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.x= element_text(size=10,color = "black"),
        axis.text.y = element_text(size=9,color = "black"),
        axis.title.x = element_text(face = "bold",size=12),
        legend.title = element_blank())


ggsave("GSE49710/3.functional analysis/Enrichment/immune resistance.pdf", width =6.5, height = 5.5)



