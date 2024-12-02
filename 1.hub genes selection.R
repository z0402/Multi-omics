
rm(list = ls())
library(tidyverse)
library(broom)
library(glmnet)
library(survival)
library(rms)
load("bulk/GSE49710/0.data/NB_exprs_list.Rda")
load("bulk/GSE49710/0.data/NB_cl_list.Rda")

# 基因的单因素/多因素回归 ----
GSE49710_cl <- cc[[4]]
GSE49710_expr <- ee[[4]]
tpm <- merge(GSE49710_cl,t(GSE49710_expr),by=0) %>% column_to_rownames(var = "Row.names")
colnames(tpm) <- str_replace_all(colnames(tpm),"-","_")



UniCox <- function(x){
  FML <- as.formula(paste0("Surv(OS_Time, OS_Event)~",x))
  Cox <- coxph(FML,data = tpm)
  Sum <- summary(Cox)
  HR <- round(Sum$coefficients[,2],2)
  PValue <- round(Sum$coefficients[,5],3)
  CI5 <- round(Sum$conf.int[,3],2)
  CI95 <- round(Sum$conf.int[,4],2)
  UniCox <- data.frame("characteristics" =x,
                       "Hazard Ratio" = HR,
                       "CI5"=CI5,
                       "CI95"=CI95,
                       "P Value" = PValue)
  
  return(UniCox)
  
} 

VarNames <- colnames(tpm)[-c(1:10)]
UniVar <- lapply(VarNames,UniCox)
UniVar <- do.call(rbind,UniVar)
table(UniVar$P.Value < 0.05 & UniVar$Hazard.Ratio > 1)

save(UniVar,file = "GSE49710/2.model construction/ML//Univar.Rda")


selected_genes1 <- UniVar[UniVar$P.Value < 0.05 & UniVar$Hazard.Ratio > 1,"characteristics"]
selected_genes2 <- UniVar[UniVar$P.Value < 0.05 & UniVar$Hazard.Ratio < 1,"characteristics"]
selected_genes <- c(selected_genes1,selected_genes2)
save(selected_genes1,selected_genes2,selected_genes,file = "GSE49710/2.model construction/ML//Univar Cox.Rda")



#hub genes KEGG----
load("GSE49710/1.WGCNA/WGCNA_hub_genes.Rda")#hub genes
load("D:/nb/NBL2024.3.4/bulk/GSE49710/6.depmap/depmap_genes.rda")#Depmap
load("GSE49710/2.model construction/ML//Univar Cox.Rda")#Uni selected genes

pre_var <- intersect(hub_genes,depmap_genes) %>% intersect(selected_genes,.)

library(clusterProfiler)
library(org.Hs.eg.db)
entrez_id <- bitr(pre_var,fromType = "SYMBOL",toType = "ENTREZID",OrgDb ="org.Hs.eg.db" )
kegg <- enrichKEGG(entrez_id$ENTREZID,organism = "hsa")

result <- kegg@result %>% filter(p.adjust < 0.05) %>% .[,c(4,8,10)]


probel2symbol <- result[,c(1,3),drop=F] 
tmp = unlist(lapply(1:nrow(probel2symbol),function(i){
  gene = trimws(unlist(strsplit(probel2symbol[,2][i],"/")))
  names(gene) = rep(probel2symbol[,1][i],length(gene))
  return(gene)
}))
probel2symbol= data.frame(ID= names(tmp),ENTREZID =tmp)
probel2symbol <- left_join(probel2symbol,entrez_id,by = "ENTREZID")
probel2symbol <- split(probel2symbol[,3],probel2symbol[,1])

dat = data.frame("PCNA, POLD1, POLE2, RFC4",
                 "AURKB, CCNA2, CDC6, CDT1, CHEK1, MAD2L1, MCM2, NDC80, PCNA, PKMYT1, PLK1, WEE1"  ,
                 "MCM2, PCNA, POLA2, POLD1, POLE2, RFC4" ,
                 "RRM1, RRM2",
                 "BIRC5, CCNA2, PCNA" ,
                 "POLD1, RAD51",
                 "PCNA, POLD1, RFC4" ,
                 "PCNA, POLD1, POLE2, RFC4" ,
                 "MAD2L1, PKMYT1, PLK1" ,
                 "CCNA2, MAD2L1, PKMYT1, PLK1"  ,
                 "RRM1, RRM2"
                 
) %>% set_names(.,names(probel2symbol)) %>% t() %>%
  as.data.frame() %>% mutate(Pathway=rownames(.)) %>%
  .[result1$Description,]

result$geneID <- dat$V1

result$Description <- factor(result$Description, levels = rev(result$Description))


p <- ggplot() +
  geom_bar(data = result,
           aes(x = -log10(p.adjust), y = Description,fill=Description),
           width = 0.4, #柱子宽度调整
           stat = 'identity') +
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  geom_text(data = result,
            aes(x = 0.1, 
                y = Description,
                label = Description),
            size = 4,
            hjust = 0) +#左对齐
  geom_text(data = result,
            aes(x = 0.1, y = Description, label = geneID),
            size = 3,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3) +
  labs(x = '-Log10 (FDR)',
       y = NULL, 
       title = 'Enriched pathways of 40 hub genes')+
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 12, face = 'bold'),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, face = 'bold.italic'),
    axis.ticks.y = element_blank())+
  scale_fill_manual(values = mycol22[-1]) +
  scale_color_manual(values = mycol22[-1])
p
ggsave("GSE49710/2.model construction/ML/KEGG.pdf",width = 4,height = 5)



library(ggvenn)
genes <- list('key drivers' = hub_genes,
              'prognostic genes' = selected_genes,
              'essential genes'=depmap_genes
)
ggvenn(genes,show_elements = F,show_percentage = F,
       fill_alpha = .8,
       fill_color = c("#33B3A6","#BCD0DE","#F8C0A6"),
       stroke_size = .4)

ggsave("GSE49710/2.model construction/ML/hub genes.pdf",width = 4,height =4)



##
load("GSE49710/2.model construction/ML//Univar.Rda")

pro.df <- UniVar[UniVar$characteristics %in% pre_var,] %>% arrange(.,desc(Hazard.Ratio))
pro.df$characteristics <- as.factor(pro.df$characteristics)
pro.df$label <- paste0(pro.df$Hazard.Ratio,"(",pro.df$CI5,"-",pro.df$CI95,")")



ggplot(pro.df)+
  geom_errorbar(aes(xmin = CI5, xmax = CI95,y=characteristics), width = 0)+
  geom_point(aes(Hazard.Ratio,characteristics,shape=characteristics),show.legend = F,color="#f06152")+
  scale_shape_manual(values = rep(15,40) )+
  scale_x_continuous(limits = c(-3.5,7.5),breaks = c(0,1,2,4,6,8),labels = c(0,1,2,4,6,8))+
  geom_vline(xintercept = 1, linetype = 2,color="black")+
  geom_text(aes(-.5,characteristics,label=label), size = 3)+
  geom_text(aes(-3,characteristics), size = 3,label="p < 0.0001")+
  theme_classic()+
  xlab("HR(95%)")+
  labs(title = "Univariate Cox analysis")+
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 10,color="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.ticks.y = element_blank())

ggsave("GSE49710/2.model construction/ML/hub genes HR.pdf",width = 5,height =9.5)


###CRISPR-----
library(data.table)
library(tidyverse)
df_effect <- fread("CRISPR/CRISPRGeneEffect.csv")
df_effect <- column_to_rownames(df_effect,"V1")
colnames(df_effect) <- str_split(colnames(df_effect)," ",simplify = T)[,1]

model <- fread("CRISPR/Model.csv")
model <- model[,c(1,3,5,10,17,18)] %>% .[.$DepmapModelType=="NBL",]

df_effect1 <- df_effect[which(rownames(df_effect) %in% model$ModelID),]


data <- df_effect1
data[data >= -1] <- 0
data[data < -1] <- 1
data1 <- data[,colSums(data,na.rm = T) >= nrow(data)*.75]
depmap_genes <- colnames(data1)



df_effect2 <- df_effect1[,pre_var] %>% rownames_to_column("Cell") %>%
  pivot_longer(cols = -1,names_to = "gene")

ggplot(df_effect2)+
  geom_boxplot(aes(gene,value),width=.7,outlier.size = .1)+
  geom_jitter(aes(gene,value),width = .3,size=1.2,color="gray60")+
  geom_hline(yintercept = -1, linetype = 2,color="#f06152")+
  labs(y="CERES")+
  theme_classic()+
  theme(
    legend.position = 'none',
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11,color="black",angle = -90,hjust = 0,vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = 'bold.italic'))


ggsave("GSE49710/2.model construction/ML/hub genes CERES.pdf",width = 10,height =4.5)




#Mal1 score----
library(tidyverse)
library(broom)
library(glmnet)
library(survival)
library(rms)
load("D:/nb/NBL2024.3.4/bulk/GSE49710/0.data/NB_exprs_list.Rda")
load("D:/nb/NBL2024.3.4/bulk/GSE49710/0.data/NB_cl_list.Rda")
Mal_sig <- read.delim("D:/nb/NBL2024.3.4/GSE192906/Mal_genes.txt")
Mal1_sig <- list(Mal1_sig=Mal_sig[Mal_sig$cluster=="1","gene"])


i=4 
names(ee)[i]
expr <- ee[[i]]
cl <- cc[[i]]
library(GSVA)
Mal_sig_score <- gsva(as.matrix(expr),Mal1_sig,method="ssgsea",
                      kcdf="Gaussian" ,
                      verbose=T)

Mal_sig_score <- as.data.frame(Mal_sig_score)
Mal_sig_score <- merge(cl,t(Mal_sig_score),by=0)


library(IOBR)
cell.cycle <- list(cell.cycle=IOBR::kegg$KEGG_CELL_CYCLE)
cell.cycle_score <- gsva(as.matrix(expr),cell.cycle,method="ssgsea",
                         kcdf="Gaussian" ,
                         verbose=T)

cell.cycle_score <- as.data.frame(cell.cycle_score)
cell.cycle_score <- merge(Mal_sig_score,t(cell.cycle_score),by.x=1,by.y=0)

ggscatter(cell.cycle_score,x="cell.cycle",y="Mal1_sig",
          size = 1.5,color ="#76afda",
          add = "reg.line",  
          add.params = list(color = "#5066a1", fill = "gray"), 
          conf.int = TRUE)+
  stat_cor(method = "spearman",size=5)+
  geom_rug(color ="#76afda")+
  labs(y='Mal1 score',x='Cell cycel score')+
  theme(legend.position = c(0.12,0.3),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'black',face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'black',face='bold'))


ggsave(filename = paste0("D:/nb/NBL2024.3.4/GSE192906/1.annotation/tumor Cell Cell cycle ",names(ee)[i],".pdf"),width = 4,height =4)







library(survival)
library(survminer)
library(ggquickeda)

Mal_sig_score$group <- factor(ifelse(Mal_sig_score[,"Mal1_sig"] > median(Mal_sig_score[,"Mal1_sig"]),"high","low"), levels = c("low","high"))


fit <- survdiff(Surv(EFS_Time, EFS_Event) ~ group, data=Mal_sig_score, na.action=na.exclude)
fit$pvalue
ggplot(Mal_sig_score, aes(time =EFS_Time, status = EFS_Event,color = group)) + 
  geom_km(size=1) + geom_kmticks()+
  theme_classic()+ 
  scale_color_manual(values = c("#ffc556","#76afda" ))+
  labs(y='Event free survival',x='Time (Days)')+#Event free   Overall
  theme(legend.position = c(0.12,0.3),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'black',face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'black',face='bold'))+
  annotate('text',x = 0.1,y=0.45,size=4,label='Log-rank\np < 0.0001',hjust=0,fontface = 'italic')


ggsave(filename = paste0("D:/nb/NBL2024.3.4/GSE192906/1.annotation/tumor Cell KM EFS ",names(ee)[i],".pdf"),width = 4,height =4)



ggplot(Mal_sig_score, aes(INSS, Mal1_sig,fill = INSS)) + 
  geom_boxplot(width=.5,outlier.shape = NA)+
  theme_classic()+ 
  scale_fill_manual(values = c("#b0d45d","#76afda","#ffc556","#eb998b","#f06152" ))+
  labs(y='Mal1 score')+
  stat_compare_means(label.x.npc = .1,size=4.5)+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'black',face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,colour = 'black',face='bold'))

ggsave(filename = paste0("D:/nb/NBL2024.3.4/GSE192906/1.annotation/tumor Cell KM INSS ",names(ee)[i],".pdf"),width = 4,height =3.5)




library(survminer)
library(survival)
library(ggquickeda)


outTab.cox <- NULL

for(i in 1:5) {
  exprSurvSub <- ee[[i]]
  clSurvSub <- cc[[i]]
  Mal_sig_score <- gsva(as.matrix(exprSurvSub),Mal1_sig,method="ssgsea",
                        kcdf="Gaussian" ,
                        verbose=T)
  Mal_sig_score <- as.data.frame(Mal_sig_score)
  Mal_sig_score <- merge(clSurvSub,t(Mal_sig_score),by=0)
  
  
  Mal_sig_score$group <- factor(ifelse(Mal_sig_score$Mal1_sig > median(Mal_sig_score$Mal1_sig),"high","low"), levels = c("low","high"))
  fit <- survdiff(Surv(OS_Time, OS_Event) ~ group, data=Mal_sig_score, na.action=na.exclude)
  
  ## OS
  coxres <- summary(coxph(as.formula(paste("Surv(OS_Time, OS_Event) ~", "group")), data = Mal_sig_score))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(Dataset = i, 
                                            event = "OS", 
                                            beta = coxres$coefficients[1,1], 
                                            hr = coxres$coefficients[1,2],
                                            lower = coxres$conf.int[1,3], 
                                            upper = coxres$conf.int[1,4], 
                                            p = coxres$coefficients[1,5],
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
  
  
}

outTab.cox$Dataset <- str_remove_all(names(cc),"_cl")


data <- outTab.cox[c(4:2,1,5),]
data$Dataset <- as.factor(data$Dataset )
data$Dataset <- toupper(data$Dataset )
data$ll <- ifelse(data$p<0.0001,'****',ifelse(data$p<0.001,'***',ifelse(data$p<0.01,'**',ifelse(data$p<0.05,'*',''))))
data$label <- paste0(sprintf("%0.2f",data$lower),"-",sprintf("%0.2f",data$upper))
data$label <- paste0(sprintf("%0.2f",data$hr),"(",data$label,")")

data$Dataset[4] <- "E-MTAB-8248"


data$upper[1] <- 29.9
data$upper[2] <- 29.9

data <- data[1:4,]

ggplot(data,aes(hr,fct_infreq(Dataset,hr)))+
  geom_vline(xintercept = 1,linetype=2,color='black')+
  geom_errorbar(aes(xmin=lower,xmax=upper),width=0,size=.5,
                position = position_dodge(width = 0.6),color = 'black')+
  geom_point(shape=15,size=4,position = position_dodge(width = 0.6),
             color="#ffc556")+
  geom_text(aes(x = min(lower)-22.5, y = Dataset, label=ll),vjust=.7,hjust = 0.5,size=4,
            position = position_dodge(width = 0.6),color = 'black')+ 
  geom_text(aes(x = min(lower)-20, y = Dataset, label=label),vjust=.5,hjust = 0,size=4,
            position = position_dodge(width = 0.6),color = 'black')+ 
  scale_x_continuous(expand = c(0,0),limits = c(-22,30),breaks = c(0,5,10,15,20,25,30))+ ## 调整
  labs(x='HR 95%CI',title =  'Univariable Cox Analysis')+
  theme_classic(base_rect_size = 1)+
  theme(                                     
    axis.text.x = element_text(size = 10,color = "black"),
    axis.text.y = element_text(size = 12,color = "black"),
    axis.title.x = element_text(size = 10,face = 'bold',hjust = 0.1,vjust = 7.5), ## 调整
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    legend.position = 'none',
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 13,face='bold'))+
  guides(color=guide_legend(nrow = 1))

ggsave(filename = "D:/nb/NBL2024.3.4/GSE192906/1.annotation/tumor Cell KM OS mutiple.pdf",width = 5.5,height =3.5)


#clinical characteristics----
dat = Mal_sig_score
dat$Age <- ifelse(dat$Age < 547,"< 18 months",">= 18 months")
colnames(dat)[8] <- "MYCN"
dat$MYCN <- ifelse(is.na(dat$MYCN),NA,ifelse(dat$MYCN == "0","Non-Amp","Amp"))
dat$Risk <- ifelse(dat$Risk == "0","low-risk","high-risk")
dat$Grade <- ifelse(is.na(dat$Grade),NA,ifelse(dat$Grade == "0","F","UF"))

dat$Age <- factor(dat$Age,levels = c("< 18 months",">= 18 months"))
dat$MYCN <- factor(dat$MYCN,levels = c("Non-Amp","Amp"))
dat$Risk <- factor(dat$Risk,levels = c("low-risk","high-risk"))
dat$INSS <- factor(dat$INSS,levels = c("1","2","3","4","4S"))
dat$Grade <- factor(dat$Grade,levels = c("F","UF"))


gname <- "group"
vname <- c("Age","MYCN","Risk","INSS","Grade")
pie.high <- pie.low <- list()
chisq.p <- c()
for (i in vname) {
  tmp <- table(dat[,gname], dat[,i])
  p <- format(chisq.test(tmp)$p.value,digits = 2)
  names(p) <- i
  chisq.p <- c(chisq.p, p)
  
  pie.dat <- 
    tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
  
  pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "high"),]
  pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "low"),]
}


green  <- "#7fb961"
blue   <- "#76afda"
orange <- "#ffc556"
cherry <- "#9e6c69"
red <- "#f06152"
mycol22 <- c("#b0d45d","#7fb961","#4c9568","#356d67","#42465c","#5066a1",
             "#76afda","#abddff","#dcf2ff","#fddbc8","#ffe788","#ffc556",
             "#e8743c","#f06152","#b20000","#eb998b","#a14462","#cca69c",
             "#9e6c69","#7d4444","#562e3c","#35212e")

# 创建颜色
Age.col <- alpha(green, c(0.4,0.8 ))
MYCN.col <- alpha(blue,  c(0.4,0.8))
Risk.col <- alpha(orange,  c(0.4,0.8))
INSS.col <- alpha(red, c(.2,.4,0.6,1 ,0.8))
Grade.col <- alpha(cherry, c(0.4,0.8))


pdf("../GSE192906/1.annotation/pieTable.pdf",width = 5, height = 7)
showLayout <- F 

{
  layout(matrix(c(   1, 7, 7, 7,  13,13,13,  19,19,  25,
                     1, 7, 7, 7,  13,13,13,  19,19,  25,
                     2, 8, 8, 8,  14,14,14,  20,20,  25,
                     2, 8, 8, 8,  14,14,14,  20,20,  25,
                     3, 9, 9, 9,  15,15,15,  21,21,  25,
                     3, 9, 9, 9,  15,15,15,  21,21,  25,
                     4, 10,10,10, 16,16,16,  22,22,  25,
                     4, 10,10,10, 16,16,16,  22,22,  25,
                     5, 11,11,11, 17,17,17,  23,23,  25,
                     5, 11,11,11, 17,17,17,  23,23,  25,
                     6, 12,12,12, 18,18,18,  24,24,  25,
                     6, 12,12,12, 18,18,18,  24,24,  25
  ),
  byrow =  T,nrow =12)
  )
  
  if(showLayout) {
    layout.show(n = 25) 
  }

  
  par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2)
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n") 
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") 
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "NB",cex = 2, col = "black") 
  
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n") 
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") 
  text((par("usr")[1]+par("usr")[2])/2,
       (par("usr")[3]+par("usr")[4])/2,
       "Age",cex = 2, col = "black")
  
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n") 
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "MYCN",cex = 2, col = "black") 
  
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "Risk",cex = 2, col = "black") 
  
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "INSS",cex = 2, col = "black") 
  
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "Grade",cex = 2, col = "black") 
  

  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "high\n(n = 249)",cex = 2, col = "black") 
  
  # High group
  pie(pie.high$Age$Pct, 
      col = Age.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$MYCN$Pct, 
      col = MYCN.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$Risk$Pct, 
      col = Risk.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$INSS$Pct, 
      col = INSS.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.high$Grade$Pct, 
      col = Grade.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  abline(h = par("usr")[3], col = "black")

  
  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "low\n(n = 249)",cex = 2, col = "black") 
  
  # Low group
  pie(pie.low$Age$Pct, 
      col = Age.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$MYCN$Pct, 
      col = MYCN.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$Risk$Pct, 
      col = Risk.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$INSS$Pct, 
      col = INSS.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  
  pie(pie.low$Grade$Pct, 
      col = Grade.col, 
      border = "white",  
      radius = 1, 
      labels = NA,
      init.angle = 90)
  symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
  abline(h = par("usr")[3], col = "black")
  

  plot(1,1,
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       "p-value",cex = 2, col = "black") 
  
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",chisq.p["Age"]),cex = 1.5, col = "black") 
  abline(v = par("usr")[4], col = "black") 
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",chisq.p["MYCN"]),cex = 1.5, col = "black") 
  abline(v = par("usr")[4], col = "black")
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n",
       ylab = "",yaxt = "n")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",chisq.p["Risk"]),cex = 1.5, col = "black") 
  abline(v = par("usr")[4], col = "black")
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",chisq.p["INSS"]),cex = 1.5, col = "black") 
  abline(v = par("usr")[4], col = "black")
  
  plot(1,1,col = "white",
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  text((par("usr")[1]+par("usr")[2])/2, 
       (par("usr")[3]+par("usr")[4])/2,
       paste0("p = ",chisq.p["Grade"]),cex = 1.5, col = "black") 
  abline(v = par("usr")[4], col = "black") 
  abline(h = par("usr")[3], col = "black")

  
  plot(0,0,col = "white",
       xlab = "",xaxt = "n", 
       ylab = "",yaxt = "n")
  legend("left",
         legend = c("< 18 months",">= 18 months",
                    "Non-Amp","Amp",
                    "low-risk","high-risk",
                    "1","2","3","4","4S",
                    "F","UF"),
         fill = c(Age.col,
                  MYCN.col,
                  Risk.col,
                  INSS.col,
                  Grade.col),
         border = NA, 
         bty = "n", 
         cex = 1.2,
         #box.lwd = 3,
         x.intersp = 0.05,
         y.intersp = 1,
         text.width = 0.075, 
         horiz = F) 
}

invisible(dev.off())























#PLK1----
load("D:/nb/NBL2024.3.4/bulk/GSE49710/0.data/NB_exprs_list.Rda")
load("GSE49710/2.model construction/ML/cl_risk_list.Rda")
load("GSE49710/2.model construction/pre_var.Rda")
library(survminer)
library(survival)
library(ggquickeda)
# minprop=0.1

OS_list <- lapply(cl_list[1],function(x){x[,c("INSS",'OS_Time','OS_Event',"RiskScore")]}) %>%
  do.call(rbind,.)  %>% 
  mutate(Dataset=str_split(str_split(rownames(.),"_cl",simplify = T)[,1],"_dataset",simplify = T)[,1],
         Sample=str_split(rownames(.),"_cl.",simplify = T)[,2]) %>%
  merge(.,t(ee[[1]])[,"PLK1",drop=F],by.x="Sample",by.y=0)



EFS_list <- lapply(cl_list[1],function(x){x[,c('EFS_Time','EFS_Event',"RiskScore")]}) %>%
  do.call(rbind,.)  %>% 
  mutate(Dataset=str_split(str_split(rownames(.),"_cl",simplify = T)[,1],"_dataset",simplify = T)[,1],
         Sample=str_split(rownames(.),"_cl.",simplify = T)[,2]) %>%
  merge(.,t(ee[[1]])[,pre_var,drop=F],by.x="Sample",by.y=0)

library(ggquickeda)


exprSurvSub <- OS_list


exprSurvSub$group <- factor(ifelse(exprSurvSub$PLK1 > median(exprSurvSub$PLK1),"high","low"),levels = c("low","high"))

fit <- survfit(Surv(OS_Time,OS_Event)~group, data=exprSurvSub, na.action=na.exclude)
message(paste0('>>> KM P = ',surv_pvalue(fit)$pval))
p <- surv_pvalue(fit)$pval
p <- formatC(p,format="e",digits=3)
ggplot(exprSurvSub, aes(time =OS_Time, status = OS_Event,color = group)) + 
  geom_km(size=1) + geom_kmticks()+
  theme_classic(base_rect_size = 1)+ 
  scale_color_manual(values = c("#ffc556","#76afda" ))+
  labs(y='Overall survival',x='Time (Days)')+
  ggtitle("E-MTAB-8248")+
  theme(legend.position = c(0.13,0.25),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 15,colour = 'black',face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,colour = 'black'))+
  annotate('text',x = 0.1,y=0.45,size=4,label=paste0('Log-rank\np = ',p),hjust=0,fontface = 'italic')
ggsave("GSE49710/clinical/E-MTAB-8248 PLK1 OS.pdf",width = 4,height = 4)




ggplot(OS_list, aes(INSS, PLK1,fill = INSS)) + 
  geom_boxplot(width=.5,outlier.shape = NA)+
  theme_classic()+ 
  scale_fill_manual(values = c("#b0d45d","#76afda","#ffc556","#eb998b","#f06152" ))+
  labs(y='PLK1 expression')+
  stat_compare_means(label.x.npc = .1,size=4)+
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 12,colour = 'black',face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,colour = 'black',face='bold'))

ggsave(filename ="GSE49710/clinical/E-MTAB-8248 PLK1 INSS.pdf",width =3.5,height =3.2)






#PLK1 BP
E_MTAB_8248_expr <- t( ee[[1]])
cor_data <- cor(E_MTAB_8248_expr,E_MTAB_8248_expr[,"PLK1"])
cor_data_genes <- rownames(cor_data)[cor_data > 0.8]

write.table(cor_data_genes,file = "PLK1 data/PLK1 NEW/PLK1 cor_data_genes.txt",sep = "\t",quote = F,row.names = F)

cor_data <- as.data.frame(cor_data)
entrez_id <- bitr(rownames(cor_data),fromType = "SYMBOL",toType = "ENTREZID",OrgDb ="org.Hs.eg.db" ) %>%
  merge(.,cor_data,by.x=1,by.y=0)

genelist <- entrez_id$V1
names(genelist) <- entrez_id$ENTREZID
genelist <- sort(genelist,decreasing = T)

gsea <- gseGO(genelist,ont = "BP","org.Hs.eg.db")
result <- gsea@result

result$Description <- factor(result$Description, levels = rev(result$Description))
ggplot(data = result[1:10,]) +
  geom_bar(aes(x =NES, y = Description,fill=Description),
           width = 0.7, #柱子宽度调整
           stat = 'identity') +
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  geom_text(
    aes(x = 0.1, 
        y = Description,
        label = Description),
    size = 4,fontface = 'italic',
    hjust = 0) +
  labs(x = 'NES',
       y = NULL,
       title = 'Biological Process')+
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10,color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12,color = "black", face = 'bold.italic'),
    axis.ticks.y = element_blank())+
  scale_fill_manual(values = mycol22[-c(5,6)]) +
  scale_color_manual(values = mycol22[-c(5,6)])

ggsave("GSE49710/clinical/E-MTAB-8248 PLK1 BP.pdf",width = 3,height = 4)


