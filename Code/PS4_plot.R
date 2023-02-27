##R script for plot in PS4 paper

setwd("F:/Manuscript_PS4/")

#load packages
library(reshape2)
library(ggsci)
library(ggplot2)
library(scales)
library(dplyr)
library(paletteer)
library(ggstatsplot)
library(gridExtra)
library(patchwork)
library(bootLR)
library(ggthemes)
library(tidyr)
library(tidyverse)
library(viridis)
library(patchwork)
library(networkD3)
library(webshot)
library(d3Network)
library(ComplexHeatmap)
library(cowplot)

*****************************************
* Figure 2, Figure S2 and Figure S3
*****************************************

########################
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )
  
pie_input<-data.frame(table(AR$Diagnosis))
pie_input$Diagnosis<-factor(pie_input$Var1,levels = c('P', 'LP', 'VUS','LB','B'))

AR1<-ggplot(pie_input, aes(x="", y=Freq, fill=Diagnosis))+geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+ blank_theme +theme(axis.text.x=element_blank())+guides(fill=FALSE)

P<-as.data.frame(table(AR[AR$Diagnosis=="P"|AR$Diagnosis=="LP",]$Gene))
P<-data.frame(AR %>%
    group_by(Diagnosis,Gene) %>%
    summarize(n=n()) %>%
    arrange(desc(n),Gene))
P<-P[P$Diagnosis=="P"|P$Diagnosis=="LP",]
data<-P
gene<-data.frame(data %>%
    group_by(Gene) %>%
    summarize(sum(n)))
gene <- arrange(gene,desc(sum.n.))
data$Gene<-factor(data$Gene,levels=as.character(unique(gene$Gene)))
data$Diagnosis<-factor(data$Diagnosis,levels=c("P","LP"))
AR2<-ggplot(data=data[data$Gene%in%gene$Gene[1:20],], mapping=aes(x=Gene,y=n,fill=Diagnosis))+
    geom_bar(stat="identity",colour = 'black', width = 1, lwd=0.1)+ theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8))+labs(fill="", x = "", y = "Number of P/LP variants")+scale_fill_manual(values=c("#780000","#c1121f"))+guides(fill=FALSE)

AR$Diagnosis<-factor(AR$Diagnosis,levels = c('P', 'LP', 'VUS','LB','B'))
p1<-ggplotGrob(ggplot(data=AR,aes(x=F_A,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.00005, 
                 color="black", alpha=0.9)+xlim(NA,0.001)+ylim(NA,6000)+xlab("Minor allele frequency in cases")+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="",y="")+guides(fill=FALSE))

AR3<-ggplot(data=AR,aes(x=F_A,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.0003, 
                 color="black", alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="Allele frequency in cases",y="Count")+annotation_custom(grob = p1,xmin = 0.0001, xmax = 0.01, ymin = 0, ymax = 7000)+xlim(NA,0.01)+guides(fill=FALSE)

p1<-ggplotGrob(ggplot(data=AR,aes(x=F_U,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.00005, 
                 color="black", alpha=0.9)+xlim(NA,0.001)+ylim(NA,6000)+xlab("MAF")+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="",y="")+guides(fill=FALSE))

AR4<-ggplot(data=AR,aes(x=F_U,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.0003, 
                 color="black", alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="Allele frequency in controls",y="Count")+annotation_custom(grob = p1,xmin = 0.0001, xmax = 0.01, ymin = 0, ymax = 7000)+xlim(NA,0.01)+guides(fill=FALSE)

cover<-data.frame(AR %>%
    group_by(Diagnosis,Diagnosis_3,VEP_region) %>%
    summarize(n=n()) %>%
    arrange(desc(n),VEP_region))
count<-data.frame(AR %>%
    group_by(Diagnosis_3) %>%
    summarize(n=n()))
cover<-full_join(cover,count,by=c("Diagnosis_3"))
cover<-cover %>% mutate(percent=n.x/n.y)
cover<-add(cover,1,1)
cover2<-data.frame(AR %>%
  group_by(Diagnosis_3,VEP_region) %>% 
  summarize(n=n()) %>%
  arrange(desc(n),VEP_region))
count<-data.frame(AR %>%
  group_by(Diagnosis_3) %>% 
  summarize(n=n()))
cover2<-full_join(cover2,count,by=c("Diagnosis_3"))
cover2<-cover2 %>% mutate(percent=n.x/n.y)
cover$VEP_region<-factor(cover$VEP_region,levels=as.character(unique(c(cover2[cover2$Diagnosis_3=="P/LP",]$VEP_region,cover2[cover2$Diagnosis_3=="VUS",]$VEP_region,cover2[cover2$Diagnosis_3=="B/LB",]$VEP_region))))
cover$Diagnosis_3<-factor(cover$Diagnosis_3,levels=c("P/LP","VUS","B/LB"))
cover$Diagnosis<-factor(cover$Diagnosis,levels = c('P', 'LP', 'VUS','LB','B'))
AR5<-ggplot(cover,aes(x=VEP_region,y=percent,fill=Diagnosis))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+guides(fill=guide_legend(title=NULL))+theme(axis.title = element_blank(),legend.position='none',axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8))+facet_grid(Diagnosis_3~.)+ylim(0,1)+guides(fill=FALSE)

pdf("Variants_basic_info.pdf",height=9,width=9)
((AR1 | AR2) / (AR5 | AR2) / (AR3|AR4))+plot_layout(ncol = 1,widths = c(1, 2),guides = 'collect')+ plot_annotation(tag_levels = 'A')
dev.off()

pdf("Fiugre2_PartB.pdf",height=6,width=5)
cover<-data.frame(AR %>%
  group_by(Diagnosis_3,VEP_region) %>% 
  summarize(n=n()) %>%
  arrange(desc(n),VEP_region))
count<-data.frame(AR %>%
  group_by(Diagnosis_3) %>% 
  summarize(n=n()))
cover<-full_join(cover,count,by=c("Diagnosis_3"))
cover<-cover %>% mutate(percent=n.x/n.y)
cover$Diagnosis_3<-factor(cover$Diagnosis_3,levels=c("P/LP","VUS","B/LB"))
cover$VEP_region<-factor(cover$VEP_region,levels=as.character(unique(c(cover[cover$Diagnosis_3=="P/LP",]$VEP_region,cover[cover$Diagnosis_3=="VUS",]$VEP_region,cover[cover$Diagnosis_3=="B/LB",]$VEP_region))))
ggplot(cover,aes(x=VEP_region,y=percent,fill=Diagnosis_3))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7"))+guides(fill=guide_legend(title=NULL))+theme(axis.title = element_blank(),legend.position='none')+ geom_text(aes(label = n.x),size=2,vjust = 0, nudge_y = 0.01)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8),plot.title = element_text(hjust = 0.5))+facet_grid(Diagnosis_3~.)+ylim(0,1)
dev.off()

pdf("Fiugre2_PartD.pdf",height=3,width=6)
p1<-ggplotGrob(ggplot(data=AR,aes(x=F_A,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.00005, 
                 color="black", alpha=0.9)+xlim(NA,0.001)+ylim(NA,6000)+xlab("Minor allele frequency in cases")+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="",y="")+guides(fill=FALSE))

ggplot(data=AR,aes(x=F_A,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.005, 
                 color="black", alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="Allele frequency in cases",y="Count")+annotation_custom(grob = p1,xmin = 0.05, xmax = 0.45, ymin = 0, ymax = 9000)+xlim(NA,0.5)+guides(fill=FALSE)
dev.off()

pdf("Fiugre2_PartE.pdf",height=3,width=6)
p1<-ggplotGrob(ggplot(data=AR,aes(x=F_U,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.00005, 
                 color="black", alpha=0.9)+xlim(NA,0.001)+ylim(NA,6000)+xlab("MAF")+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="",y="")+guides(fill=FALSE))

ggplot(data=AR,aes(x=F_U,fill=Diagnosis)) + 
  geom_histogram(binwidth=0.005, 
                 color="black", alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="Allele frequency in controls",y="Count")+annotation_custom(grob = p1,xmin = 0.05, xmax = 0.45, ymin = 0, ymax = 9000)+xlim(NA,0.5)+guides(fill=FALSE)
dev.off()

##heatmap for Figure 3
library(ComplexHeatmap)
feature<-read.table("heatmap.txt",head=T,row.names = 1)
fea_plot<-feature

row.subsections <- c(10,9)
pdf("Figure3.pdf", width = 10, height = 9)
p1<-Heatmap(fea_plot[c(1:19),c(1:4)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = circlize::colorRamp2(c(10, 20, 50, 200, 300), c("white", "#CBCDE0", "#A4ABD6","#ACA0D2", "#8076A3"),transparency=0.2),column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),row_split = data.frame(rep(c("OR and meet PM2","OR"),row.subsections)),left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(2,3)),labels = c("OR+PM2","OR"), labels_gp = gpar(col = "black", fontsize = 12))),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.0f", fea_plot[c(1:19),c(1:4)][i, j]), x, y, gp = gpar(fontsize = 14))})

p2<-Heatmap(fea_plot[c(1:19),c(5:10)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = circlize::colorRamp2(c(0.005, 0.3, 0.8, 0.99, 1), c("white", "#eff6e0", "#aec3b0", "#84a98c", "#52796f"),transparency=0.2),column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),row_split = data.frame(rep(c("OR and meet PM2","OR"),row.subsections)),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", fea_plot[c(1:19),c(5:10)][i, j]), x, y, gp = gpar(fontsize = 14))})

pch = rep("", 19)
pch[3] = "*"
pch[6] = "**"
pch[12] = "**"
pch[15] = "**"
# color mapping for -log10(pvalue)

value = rep(18.7, 19)
value[1] = 4.33
value[2] = 4.33
value[3] = 4.33
value[4] = 4.33
value[5] = 4.33
value[11] = 0
value[12] = 2.08
value[13] = 2.08
value[14] = 2.08
value[15] = 4.33
value[16] = 4.33
value[17] = 4.33
value[18] = 4.33
value[19] = 4.33

fea_plot[c(1:19),c(11)]<-value
cols=structure(c("#fef9ef","#e5989b", "#b5838d","#bf4342"),names=c("0","2.08","4.33","18.7"))
ha = rowAnnotation(
  Evidence = anno_simple(value, col = cols, pch = pch),foo = anno_mark(at = c(3, 6, 12,15), labels = c("Moderate","Strong","Supporting","Moderate"),link_width = unit(2, "mm"), #线长度
padding = unit(2, "mm"),extend = unit(0, "mm")),annotation_name_side = "bottom")

p3<-Heatmap(fea_plot[c(1:19),c(11)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = cols,column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),row_split = data.frame(rep(c("OR and meet PM2","OR"),row.subsections)),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", fea_plot[c(1:19),c(11,11)][i, j]), x, y, gp = gpar(fontsize = 14))},right_annotation = ha)

p1+p2+p3
dev.off()

##Sankey plot in Figure 4
a<-read.table("reclassification.txt",head=T,sep="\t")

a <- a %>% mutate( change = ifelse(Classification==Ref.1, "Nochange", "upgrade"))

count1<-data.frame(a %>%
  group_by(Classification,PS4_type) %>% 
  summarize(n=n()))

  
count2<-data.frame(a %>%
  group_by(PS4_type,PS4) %>% 
  summarize(n=n()))
  
count3<-data.frame(a %>%
  group_by(PS4,Reclassification_ACMG) %>% 
  summarize(n=n()))

count4<-data.frame(a %>%
  group_by(Reclassification_ACMG,change) %>% 
  summarize(n=n()))

count5<-data.frame(a %>%
  group_by(change,Gene) %>% 
  summarize(n=n()))

  
count5<-count5[count5$change=="upgrade",]

colnames(count1) <- c("source", "target", "value")
colnames(count2) <- c("source", "target", "value")
colnames(count3) <- c("source", "target", "value")
colnames(count4) <- c("source", "target", "value")
colnames(count5) <- c("source", "target", "value")
data<-rbind(count1,count2,count3,count4,count5)

nodes <- data.frame(name=c(as.character(data$source), as.character(data$target)) %>% unique())

data$IDsource=match(data$source, nodes$name)-1 
data$IDtarget=match(data$target, nodes$name)-1

data[which(data$target=="upgrade"),'direction'] <- '1'
data[which(data$target!="upgrade"),'direction'] <- '-1'
data[which(data$source=="upgrade"),'direction'] <- '1'

color <- 'd3.scaleOrdinal() .domain(["P", "LP", "VUS", "LB", "B", "P_new", "LP_new", "VUS_new", "LB_new", "B_new","PS4=Supporting","PS4=Moderate","PS4=Strong","Nochange","upgrade","DFNB31","CDH23","COL11A2","FAM65B","GJB2","LOXHD1","MYO15A","MYO7A","OTOA","S1PR2","PCDH15","TMC1","TMPRSS3","TRIOBP","1", "-1","OR","AC"]) 
.range(["#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c","#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c","#ffd166", "#e09f3e","#a98467","#F7FBFF","#5d2a42","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#bb8588","#a6808c","#dad7cd","#b8b8d1","#5b5f97"])'

Sankey.p<-sankeyNetwork(Links = data, Nodes = nodes,
             Source = "IDsource", Target = "IDtarget",
             Value = "value", NodeID = "name", LinkGroup = 'direction',
             sinksRight=FALSE, colourScale=color, nodeWidth=40, fontSize=13, nodePadding=20)
htmlwidgets::saveWidget(Sankey.p, file="Sankey.html")
webshot("Sankey.html" , "Sankey.pdf")