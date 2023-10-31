##R script for plot in PS4 paper

#load packages
library(reshape2)
library(ggsci)
library(ggplot2)
library(scales)
library(dplyr)
library(paletteer)
library(ggpie)
library(gridExtra)
library(patchwork)
library(bootLR)
library(tidyr)
library(tidyverse)
#library(viridis)
library(patchwork)
library(networkD3)
library(webshot)
library(d3Network)
library(ComplexHeatmap)
library(cowplot)

*****************************************
  * Figure 2, Figure S2 and Figure S3
*****************************************

AR<-read.table("9050.variants.allinfo.txt",header = T,sep="\t")
AR$Classification<-factor(AR$Classification,levels = c('P', 'LP', 'VUS','LB','B'))

Figure2A<-ggdonut(data = AR, group_key = "Classification", count_type = "full",
                        label_info = "all", label_type = "horizon",
                        label_size = 6, label_pos = "out")+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(fill = "Variant Classification")
##ggpie(data = AR, group_key = "Classification", count_type = "full",
       #label_info = "all", label_type = "horizon",
       #label_size = 4, label_pos = "out")+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))

#AR1<-ggplot(pie_input, aes(x="", y=Freq, fill=Classification))+geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+ blank_theme +theme(axis.text.x=element_blank())+guides(fill=FALSE)

VUS<-data.frame(table(AR[AR$VUS_class!="",]$VUS_class))
colnames(VUS)<-c("VUS_temperature_scale","Number of variants")
VUS_data<- data.frame(matrix(,nrow = 6, ncol = 0))
VUS_data$Posterior_probability<-c("[81.2%-90)","[67.5%-81.2%)","[50%-67.5%)","[32.5%-50%)","[18.8%-32.5%)","[10%-18.8%)")
VUS_data$VUS_temperature_scale<-c("Hot","Warm","Tepid","Cool","Cold","IceCold")
VUS_data<-full_join(VUS_data,VUS)
rownames(VUS_data)<-VUS_data$VUS_temperature_scale

P<-as.data.frame(table(AR[AR$Classification=="P"|AR$Classification=="LP",]$Gene))
P<-data.frame(AR %>%
                group_by(Classification,Gene) %>%
                summarize(n=n()) %>%
                arrange(desc(n),Gene))
P<-P[P$Classification=="P"|P$Classification=="LP",]
data<-P
gene<-data.frame(data %>%
                   group_by(Gene) %>%
                   summarize(sum(n)))
gene <- arrange(gene,desc(sum.n.))
data$Gene<-factor(data$Gene,levels=as.character(unique(gene$Gene)))
data$Classification<-factor(data$Classification,levels=c("P","LP"))
Figure2D<-ggplot(data=data[data$Gene%in%gene$Gene[1:20],], mapping=aes(x=Gene,y=n,fill=Classification))+
  geom_bar(stat="identity",colour = 'black', width = 1, lwd=0.1)+ theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8,face="italic"))+labs(fill="", x = "", y = "Number of P/LP variants")+scale_fill_manual(values=c("#780000","#c1121f"))+guides(fill=FALSE)


p1<-ggplotGrob(ggplot(data=AR,aes(x=AF_case,fill=Classification)) + 
                 geom_histogram(binwidth=0.00005, 
                                color="black", alpha=0.9)+xlim(NA,0.001)+ylim(NA,6000)+xlab("Minor allele frequency in cases")+
                 theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="",y="")+guides(fill=FALSE))

Figure2E<-ggplot(data=AR,aes(x=AF_case,fill=Classification)) + 
  geom_histogram(binwidth=0.0003, 
                 color="black", alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="Allele frequency in cases",y="Count")+annotation_custom(grob = p1,xmin = 0.0001, xmax = 0.01, ymin = 0, ymax = 7000)+xlim(NA,0.01)+guides(fill=FALSE)

p1<-ggplotGrob(ggplot(data=AR,aes(x=AF_control,fill=Classification)) + 
                 geom_histogram(binwidth=0.00005, 
                                color="black", alpha=0.9)+xlim(NA,0.001)+ylim(NA,6000)+xlab("MAF")+
                 theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="",y="")+guides(fill=FALSE))

Figure2F<-ggplot(data=AR,aes(x=AF_control,fill=Classification)) + 
  geom_histogram(binwidth=0.0003, 
                 color="black", alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+labs(x="Allele frequency in controls",y="Count")+annotation_custom(grob = p1,xmin = 0.0001, xmax = 0.01, ymin = 0, ymax = 7000)+xlim(NA,0.01)+guides(fill=FALSE)

cover<-data.frame(AR %>%
                    group_by(Classification,Classification_3,VEP_consequence) %>%
                    summarize(n=n()) %>%
                    arrange(desc(n),VEP_consequence))
count<-data.frame(AR %>%
                    group_by(Classification_3) %>%
                    summarize(n=n()))
cover<-full_join(cover,count,by=c("Classification_3"))
cover<-cover %>% mutate(percent=n.x/n.y)
cover2<-data.frame(AR %>%
                     group_by(Classification_3,VEP_consequence) %>% 
                     summarize(n=n()) %>%
                     arrange(desc(n),VEP_consequence))
count<-data.frame(AR %>%
                    group_by(Classification_3) %>% 
                    summarize(n=n()))
cover2<-full_join(cover2,count,by=c("Classification_3"))
cover2<-cover2 %>% mutate(percent=n.x/n.y)
cover$VEP_consequence<-factor(cover$VEP_consequence,levels=as.character(unique(c(cover2[cover2$Classification_3=="P/LP",]$VEP_consequence,cover2[cover2$Classification_3=="VUS",]$VEP_consequence,cover2[cover2$Classification_3=="B/LB",]$VEP_consequence))))
cover$Classification_3<-factor(cover$Classification_3,levels=c("P/LP","VUS","B/LB"))
cover$Classification<-factor(cover$Classification,levels = c('P', 'LP', 'VUS','LB','B'))
cover3<-data.frame(AR %>%
                     group_by(Classification_3,VEP_consequence) %>% 
                     summarize(n=n()) %>%
                     arrange(desc(n),VEP_consequence))
cover3<-full_join(cover3,count,by=c("Classification_3"))
cover3<-cover3 %>% mutate(percent=n.x/n.y)
cover3$Classification_3<-factor(cover3$Classification_3,levels=c("P/LP","VUS","B/LB"))
cover3$VEP_consequence<-factor(cover3$VEP_consequence,levels=as.character(unique(c(cover3[cover3$Classification_3=="P/LP",]$VEP_consequence,cover3[cover3$Classification_3=="VUS",]$VEP_consequence,cover3[cover3$Classification_3=="B/LB",]$VEP_consequence))))
cover4<-full_join(cover,cover3,by=c("Classification_3","VEP_consequence"))
cover4$number<-ifelse(duplicated(paste(cover4$Classification_3,cover4$VEP_consequence,sep="_")),"",cover4$n.x.y)
Figure2C<-ggplot(cover4,aes(x=VEP_consequence,y=percent.x,fill=Classification))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual(values=c("#780000","#c1121f","#bfdbf7","#cad2c5","#84a98c"))+guides(fill=guide_legend(title=NULL))+theme(axis.title = element_blank(),legend.position='none',axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10))+facet_grid(Classification_3~.)+ylim(0,1)+guides(fill=FALSE)+ geom_text(aes(x=VEP_consequence,y=percent.y,label = number),size=3,vjust = 0, nudge_y = 0.01)

pdf("Figure2.pdf",height=9,width=9)
((Figure2A+gridExtra::tableGrob(VUS_data[,c(1,3)])+plot_layout(widths = c(1, 2))) / (Figure2C | Figure2D) / (Figure2E|Figure2F))+plot_annotation(tag_levels = 'A')
dev.off()

pdf("Figure2A.pdf")
Figure2A
dev.off()

##Figure 3
#Simulation to get the allele frequency threshold that variants have a CI_low >1
simu_result<-as.data.frame(matrix(nrow=0,ncol=11))
for(sample in c(100,200,300,400,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,12000,13845,15000,20000,50000,100000)){
  for(ac_case in c(1:(sample*0.1))){
    v1<-ac_case
    v2<-sample*2
    v3<-1
    v4<-6570*2
    fi_matrix<-matrix(c(v1,v2,v3,v4),ncol=2)
    result<-fisher.test(fi_matrix)
    af_case<-v1/v2
    af_contorl<-v3/v4
    pvalue<-result$p.value
    OR<-result$estimate
    OR_CI<-result$conf.int[1]
    OR_UP<-result$conf.int[2]
    if(OR_CI>1){
      out<-cbind(sample,v1,v2,v3,v4,af_case,af_contorl,pvalue,OR,OR_CI,OR_UP)
      simu_result<-rbind(simu_result,out)
      break
    }
  }
}
Figure3A<-ggplot(data=simu_result, mapping=aes(x=sample,y=af_case))+geom_point(size=2,color="#386fa4")+theme_classic()+xlab("Sample size of case")+ylab("Allele frequency")+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+geom_hline(yintercept = 0.0005,colour="red",linetype=2)+annotate("rect",xmin=7000,xmax=20000,ymin=0,ymax=0.01,alpha=.7,fill="#cad2c5")

##heatmap for truth subset 1
subset1_result<-read.table("subset1_result.LR.txt",header = T,sep="\t",row.names = 1)
fea_plot<-subset1_result[,c(1:7,13,14,12,16)]

p1<-Heatmap(fea_plot[c(1:10),c(1:4)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = circlize::colorRamp2(c(10, 20, 50, 200, 300), c("white", "#CBCDE0", "#A4ABD6","#ACA0D2", "#8076A3"),transparency=0.2),column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.0f", fea_plot[c(1:10),c(1:4)][i, j]), x, y, gp = gpar(fontsize = 14))})

p2<-Heatmap(fea_plot[c(1:10),c(5:10)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = circlize::colorRamp2(c(0.005, 0.3, 0.8, 0.99, 1), c("white", "#eff6e0", "#aec3b0", "#84a98c", "#52796f"),transparency=0.2),column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", fea_plot[c(1:10),c(5:10)][i, j]), x, y, gp = gpar(fontsize = 14))})

value = rep(0, 10)
LR_list<-c(1.98,3.91,15.26,233)
cols=structure(rep("#fef9ef", 10),names=fea_plot[c(1:10),c(11)])
strength_cols=structure(c("#e5989b", "#b5838d","#bf4342","#780000"),names=c(1,2,3,4))
for(i in 1:length(fea_plot[,11])){
  for(k in 1:length(LR_list)){
    if(fea_plot[,11][i]>=LR_list[k])
    {
      value[i]=LR_list[k]
      cols[i]=strength_cols[k]
    }
  }
}

pch = rep("", 10)
for(k in 1:length(LR_list)){
  if(min(which(value==LR_list[k]))!="Inf"){pch[min(which(value==LR_list[k]))]="*"}
}

anno_cols=structure(c("#fef9ef","#e5989b", "#b5838d","#bf4342","#780000"),names=c(0,LR_list))
ha = rowAnnotation(
  Evidence = anno_simple(value, col = anno_cols, pch = pch),foo = anno_mark(at = c(3, 6), labels = c("Moderate","Strong"),link_width = unit(2, "mm"), #线长度
                                                                       padding = unit(2, "mm"),extend = unit(0, "mm")),annotation_name_side = "bottom")

p3<-Heatmap(fea_plot[c(1:10),c(11)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = cols,column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", fea_plot[c(1:10),c(11,11)][i, j]), x, y, gp = gpar(fontsize = 14))},right_annotation = ha)

Figure3B<-p1+p2+p3

thresh<-read.table("lr_CI.txt",header = T,sep="\t")
postp_list<-c(0.1778350,0.3129676,0.6689245,0.9754584)
Figure3C<-ggplot() + geom_line(data = thresh, mapping = aes(x = OR, y = Posterior), colour = "black",size=1)+ geom_line(data = thresh, mapping = aes(x = OR, y = Posterior1), colour = "grey",size=1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+geom_hline(yintercept = postp_list[4],colour="red",linetype=1,size=1)+geom_hline(yintercept = postp_list[3],colour="red",linetype=2,size=1)+geom_hline(yintercept = postp_list[2],colour="red",linetype=3,size=1)+geom_hline(yintercept = postp_list[1],colour="red",linetype=4,size=1)

subset3_result<-read.table("subset3_result.LR.txt",header = T,sep="\t",row.names = 1)
fea_plot<-subset3_result[c(1:9),c(1:7,13,14,12,16)]

p1<-Heatmap(fea_plot[c(1:9),c(1:4)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = circlize::colorRamp2(c(10, 20, 50, 200, 300), c("white", "#CBCDE0", "#A4ABD6","#ACA0D2", "#8076A3"),transparency=0.2),column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.0f", fea_plot[c(1:9),c(1:4)][i, j]), x, y, gp = gpar(fontsize = 14))})

p2<-Heatmap(fea_plot[c(1:9),c(5:10)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = circlize::colorRamp2(c(0.005, 0.3, 0.8, 0.99, 1), c("white", "#eff6e0", "#aec3b0", "#84a98c", "#52796f"),transparency=0.2),column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", fea_plot[c(1:9),c(5:10)][i, j]), x, y, gp = gpar(fontsize = 14))})

value = rep(0, 9)
LR_list<-c(2.11,4.46,19.90,397)
cols=structure(rep("#fef9ef", 9),names=fea_plot[c(1:9),c(11)])
strength_cols=structure(c("#e5989b", "#b5838d","#bf4342","#780000"),names=c(1,2,3,4))
for(i in 1:length(fea_plot[,11])){
  for(k in 1:length(LR_list)){
    if(fea_plot[,11][i]>=LR_list[k])
    {
      value[i]=LR_list[k]
      cols[i]=strength_cols[k]
    }
  }
}

pch = rep("", 9)
for(k in 1:length(LR_list)){
  if(min(which(value==LR_list[k]))!="Inf"){pch[min(which(value==LR_list[k]))]="*"}
}

anno_cols=structure(c("#fef9ef","#e5989b", "#b5838d","#bf4342","#780000"),names=c(0,LR_list))
ha = rowAnnotation(
  Evidence = anno_simple(value, col = anno_cols, pch = pch),foo = anno_mark(at = c(2, 3), labels = c("Moderate","Strong"),link_width = unit(2, "mm"), #线长度
                                                                            padding = unit(2, "mm"),extend = unit(0, "mm")),annotation_name_side = "bottom")

p3<-Heatmap(fea_plot[c(1:9),c(11)],cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "black", lwd = 1.5),col = cols,column_dend_height = unit(2, "cm"),row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.3f", fea_plot[c(1:9),c(11,11)][i, j]), x, y, gp = gpar(fontsize = 14))},right_annotation = ha)

Figure3D<-p1+p2+p3

pdf("Figure3.pdf",height=9,width=12)
(Figure3A+plot_spacer() ) / (Figure3C+plot_spacer() )+plot_annotation(tag_levels = 'A')+plot_layout(widths = c(1, 2))
dev.off()

pdf("Figure3B.pdf",height=8,width=12)
Figure3B+plot_annotation(tag_levels = 'A')
dev.off()

pdf("Figure3D.pdf",height=8,width=12)
Figure3D+plot_annotation(tag_levels = 'A')
dev.off()
