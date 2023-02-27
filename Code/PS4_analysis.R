##R script for analysis in PS4 paper

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

###Define functions for evidence evaluate
#Classfying VUS into six levels
add<-function(data,start,end){
    for(i in 1:nrow(data)){
        if(data$Diagnosis[i]=="P" || data$Diagnosis[i]=="LP"){data$Diagnosis_P[i]="P"}
        else{data$Diagnosis_P[i]="NonP"}
    }
    for(i in 1:nrow(data)){
        if(data$Diagnosis[i]=="P" || data$Diagnosis[i]=="LP"){data$Diagnosis_3[i]="P/LP"}
        else if(data$Diagnosis[i]=="VUS"){data$Diagnosis_3[i]="VUS"}
        else{data$Diagnosis_3[i]="B/LB"}
    }
    for(i in 1:nrow(data)){
        if(data$Diagnosis[i]=="P" || data$Diagnosis[i]=="LP"){data$Real[i]=1}
        else{data$Real[i]=0}
    }
    for(k in 1:nrow(data)){
        data$VeryStrong[k]=0
        data$Strong[k]=0
        data$Moderate[k]=0
        data$Supporting[k]=0
        for(i in c(start:end)){
            if(data[k,i] == 1){data$VeryStrong[k]=data$VeryStrong[k]+1}
            if(data[k,i] == 2){data$Strong[k]=data$Strong[k]+1}
            if(data[k,i] == 3){data$Moderate[k]=data$Moderate[k]+1}
            if(data[k,i] == 4){data$Supporting[k]=data$Supporting[k]+1}
        }
    }
    for(i in c(start:end)){
        data <- data %>% mutate( test = ifelse(data[,i] == 1, 1, 0)) %>% plyr::rename(c("test" = paste(colnames(data)[i],"VeryStrong",sep="_")))
        data <- data %>% mutate( test = ifelse(data[,i] == 2, 1, 0)) %>% plyr::rename(c("test" = paste(colnames(data)[i],"Strong",sep="_")))
        data <- data %>% mutate( test = ifelse(data[,i] == 3, 1, 0)) %>% plyr::rename(c("test" = paste(colnames(data)[i],"Moderate",sep="_")))
        data <- data %>% mutate( test = ifelse(data[,i] == 4, 1, 0)) %>% plyr::rename(c("test" = paste(colnames(data)[i],"Supporting",sep="_")))
    }
    return(data)
}

VUS_classify<-function(data){
    for(i in 1:nrow(data)){
        if((data$VeryStrong[i]==1) || (data$Strong[i]==1 && data$Supporting[i]==1) || (data$Moderate[i]==2 && data$Supporting[i]==1) || (data$Moderate[i]==1 && data$Supporting[i]==3)){data$VUS_class[i]<-"Hot"}
        else if((data$Strong[i]==1) || (data$Moderate[i]==2) || (data$Moderate[i]==1 && data$Supporting[i]==2) || (data$Supporting[i]==4)){data$VUS_class[i]<-"Warm"}
        else if((data$Moderate[i]==1 && data$Supporting[i]==1) | (data$Supporting[i]==3)){data$VUS_class[i]<-"Tepid"}
        else if((data$Moderate[i]==1) || (data$Supporting[i]==2)){data$VUS_class[i]<-"Cool"}
        else if(data$Supporting[i]==1){data$VUS_class[i]<-"Cold"}
        else if((data$Supporting[i]==0) && (data$Moderate[i]==0)&& (data$Strong[i]==0)&& (data$VeryStrong[i]==0)){data$VUS_class[i]<-"IceCold"}
        else{data$VUS_class[i]<-"NA"}
    }
    return(data)
}


#Calculate positive Likelihood ratio and other evaluation metrics
evaluate<-function(data,start,end,dig_p){
    cover<-as.data.frame(colnames(data[start:end]))
    for(i in start:end){cover$P[i-start+1]=sum(data[,i] != "0" & data[,"Diagnosis_3"] == "P/LP")}
    for(i in start:end){cover$VUS[i-start+1]=sum(data[,i] != "0" & data[,"Diagnosis_3"] == "VUS")}
    for(i in start:end){cover$B[i-start+1]=sum(data[,i] != "0" & data[,"Diagnosis_3"] == "B/LB")}
    #info<-information.gain(Diagnosis_P~.,data=data[,c(start:end,dig_p)],unit="log2")
    #info_ratio<-gain.ratio(Diagnosis_P~.,data=data[,c(start:end,dig_p)],unit="log2")
    output<-data.frame(cover)
    output$Evidence<-rownames(output)
    colnames(output)<-c("Evidence","P/LP","VUS","B/LB")
    feature<-as.data.frame(matrix(nrow=0,ncol=21))
    for(i in 1:nrow(data)){
        for(k in c(start:end)){
            if(data[i,k]!=0){data[i,k]=1}
        }
    }
    for(i in c(start:end)){
        if(unlist(strsplit(colnames(data)[i],split="_"))[1] %in% c("PM5","PS1")){
            data2<-data[data$VEP_region=="missense_variant",]
        }
        else if(unlist(strsplit(colnames(data)[i],split="_"))[1] %in% c("PP3")){data2<-data[data$VEP_region=="missense_variant" | grepl("splice",data$VEP_region),]}
        else if(unlist(strsplit(colnames(data)[i],split="_"))[1] %in% c("PM4")){data2<-data[grepl("inframe",data$VEP_region),]}
        else if(unlist(strsplit(colnames(data)[i],split="_"))[1] %in% c("PVS1")){data2<-data[data$VEP_impact=="MODERATE" | data$VEP_impact=="HIGH" | grepl("splice",data$VEP_region),]}
        else{data2<-data}
        if(length(table(data2$Diagnosis_P,data2[,i]))!=2){
            if(table(data2$Diagnosis_P,data2[,i])[2] !=0 & table(data2$Diagnosis_P,data2[,i])[3] !=0 & table(data2$Diagnosis_P,data2[,i])[4] !=0 & table(data2$Diagnosis_P,data2[,i])[1] !=0){
                a<-BayesianLR.test(table(data2$Diagnosis_P,data2[,i])[4],table(data2$Diagnosis_P,data2[,i])[4]+table(data2$Diagnosis_P,data2[,i])[2],table(data2$Diagnosis_P,data2[,i])[1],table(data2$Diagnosis_P,data2[,i])[1]+table(data2$Diagnosis_P,data2[,i])[3],R=1000)
                TP<-table(data2$Diagnosis_P,data2[,i])[4]
                FN<-table(data2$Diagnosis_P,data2[,i])[2]
                FP<-table(data2$Diagnosis_P,data2[,i])[3]
                TN<-table(data2$Diagnosis_P,data2[,i])[1]
                total<-TP+FN+TN+FP
                #confu_final<-cbind(a$truePos,a$totalDzPos,a$trueNeg,a$totalDzNeg,a[1],a[2],a[3],a[4])
                Prevalence<- (TP + FN) / total
                Accuracy<- (TP + TN) / total
                PPV <- TP / (TP + FP)
                NPV <- TN/(FN+TN)
                F1_score<- 2*(a$statistics[1]*PPV/(a$statistics[1]+PPV))
                FNR <- FN / (FN+TP)
                FPR <- FP / (FP+TN)
                FOR <- FN / (FN+TN)
                FDR <- FP / (FP+TP)
                confu_final<-cbind(colnames(data2)[i],TP,FN,FP,TN,Accuracy,PPV,NPV,FNR,FPR,FOR,FDR,F1_score,a$statistics[1],a$statistics[2],a$posLR,a$posLR.ci[1],a$posLR.ci[2],a$negLR,a$negLR.ci[1],a$negLR.ci[2])
                feature<-rbind(feature,confu_final)
            }
        }
    }
    rownames(feature)<-feature[,1]
    #colnames(feature)<-c("case_p","all_p","control_n","all_n","posLR","posLR_LB","posLR_UB","negLR","negLR_LB","negLR_UB")
    colnames(feature)<-c("name","TP","FN","FP","TN","Accuracy","PPV","NPV","FNR","FPR","FOR","FDR","F1","Sensitivity","Specificity","posLR","posLR_LB","posLR_UB","negLR","negLR_LB","negLR_UB")
    output<-data.frame(c(output[which(output$Evidence %in% rownames(feature)),],feature))
    return(output)
}

****************************************
* Set1 analysis
* AF_case >0 and AF_control >0
****************************************
AR<-read.table("9050_variants.txt",head=T,sep="\t")
AR<-add(AR,25,25)
#Read data
ARall<-AR
#Remove feature which level=0
AR<-AR[AR$F_U!=0,]

AR<-VUS_classify(AR)
write.table(AR, file = "AR.OR.variant.txt", sep = "\t", col.names = T, row.names = F,quote=F)

for(k in c(1:20)){
    AR <- AR %>% mutate( test = ifelse(CI_LOW_gnomAD>1 & P_gnomAD<0.05 & OR_gnomAD>k, 1, 0)) %>% plyr::rename(c("test" = paste("OR_gnomAD_",k,sep="")))
}

for(k in c(1:20)){
    AR <- AR %>% mutate( test = ifelse(CI_LOW>1 & P<0.05 & OR>k, 1, 0)) %>% plyr::rename(c("test" = paste("OR_",k,sep="")))
}

AR_1<-AR[(AR$Diagnosis_3=="P/LP")|(AR$Diagnosis_3=="B/LB")|(AR$VUS_class=="Cold" | AR$VUS_class=="IceCold" | AR$VUS_class=="Cool"),]


AR_result<-evaluate(AR_1,42,61,35)

write.table(AR_result, file = "AR_result.OR.txt", sep = "\t", col.names = T, row.names = F,quote=F)

##Correlation between gnomAD and CDGC

cor.test(AR[AR$F_U!=0 & AR$gnomAD_AF_EAS!=0,]$F_U,AR[AR$F_U!=0 & AR$gnomAD_AF_EAS!=0,]$gnomAD_AF_EAS) #0.517

AR<-AR[AR$F_U!=0 & AR$gnomAD_AF_EAS!=0,]

cor.test(AR$OR_CDGC,AR$OR_gnomAD) #0.726

p1<-ggplot(AR,mapping=aes(x=F_U,y=gnomAD_AF_EAS,colour="darkred"))+xlab("Allele frequency in CDGC controls")+ylab("Allele frequency in gnomAD East Asian population")+geom_point(size=3)+theme_classic()+theme(legend.position="none")+scale_color_paletteer_d("jcolors::pal5")+annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1)+annotate("text", x = 0.1, y = 1, label = "R = 0.996",color="black",size = 6)

p2<-ggplot(AR[AR$gnomAD_AF_EAS!=1,],mapping=aes(x=OR,y=OR_gnomAD,colour="darkred"))+xlab("OR by comparing with CDGC controls")+ylab("OR by comparing with gnomAD East Asian population")+geom_point(size=3)+theme_classic()+theme(legend.position="none")+scale_color_paletteer_d("ggthemes::excel_Badge")+annotate(geom = "segment", x = 0, y = 0, xend = 150, yend = 150)+annotate("text", x = 15, y = 150, label = "R = 0.719",color="black",size = 6)+xlim(0,150)

pdf("CDGC_gnomAD.pdf",height=5,width=10)
(p1|p2)+ plot_annotation(tag_levels = 'A')
dev.off()

****************************************
* Set2 analysis
* AF_case >0 and AF_control =0
****************************************
#Read data
AR<-ARall[ARall$F_U==0,]

AR<-VUS_classify(AR)
write.table(AR, file = "AR.AC.variant.txt", sep = "\t", col.names = T, row.names = F,quote=F)


for(ac in c(2:20)){
    AR <- AR %>% mutate( test = ifelse(AC>=ac, 1, 0)) %>% plyr::rename(c("test" = paste("AC",ac,sep="_")))
}

AR_1<-AR[(AR$Diagnosis_3=="P/LP")|(AR$Diagnosis_3=="B/LB")|(AR$VUS_class=="Cold" | AR$VUS_class=="IceCold" | AR$VUS_class=="Cool"),]
AR_result<-evaluate(AR_1,41,60,35)

write.table(AR_result, file = "AR_result.AC.txt", sep = "\t", col.names = T, row.names = F,quote=F)


****************************************
* AD set analysis
* case >0 and control =0
****************************************
#Read data
ADall<-AD
#Remove feature which level=0
AD<-AD[AD$AF_ctl==0,]

for(k in c(2:10)){
    AD <- AD %>% mutate( test = ifelse(AC_case>=k, 1, 0)) %>% plyr::rename(c("test" = paste("AC_",k,sep="")))
}


AD_1<-AD[(AD$Diagnosis_3=="P/LP")|(AD$Diagnosis_3=="B/LB")|(AD$VUS_class=="Cold" | AD$VUS_class=="IceCold" | AD$VUS_class=="Cool"),]

AD_result<-evaluate(AD_1,27,35,24)

write.table(AD_result, file = "AD_result.AC.txt", sep = "\t", col.names = T, row.names = F,quote=F)
