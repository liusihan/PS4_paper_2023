##R script for analysis in PS4 paper

#load packages
library(reshape2)
library(ggsci)
library(stringr)
library(ggplot2)
library(scales)
library(dplyr)
library(paletteer)
#library(ggstatsplot)
library(gridExtra)
library(patchwork)
library(bootLR)
#library(ggthemes)
library(tidyr)
library(tidyverse)
library(viridis)
library(networkD3)
library(webshot)
library(d3Network)
library(ComplexHeatmap)
library(cowplot)

###Define functions for evidence evaluate
#Classfying VUS into six levels
add_info<-function(data){
    for(i in 1:nrow(data)){
        if(data$Classification[i]=="P" || data$Classification[i]=="LP"){data$Classification_P[i]="P"}
        else{data$Classification_P[i]="NonP"}
    }
    for(i in 1:nrow(data)){
        if(data$Classification[i]=="P" || data$Classification[i]=="LP"){data$Classification_3[i]="P/LP"}
        else if(data$Classification[i]=="VUS"){data$Classification_3[i]="VUS"}
        else{data$Classification_3[i]="B/LB"}
    }
    for(i in 1:nrow(data)){
        if(data$Classification[i]=="P" || data$Classification[i]=="LP"){data$Real[i]=1}
        else{data$Real[i]=0}
    }
    for(k in 1:nrow(data)){
        data$VeryStrong[k]=0
        data$Strong[k]=0
        data$Moderate[k]=0
        data$Supporting[k]=0
        #remove evidence for benign
        evidence<-gsub(pattern = "[B][A-Z]\\d[=][A-Z]", replacement = "", x = data$Evidence[k])
        data$VeryStrong[k]=str_count(string = evidence, pattern = "=VS")
        data$Strong[k]=str_count(string = evidence, pattern = "=S")
        data$Moderate[k]=str_count(string = evidence, pattern = "=M")
        data$Supporting[k]=str_count(string = evidence, pattern = "=P")
        }
    return(data)
}

VUS_classify<-function(data){
    for(i in 1:nrow(data)){
        if(data$Classification[i]=="VUS"){
            if((data$VeryStrong[i]==1) || (data$Strong[i]==1 && data$Supporting[i]==1) || (data$Moderate[i]==2 && data$Supporting[i]==1) || (data$Moderate[i]==1 && data$Supporting[i]==3)){data$VUS_class[i]<-"Hot"}
            else if((data$Strong[i]==1) || (data$Moderate[i]==2) || (data$Moderate[i]==1 && data$Supporting[i]==2) || (data$Supporting[i]==4)){data$VUS_class[i]<-"Warm"}
            else if((data$Moderate[i]==1 && data$Supporting[i]==1) | (data$Supporting[i]==3)){data$VUS_class[i]<-"Tepid"}
            else if((data$Moderate[i]==1) || (data$Supporting[i]==2)){data$VUS_class[i]<-"Cool"}
            else if(data$Supporting[i]==1){data$VUS_class[i]<-"Cold"}
            else if((data$Supporting[i]==0) && (data$Moderate[i]==0)&& (data$Strong[i]==0)&& (data$VeryStrong[i]==0)){data$VUS_class[i]<-"IceCold"}
            else{data$VUS_class[i]<-"NA"}
        }
        else{data$VUS_class[i]<-""}
    }
    return(data)
}

#Calculate positive Likelihood ratio and other evaluation metrics
LR<-function(data,start,end){
    output<-as.data.frame(matrix(nrow=0,ncol=21))
    for(i in start:end){
        if(length(table(data$Classification_P,data[,i]))!=2){
            if(table(data$Classification_P,data[,i])[2] !=0 & table(data$Classification_P,data[,i])[3] !=0 & table(data$Classification_P,data[,i])[4] !=0 & table(data$Classification_P,data[,i])[1] !=0){
                a<-BayesianLR.test(table(data$Classification_P,data[,i])[4],table(data$Classification_P,data[,i])[4]+table(data$Classification_P,data[,i])[2],table(data$Classification_P,data[,i])[1],table(data$Classification_P,data[,i])[1]+table(data$Classification_P,data[,i])[3],R=10000)
                TP<-table(data$Classification_P,data[,i])[4]
                FN<-table(data$Classification_P,data[,i])[2]
                FP<-table(data$Classification_P,data[,i])[3]
                TN<-table(data$Classification_P,data[,i])[1]
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
                confu_final<-cbind(colnames(data)[i],TP,FN,FP,TN,Accuracy,PPV,NPV,FNR,FPR,FOR,FDR,F1_score,a$statistics[1],a$statistics[2],a$posLR,a$posLR.ci[1],a$posLR.ci[2],a$negLR,a$negLR.ci[1],a$negLR.ci[2])
                output<-rbind(output,confu_final)
            }
        }
    }
    colnames(output)<-c("name","TP","FN","FP","TN","Accuracy","PPV","NPV","FNR","FPR","FOR","FDR","F1","Sensitivity","Specificity","posLR","posLR_LB","posLR_UB","negLR","negLR_LB","negLR_UB")
    return(output)
}

###define functions for calculate local LR
local_lr<-function(input_data,feature,op_vs,alpha,minpoints,increment){
    op_s<-op_vs^(2^-1)
    op_m<-op_vs^(2^-2)
    op_p<-op_vs^(2^-3)
    postp_vs<-op_vs*alpha/((op_vs-1)*alpha+1)
    postp_s<-op_s*alpha/((op_s-1)*alpha+1)
    postp_m<-op_m*alpha/((op_m-1)*alpha+1)
    postp_p<-op_p*alpha/((op_p-1)*alpha+1)
    print(c(op_vs,op_s,op_m,op_p,postp_vs,postp_s,postp_m,postp_p))
    w<-(1-alpha)*(sum(input_data$Classification_P=="P"))/(alpha*(sum(input_data$Classification_P=="NonP")))
    thrs<-sort(unique(c(input_data[[feature]],floor(min(input_data[[feature]])),ceiling(max(input_data[[feature]])))))
    post <- rep(0, times=length(thrs))
    finalwindow <- rep(0, times=length(thrs))
    maxthrs<-ceiling(max(input_data[[feature]]))
    minthrs<-floor(min(input_data[[feature]]))
    output<-as.data.frame(matrix(nrow=0,ncol=3))
    for(i in c(1:length(thrs))){
        halfwindow = 0
        lo <- thrs[i] - halfwindow
        hi <- thrs[i] + halfwindow
        while (1){
            pos <- length(input_data[input_data[[feature]]>lo&input_data[[feature]]<hi&input_data$Classification_P=="P",][[feature]])
            neg <- length(input_data[input_data[[feature]]>lo&input_data[[feature]]<hi&input_data$Classification_P=="NonP",][[feature]])
            if (hi > maxthrs){c <- (maxthrs - lo) / (hi - lo)}
            else if(lo < minthrs){c <- (hi - minthrs) / (hi - lo)}
            else{c<-1}
            if(c<=0){break}
            if(pos + neg < c * minpoints){
                halfwindow <- halfwindow + increment
                lo <- thrs[i] - halfwindow
                hi <- thrs[i] + halfwindow
            }
            if(pos + neg >= c * minpoints){
                break
            }
        }
        post[i] <- pos / (pos + w * neg)
        finalwindow[i] <- halfwindow
        temp<-cbind(thrs[i],post[i],finalwindow[i])
        output<-rbind(output,temp)
    }
    colnames(output)<-c("thrs","post_p","finalwindow")
    return(output)
}

local_bootstrapped_lr<-function(input_data,feature,op_vs,alpha,bootstrap,minpoints,increment){
    for(i in c(1:bootstrap)){
        idx=sample(1:nrow(input_data),replace=T)
        sample_data<-input_data[idx,]
        bootstrap_i<-local_lr(sample_data,feature,op_vs,alpha,minpoints,increment)
        write.table(bootstrap_i, file = paste("bootstrap_",i,".txt",sep=""), sep = "\t", col.names = T, row.names = F,quote=F)
    }
}

get_threshold<-function(postp_list,discountonesided,bootstrap){
    thresh<-as.data.frame(matrix(nrow=bootstrap,ncol=4))
    DiscountedThreshold<-c(0,0,0,0)
    for(i in c(1:bootstrap)){
        input<-read.table(paste("bootstrap_",i,".txt",sep=""),sep = "\t",header=T)
        for(j in c(1:length(postp_list))){
            ind = min(which(input$post_p>=postp_list[j]))
            if(ind==1){thresh[i,j] = input$thrs[2]}
            else if(ind>1){thresh[i,j] = input$thrs[ind]}
            else(thresh[i,j] = NA)
        }
    }
    for(j in c(1:length(postp_list))){
        invalids<-sum(is.na(thresh[,j]))
        if(invalids > (discountonesided * bootstrap)){DiscountedThreshold[j] = NA}
        else{
            t=sort(thresh[which(!is.na(thresh[,j])),j],decreasing = TRUE)
            DiscountedThreshold[j] = t[floor(discountonesided * (bootstrap-invalids)) + 1]
        }
    }
    return(DiscountedThreshold)
}

CI_95<-function(bootstrap){
    postp<-read.table("OR_raw.txt",sep = "\t",header=T)
    postp_matrix<-postp[,c(1:2)]
    for(i in 1:bootstrap){
        input<-read.table(paste("bootstrap_",i,".txt",sep=""),sep = "\t",header=T)
        postp_matrix<-left_join(postp_matrix,input[,c(1:2)],by="thrs")
        colnames(postp_matrix)[i+2]<-paste("postp_",i,sep="")
    }
    output<-as.data.frame(matrix(nrow=nrow(postp_matrix),ncol=3))
    colnames(output)<-c("OR","Posterior","Posterior1")
    output[,1]<-postp_matrix[,1]
    output[,2]<-postp_matrix[,2]
    for(j in 1:nrow(postp_matrix)){
        a<-t.test(postp_matrix[j,c(3:ncol(postp_matrix))])
        output[j,3]<-a$conf.int[1]
    }
    return(output)
}

Classification<-function(data){
  data$Classification_new=""
  for(k in 1:nrow(data)){
    #remove evidence for benign
    p_evidence<-gsub(pattern = "[B][A-Z]\\d[=][A-Z]", replacement = "", x = data$Evidence[k])
    #remove evidence for pathogenic
    b_evidence<-gsub(pattern = "[P][A-Z]\\d[=][A-Z]", replacement = "", x = data$Evidence[k])
    p_VeryStrong=str_count(string = p_evidence, pattern = "=VS")
    p_Strong=str_count(string = p_evidence, pattern = "=S")
    p_Moderate=str_count(string = p_evidence, pattern = "=M")
    p_Supporting=str_count(string = p_evidence, pattern = "=P")
    Stand_alone=str_count(string = b_evidence, pattern = "=A")
    b_Strong=str_count(string = b_evidence, pattern = "=S")
    b_Moderate=str_count(string = b_evidence, pattern = "=M")
    b_Supporting=str_count(string = b_evidence, pattern = "=P")
    PVS1=str_count(string = p_evidence, pattern = "PVS1=VS")
    PM2=str_count(string = p_evidence, pattern = "PM2=P")
    p_class=""
    b_class=""
    if(p_VeryStrong==1 && p_Moderate==1){p_class="LP"}
    if(p_Strong==1 && p_Moderate>=1){p_class="LP"}
    if(p_Strong==1 && p_Supporting>=2){p_class="LP"}
    if(p_Moderate>=3){p_class="LP"}
    if(p_Moderate==2 && p_Supporting>=2){p_class="LP"}
    if(p_Moderate==1 && p_Supporting>=4){p_class="LP"}
    if(PVS1==1 && PM2==1){p_class="LP"}
    if(p_VeryStrong==1){
      if(p_Strong>=1 || p_Moderate>=2 || p_Supporting>=2 || (p_Moderate==1 && p_Supporting==1)){
        p_class="P"
      }
    }
    if(p_Strong>=2){p_class="P"}
    if(p_Strong>=1){
      if(p_Moderate>=3 || (p_Moderate>=2 && p_Supporting>=2) || (p_Moderate>=1 && p_Supporting>=4)){
        p_class="P"
      }
    }
    if((b_Strong==1 && b_Supporting>=1) || (b_Supporting>=2) || (b_Strong==1 && b_Moderate==1) || ((b_Supporting>=1 && b_Moderate==1))){b_class="LB"}
    if(Stand_alone==1 || b_Strong>=2){b_class="B"}
    if(b_class=="" && p_class==""){data$Classification_new[k]="VUS"}
    else if(b_class!="" && p_class!=""){data$Classification_new[k]="VUS"}
    else if(b_class=="" && p_class!=""){data$Classification_new[k]=p_class}
    else{data$Classification_new[k]=b_class}
  }
  return(data)
}

reclassification<-function(data){
  data$PS4_level=""
  data$PS4_type=""
  data$Reclassification=""
  data$Evidence_combine=""
  for(k in 1:nrow(data)){
    if(data$AF_case[k]>=0.0005 & data$AF_control[k]>0){
      if(data$OR[k]>6 & data$CI_LOW[k]>=1){data$PS4_level[k]="PS4_Strong";data$PS4_type[k]="PS4_OR";data$Evidence_combine[k]=paste(data$Evidence[k],";PS4=S",sep="")}
      else if(data$OR[k]>3 & data$CI_LOW[k]>=1){data$PS4_level[k]="PS4_Moderate";data$PS4_type[k]="PS4_OR";data$Evidence_combine[k]=paste(data$Evidence[k],";PS4=M",sep="")}
      else{data$Evidence_combine[k]=data$Evidence[k]}
    }
    else if(data$AF_case[k]<0.0005 & data$AF_control[k]>0){
      if(data$OR[k]>=2.27){data$PS4_level[k]="PS4_Supporting";data$PS4_type[k]="PS4_OR";data$Evidence_combine[k]=paste(data$Evidence[k],";PS4=P",sep="")}
      else{data$Evidence_combine[k]=data$Evidence[k]}
    }
    else{
      if(data$AC_case[k]>=6){data$PS4_level[k]="PS4_Moderate";data$PS4_type[k]="PS4_AC";data$Evidence_combine[k]=paste(data$Evidence[k],";PS4=M",sep="")}
      else if(data$AC_case[k]>=3){data$PS4_level[k]="PS4_Supporting";data$PS4_type[k]="PS4_AC";data$Evidence_combine[k]=paste(data$Evidence[k],";PS4=P",sep="")}
      else{data$Evidence_combine[k]=data$Evidence[k]}
    }
    #remove evidence for benign
    p_evidence<-gsub(pattern = "[B][A-Z]\\d[=][A-Z]", replacement = "", x = data$Evidence_combine[k])
    #remove evidence for pathogenic
    b_evidence<-gsub(pattern = "[P][A-Z]\\d[=][A-Z]", replacement = "", x = data$Evidence_combine[k])
    p_VeryStrong=str_count(string = p_evidence, pattern = "=VS")
    p_Strong=str_count(string = p_evidence, pattern = "=S")
    p_Moderate=str_count(string = p_evidence, pattern = "=M")
    p_Supporting=str_count(string = p_evidence, pattern = "=P")
    Stand_alone=str_count(string = b_evidence, pattern = "=A")
    b_Strong=str_count(string = b_evidence, pattern = "=S")
    b_Moderate=str_count(string = b_evidence, pattern = "=M")
    b_Supporting=str_count(string = b_evidence, pattern = "=P")
    PVS1=str_count(string = p_evidence, pattern = "PVS1=VS")
    PM2=str_count(string = p_evidence, pattern = "PM2=P")
    p_class=""
    b_class=""
    if(p_VeryStrong==1 && p_Moderate==1){p_class="LP"}
    if(p_Strong==1 && p_Moderate>=1){p_class="LP"}
    if(p_Strong==1 && p_Supporting>=2){p_class="LP"}
    if(p_Moderate>=3){p_class="LP"}
    if(p_Moderate==2 && p_Supporting>=2){p_class="LP"}
    if(p_Moderate==1 && p_Supporting>=4){p_class="LP"}
    if(PVS1==1 && PM2==1){p_class="LP"}
    if(p_VeryStrong==1){
      if(p_Strong>=1 || p_Moderate>=2 || p_Supporting>=2 || (p_Moderate==1 && p_Supporting==1)){
        p_class="P"
      }
    }
    if(p_Strong>=2){p_class="P"}
    if(p_Strong>=1){
      if(p_Moderate>=3 || (p_Moderate>=2 && p_Supporting>=2) || (p_Moderate>=1 && p_Supporting>=4)){
        p_class="P"
      }
    }
    if((b_Strong==1 && b_Supporting>=1) || (b_Supporting>=2) || (b_Strong==1 && b_Moderate==1) || ((b_Supporting>=1 && b_Moderate==1))){b_class="LB"}
    if(Stand_alone==1 || b_Strong>=2){b_class="B"}
    if(b_class=="" && p_class==""){data$Reclassification[k]="VUS"}
    else if(b_class!="" && p_class!=""){data$Reclassification[k]="VUS"}
    else if(b_class=="" && p_class!=""){data$Reclassification[k]=p_class}
    else{data$Reclassification[k]=b_class}
  }
  return(data)
}

****************************************
* Correlation between gnomAD_EAS data 
* and CDGC controls
****************************************
AR<-read.table("9050.variants.txt",head=T,sep="\t")
AR<-add_info(AR)
AR<-VUS_classify(AR)
write.table(AR, file = "9050.variants.allinfo.txt", sep = "\t", col.names = T, row.names = F,quote=F)
AR$Classification_3<-factor(AR$Classification_3,levels=c("P/LP","VUS","B/LB"))

cor.test(AR[AR$AF_control!=0 & AR$gnomAD_AF_EAS!=0,]$AF_control,AR[AR$AF_control!=0 & AR$gnomAD_AF_EAS!=0,]$gnomAD_AF_EAS) #0.517

cor.test(AR[AR$AF_control!=0 & AR$gnomAD_AF_EAS!=0,]$OR,AR[AR$AF_control!=0 & AR$gnomAD_AF_EAS!=0,]$OR_gnomAD) #0.726

****************************************
* PM2 and PP3 evaluation
****************************************
AR<-mutate(AR,PM2=ifelse(str_detect(AR$Evidence,"PM2"),1,0),PP3_P=ifelse(str_detect(AR$Evidence,"PP3"),1,0),PP3_M=ifelse(str_detect(AR$Evidence,"PP3=M"),1,0))
#PM2
truth_set<-AR[AR$VUS_class!="Hot" & AR$VUS_class!="Warm" & AR$VUS_class!="Tepid" & AR$Database!="Conflicting",]

PM2_result<-LR(truth_set,38,38)
write.table(PM2_result, file = "PM2.LR.txt", sep = "\t", col.names = T, row.names = F,quote=F)

PP3<-truth_set[(truth_set$VEP_consequence=="missense_variant"),]
PP3<-PP3[(PP3$Database!="" | PP3$Classification_P=="NonP"),]
PP3_result<-LR(PP3,39,40)

write.table(PP3_result, file = "PP3.LR.txt", sep = "\t", col.names = T, row.names = F,quote=F)

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

****************************************
* Truth set1 analysis
* AF_case >=0.0005 and AF_control >0
****************************************
#Extracted variants with AF_case >=0.0005 and AF_control >0
truth_subset1<-truth_set[truth_set$AF_control!=0 & truth_set$AF_case>=0.0005,]
write.table(truth_subset1, file = "truth_subset1.txt", sep = "\t", col.names = T, row.names = F,quote=F)

#Adding the tag if a variant passed the OR cutoff
for(k in c(1:10)){
  truth_subset1 <- truth_subset1 %>% mutate( test = ifelse(CI_LOW>1 & OR>k, 1, 0)) %>% plyr::rename(c("test" = paste("OR_",k,sep="")))
}

# calculate positive LR
subset1_result<-LR(truth_subset1,41,50)

write.table(subset1_result, file = "subset1_result.LR.txt", sep = "\t", col.names = T, row.names = F,quote=F)

****************************************
* Truth set 2 analysis
* AF_case <0.0005 and AF_control >0
****************************************
truth_subset2<-truth_set[truth_set$AF_control!=0 & truth_set$AF_case<0.0005,]
write.table(truth_subset2, file = "truth_subset2.txt", sep = "\t", col.names = T, row.names = F,quote=F)

truth_subset2$Classification_P<-factor(truth_subset2$Classification_P,levels=c("P","NonP"))

setwd("./PS4_lr")
OR_raw<-local_lr(truth_subset2,"OR",387,186/1997,100,0.01)
write.table(OR_raw, file = "OR_raw.txt", sep = "\t", col.names = T, row.names = F,quote=F)

local_bootstrapped_lr(truth_subset2,"OR",378,194/2032,10000,100,0.01)
postp_list<-c(0.1778350,0.3129676,0.6689245,0.9754584)
get_threshold(postp_list,0.05,10000)

thresh<-CI_95(10000)
write.table(thresh, file = "lr_CI.txt", sep = "\t", col.names = T, row.names = F,quote=F)

****************************************
* Truth set 3 analysis
* AF_control =0
****************************************
truth_subset3<-truth_set[truth_set$AF_control==0,]
write.table(truth_subset3, file = "truth_subset3.txt", sep = "\t", col.names = T, row.names = F,quote=F)
for(ac in c(2:10)){
  truth_subset3 <- truth_subset3 %>% mutate( test = ifelse(AC_case>=ac, 1, 0)) %>% plyr::rename(c("test" = paste("AC",ac,sep="_")))
}

truth_subset3_result<-LR(truth_subset3,41,49)

write.table(truth_subset3_result, file = "subset3_result.LR.txt", sep = "\t", col.names = T, row.names = F,quote=F)

##reclassification
AR_reclass<-reclassification(AR)
write.table(AR_reclass, file = "9050.variants.reclassification.txt", sep = "\t", col.names = T, row.names = F,quote=F)
