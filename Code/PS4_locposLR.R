#load package
library(reshape2)
library(ggsci)
library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(gridExtra)
library(patchwork)
library(FSelector)
library(caret)
library(bootLR)
library(ggthemes)
library(tidyr)
library(yardstick)
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
    w<-(1-alpha)*(sum(input_data$Diagnosis=="P"))/(alpha*(sum(input_data$Diagnosis=="NonP")))
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
            pos <- length(input_data[input_data[[feature]]>lo&input_data[[feature]]<hi&input_data$Diagnosis=="P",][[feature]])
            neg <- length(input_data[input_data[[feature]]>lo&input_data[[feature]]<hi&input_data$Diagnosis=="NonP",][[feature]])
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
postp_list<-c(0.182844,0.318783,0.671789,0.975098)
get_threshold<-function(postp_list,discountonesided,bootstrap){
    thresh<-as.data.frame(matrix(nrow=bootstrap,ncol=4))
    DiscountedThreshold<-c(0,0,0,0)
    for(i in c(1:bootstrap)){
        input<-read.table(paste("bootstrap_",i,".txt",sep=""),sep = "\t",header=T)
        for(j in c(1:length(postp_list))){
            ind = min(which(input$post_p>postp_list[j]))
            if(ind>1){thresh[i,j] = input$thrs[ind]}
            else(thresh[i,j] = NA)
        }
    }
    for(j in c(1:length(postp_list))){
        invalids<-sum(is.na(thresh[,j]))
        if(invalids > (discountonesided * bootstrap)){DiscountedThreshold[j] = NA}
        else{
            t=sort(thresh[which(!is.na(thresh[,j])),j],decreasing = TRUE)
            DiscountedThreshold[j] = t[floor(discountonesided * bootstrap) - invalids + 1]
        }
    }
    return(DiscountedThreshold)
}

CI_95<-function(postp_list,CI,bootstrap){
    CI<-as.data.frame(matrix(nrow=955,ncol=bootstrap))
    CI_95<-as.data.frame(matrix(nrow=955,ncol=1))
    for(i in c(1:bootstrap)){
        input<-read.table(paste("bootstrap_",i,".txt",sep=""),sep = "\t",header=T)
        for(j in c(1:955)){
            CI[j,i]=input[j,2]
        }
    }
    for(j in c(1:955)){
        CI_95[j,1]<-unname(quantile(CI[j,], c(.95)))
    }
    return(CI_95)
}

****************************************
* Truth set 2 analysis
* AF_case <0.0005 and AF_control >0
****************************************
AR<-read.table("9050.variants.txt",head=T,sep="\t")
AR<-AR[AR$F_U!=0 & AR$F_A<0.0005,]
AR<-VUS_classify(AR)
AR2<-AR[(AR$Diagnosis_3=="P/LP")|(AR$Diagnosis_3=="B/LB")|(AR$VUS_class=="Cold" | AR$VUS_class=="IceCold" | AR$VUS_class=="Cool"),]
local_bootstrapped_lr(AR2,"OR",366,199/2059,10000,100,0.01)
get_threshold(postp_list,0.05,1000)
OR_raw<-local_lr(AR2,"OR",366,199/2059,100,0.01)
or<-read.table("OR.result.txt",sep="\t",head=T)
thresh<-as.data.frame(matrix(nrow=nrow(or),ncol=3))
thresh[,1]<-or[,1]
thresh[,2]<-or[,2]
for(i in c(1:nrow(or))){
    a<-t.test(or[i,c(3:ncol(or))])
    thresh[i,3]<-a$conf.int[1]
}
colnames(thresh)<-c("OR","Posterior","Posterior1")

p1<-ggplot() + geom_line(data = thresh, mapping = aes(x = OR, y = Posterior), colour = "black",size=1)+ geom_line(data = thresh, mapping = aes(x = OR, y = Posterior1), colour = "grey",size=1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+geom_hline(yintercept = postp_list[4],colour="red",linetype=1,size=1)+geom_hline(yintercept = postp_list[3],colour="red",linetype=2,size=1)+geom_hline(yintercept = postp_list[2],colour="red",linetype=3,size=1)+geom_hline(yintercept = postp_list[1],colour="red",linetype=4,size=1)

ef<-read.table("EF.result.txt",sep="\t",head=T)
thresh<-as.data.frame(matrix(nrow=nrow(ef),ncol=3))
thresh[,1]<-ef[,1]
thresh[,2]<-ef[,2]
for(i in c(1:nrow(ef))){
    a<-t.test(ef[i,c(3:ncol(ef))])
    thresh[i,3]<-a$conf.int[1]
}
colnames(thresh)<-c("EF","Posterior","Posterior1")

p2<-ggplot() + geom_line(data = thresh, mapping = aes(x = EF, y = Posterior), colour = "black",size=1)+ geom_line(data = thresh, mapping = aes(x = EF, y = Posterior1), colour = "grey",size=1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+geom_hline(yintercept = postp_list[4],colour="red",linetype=1,size=1)+geom_hline(yintercept = postp_list[3],colour="red",linetype=2,size=1)+geom_hline(yintercept = postp_list[2],colour="red",linetype=3,size=1)+geom_hline(yintercept = postp_list[1],colour="red",linetype=4,size=1)

pdf("OR_local_lr.pdf",height=5,width=12)
(p1|p2)+ plot_annotation(tag_levels = 'A')
dev.off()
