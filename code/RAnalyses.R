library(AUCRF)
library(pROC)
library(dplyr)
library(ggplot2)

pval <- function(p){
  if(p < 0.001){p <- 'p<0.001'}
  else{p <- sprintf('p=%.1g', p)}
  return(p)
}

setwd("/Users/cathylicioussmith/Documents/Bioinf 575/Final Project")

#read in data
meta = read.table('metadata', header=T, sep='\t')
shared_clean = read.table('clean.shared', header=T)
shared_dirty = read.table('dirty.shared', header=T)
#shared files only have 141 observations because samples with low diversity were removed previously

#select columns to use
shared_clean = select(shared_clean,-label,-numOtus)
shared_dirty = select(shared_dirty,-label,-numOtus)

#classify patients in meta data file
meta$dx[(meta$diagnosis_s=='Normal' | meta$diagnosis_s=='High Risk Normal')] = 'normal'
meta$dx[(meta$diagnosis_s=='Adenoma' | meta$diagnosis_s=='adv Adenoma')] = 'adenoma'
meta$dx[meta$diagnosis_s=='Cancer'] = 'cancer'
meta$dx=as.factor(meta$dx)
meta$lesion=1
meta$lesion[meta$dx=='normal'] = 0
meta$lesion = as.factor(meta$lesion)

#create all data dataset
all_data=merge(meta,shared_clean,by.x = 'Sample_Name_s',by.y = 'Group')
all_data=distinct(select(all_data,Sample_Name_s,lesion,dx,fit_result_s,contains("Otu")))
#relevels factor for future figure
all_data$dx=relevel(all_data$dx,"normal")

#merge meta and microbiome data
clean_ade_data=merge(meta,shared_clean,by.x = 'Sample_Name_s',by.y = 'Group')
dirty_ade_data=merge(meta,shared_dirty,by.x = 'Sample_Name_s',by.y = 'Group')
#selects columns for random forest and removes rows duplicating from multiple sequencing runs
clean_ade_data=filter(clean_ade_data,dx!="cancer")
clean_ade_data=distinct(select(clean_ade_data,lesion,contains("Otu")))
dirty_ade_data=filter(dirty_ade_data,dx!="cancer")
dirty_ade_data=distinct(select(dirty_ade_data,lesion,contains("Otu")))

### Adenoma vs Normal models
set.seed(021316)
clean_ade_model = AUCRF(lesion~., data=clean_ade_data, pdel=0.05, ntree=500, ranking='MDA')
clean_ade_probs = predict(clean_ade_model$RFopt, type='prob')[,2]
clean_ade_probs = (clean_ade_probs-0.9*min(clean_ade_probs))/(1-min(clean_ade_probs)) #normalization to spread out probabilities
clean_ade_roc = roc(clean_ade_data$lesion~clean_ade_probs)
set.seed(021316)
dirty_ade_model = AUCRF(lesion~., data=dirty_ade_data, pdel=0.05, ntree=500, ranking='MDA')
dirty_ade_probs = predict(dirty_ade_model$RFopt, type='prob')[,2]
dirty_ade_probs = (dirty_ade_probs-0.9*min(dirty_ade_probs))/(1-min(dirty_ade_probs)) #normalization to spread out probabilities
dirty_ade_roc = roc(dirty_ade_data$lesion~dirty_ade_probs)
ade_roc_test=roc.test(dirty_ade_roc,clean_ade_roc)

### Cancer vs normal models
clean_canc_data=merge(meta,shared_clean,by.x = 'Sample_Name_s',by.y = 'Group')
dirty_canc_data=merge(meta,shared_dirty,by.x = 'Sample_Name_s',by.y = 'Group')
#selects columns for random forest and removes rows duplicating from multiple sequencing runs
clean_canc_data=filter(clean_canc_data,dx!="adenoma")
clean_canc_data=distinct(select(clean_canc_data,lesion,contains("Otu")))
dirty_canc_data=filter(dirty_canc_data,dx!="adenoma")
dirty_canc_data=distinct(select(dirty_canc_data,lesion,contains("Otu")))

set.seed(021316)
clean_canc_model = AUCRF(lesion~., data=clean_canc_data, pdel=0.05, ntree=500, ranking='MDA')
clean_canc_probs = predict(clean_canc_model$RFopt, type='prob')[,2]
clean_canc_roc = roc(clean_canc_data$lesion~clean_canc_probs)
set.seed(021316)
dirty_canc_model = AUCRF(lesion~., data=dirty_canc_data, pdel=0.05, ntree=500, ranking='MDA')
dirty_canc_probs = predict(dirty_canc_model$RFopt, type='prob')[,2]
dirty_canc_roc = roc(dirty_canc_data$lesion~dirty_canc_probs)
canc_roc_test=roc.test(dirty_canc_roc,clean_canc_roc)

### Lesion vs normal models
clean_les_data=merge(meta,shared_clean,by.x = 'Sample_Name_s',by.y = 'Group')
dirty_les_data=merge(meta,shared_dirty,by.x = 'Sample_Name_s',by.y = 'Group')
#selects columns for random forest and removes rows duplicating from multiple sequencing runs
clean_les_data=distinct(select(clean_les_data,lesion,contains("Otu")))
dirty_les_data=distinct(select(dirty_les_data,lesion,contains("Otu")))

set.seed(021316)
clean_les_model = AUCRF(lesion~., data=clean_les_data, pdel=0.05, ntree=500, ranking='MDA')
clean_les_probs = predict(clean_les_model$RFopt, type='prob')[,2]
clean_les_roc = roc(clean_les_data$lesion~clean_les_probs)
set.seed(021316)
dirty_les_model = AUCRF(lesion~., data=dirty_les_data, pdel=0.05, ntree=500, ranking='MDA')
dirty_les_probs = predict(dirty_les_model$RFopt, type='prob')[,2]
dirty_les_roc = roc(dirty_les_data$lesion~dirty_les_probs)
les_roc_test=roc.test(dirty_les_roc,clean_les_roc)

#figures
plot(clean_les_model, cex.lab=1.25, cex.axis=1.25)
title("Prediction Accuracy for Chimera- Lesion Classification \nAcross Number of Variables Selected", cex.main=1.5)

#compares dirty and clean for cancer
F2D <- ggroc(list(clean_canc_roc,dirty_canc_roc),size=2.5) +
  geom_abline(intercept = 1, slope = 1, linetype="dashed", size=1.5) +
  labs(title="Classification of Cancer Patients", 
       y = "Sensitivity", x="Specificity", subtitle=paste("Difference in AUC: z=",paste(abs(round(canc_roc_test$statistic,3)),pval(canc_roc_test$p.value))))+
  scale_color_manual(name="Dataset", values = c("red","blue"),
                     labels=c(sprintf('Chimera-: AUC: %.3g', clean_canc_roc$auc), 
                              sprintf('Chimera+: AUC: %.3g', dirty_canc_roc$auc)))+
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 20),
        axis.title.x =  element_text(margin = margin(t=15,unit = "pt"),face = "bold"), 
        legend.background = element_rect(linetype = 'solid',color = 'black'),
        plot.title = element_text(hjust = .5,size = 28, face = "bold"), legend.text = element_text(size=24,margin = margin(b=10,unit="pt")),
        legend.position=c(.7,.2),legend.title = element_text(size=24, face = "bold"), strip.text = element_text(size = 24, face = "bold"),
        legend.title.align = .5, legend.key.size = unit(40,"pt"), axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"), plot.subtitle = element_text(hjust = .5,size = 20),
        aspect.ratio = 1) +
  coord_fixed()

#compares dirty and clean for adenomas
F2A <- ggroc(list(clean_ade_roc,dirty_ade_roc),size=2.5) +
  geom_abline(intercept = 1, slope = 1, linetype="dashed", size=1.5) +
  labs(title="Classification of Adenoma Patients", 
       y = "Sensitivity", x="Specificity", subtitle=paste("Difference in AUC: z=",paste(abs(round(ade_roc_test$statistic,3)),pval(ade_roc_test$p.value))))+
  scale_color_manual(name="Dataset", values = c("red","blue"),
                     labels=c(sprintf('Chimera-: AUC: %.3g', clean_ade_roc$auc), 
                              sprintf('Chimera+: AUC: %.3g', dirty_ade_roc$auc)))+
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 20),
        axis.title.x =  element_text(margin = margin(t=15,unit = "pt"),face = "bold"), 
        legend.background = element_rect(linetype = 'solid',color = 'black'),
        plot.title = element_text(hjust = .5,size = 28, face = "bold"), legend.text = element_text(size=24,margin = margin(b=10,unit="pt")),
        legend.position=c(.7,.2),legend.title = element_text(size=24, face = "bold"), strip.text = element_text(size = 24, face = "bold"),
        legend.title.align = .5, legend.key.size = unit(40,"pt"), axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"), plot.subtitle = element_text(hjust = .5,size = 20),
        aspect.ratio = 1) +
  coord_fixed()

#compares dirty and clean for lesion
F2G <- ggroc(list(clean_les_roc,dirty_les_roc),size=2.5) +
  geom_abline(intercept = 1, slope = 1, linetype="dashed", size=1.5) +
  labs(title="Classification of Lesion Patients", 
       y = "Sensitivity", x="Specificity", subtitle=paste("Difference in AUC: z=",paste(abs(round(les_roc_test$statistic,3)),pval(les_roc_test$p.value))))+
  scale_color_manual(name="Dataset", values = c("red","blue"),
                     labels=c(sprintf('Chimera-: AUC: %.3g', clean_les_roc$auc), 
                              sprintf('Chimera+: AUC: %.3g', dirty_les_roc$auc)))+
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 20),
        axis.title.x =  element_text(margin = margin(t=15,unit = "pt"),face = "bold"), 
        legend.background = element_rect(linetype = 'solid',color = 'black'),
        plot.title = element_text(hjust = .5,size = 28, face = "bold"), legend.text = element_text(size=24,margin = margin(b=10,unit="pt")),
        legend.position=c(.7,.2),legend.title = element_text(size=24, face = "bold"), strip.text = element_text(size = 24, face = "bold"),
        legend.title.align = .5, legend.key.size = unit(40,"pt"), axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"), plot.subtitle = element_text(hjust = .5,size = 20),
        aspect.ratio = 1) +
  coord_fixed()

#calculate Youden's J statistic for optimal cut off
canc_roc=roc(clean_canc_data$lesion~clean_canc_probs)
canc_cutoff <- coords(roc=canc_roc, x='best', best.method='y', ret=c('threshold'))
ade_roc=roc(clean_ade_data$lesion~clean_ade_probs)
ade_cutoff <- coords(roc=ade_roc, x='best', best.method='y', ret=c('threshold'))
les_roc=roc(clean_les_data$lesion~clean_les_probs)
les_cutoff <- coords(roc=les_roc, x='best', best.method='y', ret=c('threshold'))

#classify cancer model
clean_canc_TP=clean_canc_probs[clean_canc_data$lesion==1 & clean_canc_probs>=canc_cutoff]
clean_canc_sens=length(clean_canc_TP)/length(clean_canc_probs[clean_canc_data$lesion==1])
clean_canc_TN=clean_canc_probs[clean_canc_data$lesion==0 & clean_canc_probs<canc_cutoff]
clean_canc_spec=length(clean_canc_TN)/length(clean_canc_probs[clean_canc_data$lesion==0])
clean_canc_correct=(length(clean_canc_TP)+length(clean_canc_TN))/length(clean_canc_data$lesion)
dirty_canc_TP=dirty_canc_probs[dirty_canc_data$lesion==1 & dirty_canc_probs>=canc_cutoff]
dirty_canc_sens=length(dirty_canc_TP)/length(dirty_canc_probs[dirty_canc_data$lesion==1])
dirty_canc_TN=dirty_canc_probs[dirty_canc_data$lesion==0 & dirty_canc_probs<canc_cutoff]
dirty_canc_spec=length(dirty_canc_TN)/length(dirty_canc_probs[dirty_canc_data$lesion==0])
dirty_canc_correct=(length(dirty_canc_TP)+length(dirty_canc_TN))/length(dirty_canc_data$lesion)

#McNemar to see which if there's a difference in who the test has better sensitivity (diseased patients only) & specificity (normal patients only)
canc_probs=as.data.frame(cbind(clean_canc_data$lesion,clean_canc_probs,dirty_canc_probs))
canc_mcnemar_b_sens=length(which(canc_probs$V1==2 & canc_probs$clean_canc_probs>=canc_cutoff & canc_probs$dirty_canc_probs<canc_cutoff))
canc_mcnemar_c_sens=length(which(canc_probs$V1==2 & canc_probs$clean_canc_probs<canc_cutoff & canc_probs$dirty_canc_probs>=canc_cutoff))
canc_binom_sens=binom.test(canc_mcnemar_b_sens,canc_mcnemar_b_sens+canc_mcnemar_c_sens,.5)
canc_mcnemar_b_spec=length(which(canc_probs$V1==1 & canc_probs$clean_canc_probs<canc_cutoff & canc_probs$dirty_canc_probs>=canc_cutoff))
canc_mcnemar_c_spec=length(which(canc_probs$V1==1 & canc_probs$clean_canc_probs>=canc_cutoff & canc_probs$dirty_canc_probs<canc_cutoff))
canc_binom_spec=binom.test(canc_mcnemar_b_spec,canc_mcnemar_b_spec+canc_mcnemar_c_spec,.5)
#didn't use mcnemar bc low cell counts

#classify ade model
clean_ade_TP=clean_ade_probs[clean_ade_data$lesion==1 & clean_ade_probs>=ade_cutoff]
clean_ade_sens=length(clean_ade_TP)/length(clean_ade_probs[clean_ade_data$lesion==1])
clean_ade_TN=clean_ade_probs[clean_ade_data$lesion==0 & clean_ade_probs<ade_cutoff]
clean_ade_spec=length(clean_ade_TN)/length(clean_ade_probs[clean_ade_data$lesion==0])
clean_ade_correct=(length(clean_ade_TP)+length(clean_ade_TN))/length(clean_ade_data$lesion)
dirty_ade_TP=dirty_ade_probs[dirty_ade_data$lesion==1 & dirty_ade_probs>=ade_cutoff]
dirty_ade_sens=length(dirty_ade_TP)/length(dirty_ade_probs[dirty_ade_data$lesion==1])
dirty_ade_TN=dirty_ade_probs[dirty_ade_data$lesion==0 & dirty_ade_probs<ade_cutoff]
dirty_ade_spec=length(dirty_ade_TN)/length(dirty_ade_probs[dirty_ade_data$lesion==0])
dirty_ade_correct=(length(dirty_ade_TP)+length(dirty_ade_TN))/length(dirty_ade_data$lesion)

#McNemar to see which if there's a difference in who the test has better sensitivity (diseased patients only) & specificity (normal patients only)
ade_probs=as.data.frame(cbind(clean_ade_data$lesion,clean_ade_probs,dirty_ade_probs))
ade_mcnemar_b_sens=length(which(ade_probs$V1==2 & ade_probs$clean_ade_probs>=ade_cutoff & ade_probs$dirty_ade_probs<ade_cutoff))
ade_mcnemar_c_sens=length(which(ade_probs$V1==2 & ade_probs$clean_ade_probs<ade_cutoff & ade_probs$dirty_ade_probs>=ade_cutoff))
ade_binom_sens=binom.test(ade_mcnemar_b_sens,ade_mcnemar_b_sens+ade_mcnemar_c_sens,.5)
ade_mcnemar_b_spec=length(which(ade_probs$V1==1 & ade_probs$clean_ade_probs<ade_cutoff & ade_probs$dirty_ade_probs>=ade_cutoff))
ade_mcnemar_c_spec=length(which(ade_probs$V1==1 & ade_probs$clean_ade_probs>=ade_cutoff & ade_probs$dirty_ade_probs<ade_cutoff))
ade_binom_spec=binom.test(ade_mcnemar_b_spec,ade_mcnemar_b_spec+ade_mcnemar_c_spec,.5)
#didn't use mcnemar bc low cell counts

#classify les model
clean_les_TP=clean_les_probs[clean_les_data$lesion==1 & clean_les_probs>=les_cutoff]
clean_les_sens=length(clean_les_TP)/length(clean_les_probs[clean_les_data$lesion==1])
clean_les_TN=clean_les_probs[clean_les_data$lesion==0 & clean_les_probs<les_cutoff]
clean_les_spec=length(clean_les_TN)/length(clean_les_probs[clean_les_data$lesion==0])
clean_les_correct=(length(clean_les_TP)+length(clean_les_TN))/length(clean_les_data$lesion)
dirty_les_TP=dirty_les_probs[dirty_les_data$lesion==1 & dirty_les_probs>=les_cutoff]
dirty_les_sens=length(dirty_les_TP)/length(dirty_les_probs[dirty_les_data$lesion==1])
dirty_les_TN=dirty_les_probs[dirty_les_data$lesion==0 & dirty_les_probs<les_cutoff]
dirty_les_spec=length(dirty_les_TN)/length(dirty_les_probs[dirty_les_data$lesion==0])
dirty_les_correct=(length(dirty_les_TP)+length(dirty_les_TN))/length(dirty_les_data$lesion)

#McNemar to see which if there's a difference in who the test has better sensitivity (diseased patients only) & specificity (normal patients only)
les_probs=as.data.frame(cbind(clean_les_data$lesion,clean_les_probs,dirty_les_probs,all_data$dx))
les_mcnemar_b_sens=length(which(les_probs$V1==2 & les_probs$clean_les_probs>=les_cutoff & les_probs$dirty_les_probs<les_cutoff))
les_mcnemar_c_sens=length(which(les_probs$V1==2 & les_probs$clean_les_probs<les_cutoff & les_probs$dirty_les_probs>=les_cutoff))
les_binom_sens=binom.test(les_mcnemar_b_sens,les_mcnemar_b_sens+les_mcnemar_c_sens,.5)
les_mcnemar_b_spec=length(which(les_probs$V1==1 & les_probs$clean_les_probs<les_cutoff & les_probs$dirty_les_probs>=les_cutoff))
les_mcnemar_c_spec=length(which(les_probs$V1==1 & les_probs$clean_les_probs>=les_cutoff & les_probs$dirty_les_probs<les_cutoff))
les_binom_spec=binom.test(les_mcnemar_b_spec,les_mcnemar_b_spec+les_mcnemar_c_spec,.5)
#didn't use mcnemar bc low cell counts

#create Youden's J figures
Fig2B=ggplot(ade_probs, aes(x=as.factor(V1), y=clean_ade_probs, color=as.factor(V1))) +
  geom_point(position=position_jitter(),size=5)+
  geom_hline(yintercept =ade_cutoff, linetype="dashed", size=1.5) +
  labs(title="Chimera- Classification of Adenoma \n at Optimal Probability Cut Off", 
       subtitle=paste(paste(paste(sprintf('Sensitivity: %.3g', 100*clean_ade_sens), '%; ', sep=""),
                      paste(sprintf('Specificity: %.3g', 100*clean_ade_spec), '%; ', sep="")),
                      paste(sprintf('Accuracy: %.3g', 100*clean_ade_correct), '%', sep="")),
       y = "Probability of Adenoma")+
  scale_x_discrete(name="Dataset", labels=c("Normal","Adenoma")) +
  scale_color_manual(values = c("blue","red")) +
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 24, face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none", plot.title = element_text(hjust = .5,size = 27, face="bold"), 
        strip.text = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"),
        plot.subtitle = element_text(hjust = .5,size = 18))

Fig2C=ggplot(ade_probs, aes(x=as.factor(V1), y=dirty_ade_probs, color=as.factor(V1))) +
  geom_point(position=position_jitter(),size=5)+
  geom_hline(yintercept =ade_cutoff, linetype="dashed", size=1.5) +
  labs(title="Chimera+ Classification of Adenoma \n at Optimal Probability Cut Off", 
       subtitle=paste(paste(paste(sprintf('Sensitivity: %.3g', 100*dirty_ade_sens), '%; ', sep=""),
                            paste(sprintf('Specificity: %.3g', 100*dirty_ade_spec), '%; ', sep="")),
                      paste(sprintf('Accuracy: %.3g', 100*dirty_ade_correct), '%', sep="")),
       y = "Probability of Adenoma")+
  scale_x_discrete(name="Dataset", labels=c("Normal","Adenoma")) +
  scale_color_manual(values = c("blue","red")) +
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 24, face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none", plot.title = element_text(hjust = .5,size = 27, face="bold"), 
        strip.text = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"),
        plot.subtitle = element_text(hjust = .5,size = 18))

Fig2E=ggplot(canc_probs, aes(x=as.factor(V1), y=clean_canc_probs, color=as.factor(V1))) +
  geom_point(position=position_jitter(),size=5)+
  geom_hline(yintercept =canc_cutoff, linetype="dashed", size=1.5) +
  labs(title="Chimera- Classification of Cancer \n at Optimal Probability Cut Off", 
       subtitle=paste(paste(paste(sprintf('Sensitivity: %.3g', 100*clean_canc_sens), '%; ', sep=""),
                            paste(sprintf('Specificity: %.3g', 100*clean_canc_spec), '%; ', sep="")),
                      paste(sprintf('Accuracy: %.3g', 100*clean_canc_correct), '%', sep="")),
       y = "Probability of Cancer")+
  scale_x_discrete(name="Dataset", labels=c("Normal","Cancer")) +
  scale_color_manual(values = c("blue","red")) +
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 24, face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none", plot.title = element_text(hjust = .5,size = 27, face="bold"), 
        strip.text = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"),
        plot.subtitle = element_text(hjust = .5,size = 18))

Fig2F=ggplot(canc_probs, aes(x=as.factor(V1), y=dirty_canc_probs, color=as.factor(V1))) +
  geom_point(position=position_jitter(),size=5)+
  geom_hline(yintercept =canc_cutoff, linetype="dashed", size=1.5) +
  labs(title="Chimera+ Classification of Cancer \n at Optimal Probability Cut Off", 
       subtitle=paste(paste(paste(sprintf('Sensitivity: %.3g', 100*dirty_canc_sens), '%; ', sep=""),
                            paste(sprintf('Specificity: %.3g', 100*dirty_canc_spec), '%; ', sep="")),
                      paste(sprintf('Accuracy: %.3g', 100*dirty_canc_correct), '%', sep="")),
       y = "Probability of Cancer")+
  scale_x_discrete(name="Dataset", labels=c("Normal","Cancer")) +
  scale_color_manual(values = c("blue","red")) +
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 24, face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none", plot.title = element_text(hjust = .5,size = 27, face="bold"), 
        strip.text = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"),
        plot.subtitle = element_text(hjust = .5,size = 18))

Fig2H=ggplot(les_probs, aes(x=as.factor(V4), y=clean_les_probs, color=as.factor(V4))) +
  geom_point(position=position_jitter(),size=5)+
  geom_hline(yintercept =les_cutoff, linetype="dashed", size=1.5) +
  labs(title="Chimera- Classification of Lesions \n at Optimal Probability Cut Off", 
       subtitle=paste(paste(paste(sprintf('Sensitivity: %.3g', 100*clean_les_sens), '%; ', sep=""),
                            paste(sprintf('Specificity: %.3g', 100*clean_les_spec), '%; ', sep="")),
                      paste(sprintf('Accuracy: %.3g', 100*clean_les_correct), '%', sep="")),
       y = "Probability of Lesion")+
  scale_x_discrete(name="Dataset", labels=c("Normal","Adenoma","Cancer")) +
  scale_color_manual(values = c("blue","red","gold")) +
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 24, face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none", plot.title = element_text(hjust = .5,size = 27, face="bold"), 
        strip.text = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"),
        plot.subtitle = element_text(hjust = .5,size = 18))

Fig2I=ggplot(les_probs, aes(x=as.factor(V4), y=dirty_les_probs, color=as.factor(V4))) +
  geom_point(position=position_jitter(),size=5)+
  geom_hline(yintercept =les_cutoff, linetype="dashed", size=1.5) +
  labs(title="Chimera+ Classification of Lesion \n at Optimal Probability Cut Off", 
       subtitle=paste(paste(paste(sprintf('Sensitivity: %.3g', 100*dirty_les_sens), '%; ', sep=""),
                            paste(sprintf('Specificity: %.3g', 100*dirty_les_spec), '%; ', sep="")),
                      paste(sprintf('Accuracy: %.3g', 100*dirty_les_correct), '%', sep="")),
       y = "Probability of Lesion")+
  scale_x_discrete(name="Dataset",labels=c("Normal","Adenoma", "Cancer")) +
  scale_color_manual(values = c("blue","red","gold")) +
  theme_classic() +
  theme(text=element_text(size=24),axis.text.x=element_text(size = 24, face = "bold"),
        axis.title.x =  element_blank(),
        legend.position = "none", plot.title = element_text(hjust = .5,size = 27, face="bold"), 
        strip.text = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(size = 26),
        axis.title.y = element_text(margin = margin(r=15,unit = "pt"),face = "bold"),
        plot.subtitle = element_text(hjust = .5,size = 18))

#les_mcnemar=(les_mcnemar_b-les_mcnemar_c)^2/(les_mcnemar_b+les_mcnemar_c)