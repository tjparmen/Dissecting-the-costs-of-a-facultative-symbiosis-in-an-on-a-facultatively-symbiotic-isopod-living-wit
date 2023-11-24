rm(list=ls())

library(lme4)
library(emmeans)
library(ggplot2)
library(RVAideMemoire)
library(afex)
library(DHARMa)

datasetsoil<-read.table("notodef.txt",header=TRUE)  
#edit(datasetsoil)
datasetsoil$soil <- factor(datasetsoil$soil, levels = c("plaster", "sand", "organic"))
levels(datasetsoil$soil)
datasetsoil$obs<-c(1:nrow(datasetsoil))

#aggressive interactions (biting, opening mandibles, chasing)
#####
proportion_aggression<-cbind(datasetsoil$Agr,datasetsoil$Nagr)
fitsoil<-glmer(proportion_aggression~soil*semicochemicals+colony+(1|id),binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),data= datasetsoil)
testDispersion(simulateResiduals(fitsoil, refit=T)) #overdispersion
summary(fitsoil)
fitsoilobs<-glmer(proportion_aggression~soil*semicochemicals+colony+(1|id)+(1|obs),binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),data= datasetsoil)
summary(fitsoilobs)
fitsoiloverdispmix<-mixed(cbind(Agr,Nagr)~soil*semicochemicals+colony+(1|id)+(1|obs),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),method="LRT",family="binomial",data= datasetsoil)

em<-(emmeans(fitsoiloverdispmix, pairwise ~soil|semicochemicals,type="response"))
emsoil<-(emmeans(fitsoiloverdispmix, pairwise ~soil,type="response"))
emcolony<-(emmeans(fitsoiloverdispmix, pairwise ~colony,type="response"))


soil<-as.data.frame(em$emmeans)$'soil'
semichemicals<-as.data.frame(em$emmeans)$'semicochemicals'
prob<-as.data.frame(em$emmeans)$'prob'
asymp.UCL<-as.data.frame(em$emmeans)$'asymp.UCL'
asymp.LCL<-as.data.frame(em$emmeans)$'asymp.LCL'
df<-data.frame(soil,semichemicals,prob,asymp.UCL,asymp.LCL)
df

p<-ggplot(df, aes(soil, prob, color=semichemicals)) + geom_point(position=position_dodge(0.6)) +xlab("")+ geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=0.4,position=position_dodge(0.6))+theme_bw()+theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(axis.text.y=element_text(size=18),axis.text.x = element_text(size=18),axis.title.y=element_text(size=18))+scale_y_continuous(breaks = scales::pretty_breaks(n = 5),name = "Proportion aggressive interactions", limits = c(0, 0.65),expand =c(0,0),labels = scales::number_format(accuracy = 0.01))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
p2<-p+scale_color_manual(labels=c("-","+"),guide=guide_legend(direction="horizontal",title.position = "top"),values=c("#999999", "black"))
p2e<-p2+ theme(legend.position = c(0.75, 0.85))+labs(color = "semiochemicals")
p2e  

#biting
#####
Bproportie<-cbind(datasetsoil$B,20-(datasetsoil$B))
fitsoilB<-glmer(Bproportie~soil*semicochemicals+colony+(1|id),binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),data= datasetsoil)
testOverdispersion(simulateResiduals(fitsoilB, refit=T)) #no overdispersion

datasetsoil$NB<-20-(datasetsoil$B)
fitsoiloverdisp_Bmix<-mixed(cbind(B,NB)~soil*semicochemicals+colony+(1|id),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),method="LRT",family="binomial",data= datasetsoil)
emB<-(emmeans(fitsoiloverdisp_Bmix, pairwise ~semicochemicals|soil,adjust="none",type="response"))
emBsoil<-(emmeans(fitsoiloverdisp_Bmix, pairwise ~soil,adjust="none",type="response"))

soilB<-as.data.frame(emB$emmeans)$'soil'
semichemicalsB<-as.data.frame(emB$emmeans)$'semicochemicals'
probB<-as.data.frame(emB$emmeans)$'prob'
asymp.UCLB<-as.data.frame(emB$emmeans)$'asymp.UCL'
asymp.LCLB<-as.data.frame(emB$emmeans)$'asymp.LCL'
dfB<-data.frame(soilB,semichemicalsB,probB,asymp.UCLB,asymp.LCLB)

p_B<-ggplot(dfB, aes(soilB, probB, color=semichemicalsB)) + geom_point(position=position_dodge(0.6)) +xlab("")+ geom_errorbar(aes(ymin=asymp.UCLB, ymax=asymp.LCLB), width=0.4,position=position_dodge(0.6))+theme_bw()+theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(axis.text.y=element_text(size=18),axis.text.x = element_text(size=18),axis.title.y=element_text(size=18))+scale_y_continuous(breaks = scales::pretty_breaks(n = 5),name = "Proportion biting attempts", limits = c(0, 0.12),expand =c(0,0) ) 
p2_B<-p_B+scale_color_manual(labels=c("-","+"),guide=guide_legend(direction="horizontal",title.position = "top"),values=c("#999999", "black"))
p2_B2<-p2_B+ theme(legend.position = c(0.75, 0.85))+labs(color = "semiochemicals")
p2_B2


#####
#beetle bending abdomen
#

#####
datasetsoil$notbending<-20-datasetsoil$bending
fitsoilbending<-glmer(cbind(bending,notbending)~soil*semicochemicals+colony+(1|id),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),family="binomial",data= datasetsoil)
testOverdispersion(simulateResiduals(fitsoilbending, refit=T)) #overdispersion

fitsoilbendingmix<-mixed(cbind(bending,notbending)~soil*semicochemicals+colony+(1|id)+(1|obs),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),method="LRT",family="binomial",data= datasetsoil)
emmeans(fitsoilbendingmix, "soil", contr = "pairwise", type="response")


