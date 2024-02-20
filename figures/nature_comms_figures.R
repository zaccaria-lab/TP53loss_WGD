#Figure 1f, Supplementary Figure 3a, 

library(tidyverse)
library(ggpubr)
library(cowplot)
library(dplyr)
library(arsenal)
library(forcats)


defaultW <- getOption("warn") 
options(warn = -1) 


#Data_tables and data format
TP53_lesion_specifc_first<-read.csv("~/figures/TP53_fig1f.csv")
TP53_lesion_specifc_first$p53.MDM2<-fct_rev(as.factor(TP53_lesion_specifc_first$p53.MDM2))

TP53_imaging_first<-read.csv("~/figures/TP53_fig1f_supp.csv")
TP53_imaging_first$p53.MDM2<-fct_rev(as.factor(TP53_imaging_first$p53.MDM2))

TP53_lesion_specifc_max <-read.csv("~/figures/TP53_supp3a.csv")
TP53_lesion_specifc_max$p53.MDM2<-fct_rev(as.factor(TP53_lesion_specifc_max$p53.MDM2))

TP53_imaging_max<-read.csv("~/figures/TP53_supp3a_supp.csv")
TP53_imaging_max$p53.MDM2<-fct_rev(as.factor(TP53_imaging_max$p53.MDM2))

###############
##
#Figure 1f
##
###############


#Plot fold change from baseline to first scan of individual lesions
tmp<-ggplot(TP53_lesion_specifc_first, aes(x = reorder(patient_id, foldchange_mean),Percent.change.at.first.scan, fill=p53.MDM2))

plot_first<-tmp+stat_summary(fun.data=median_mad, geom= "errorbar", width =1, color="grey48")+
  stat_summary(fun = median,
               size = 0.1, colour ="grey48") +
  geom_dotplot(binaxis = "y", stackdir = 'center', stackratio = 0.1, dotsize=1.4,position = position_jitter(height = 2, width = 0) ,binwidth=5,stroke = 0.5)+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(
    axis.title.y = element_text( size=20),axis.text = element_text(size = 15) )+
  scale_fill_manual(values=c("royalblue", "red"))+
  theme(legend.position = "none", axis.title.x=element_blank())+
  theme(plot.margin = margin(0, 0.1, 0, 0.1, "cm"))+
  theme(axis.title.y=element_blank() )+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
  geom_hline(yintercept=c(10), linetype="dotted", colour = "gray18")+
  geom_hline(yintercept=c(-30), linetype="dotted", colour="gray18")+
  theme(
    strip.background = element_blank()
  )+scale_x_discrete(expand = c(.05, .05))

ordering<-as.vector(reorder(TP53_lesion_specifc_first$patient_id, TP53_lesion_specifc_first$foldchange_mean)%>%unique())

TP53_imaging_first$patient_id<-factor(TP53_imaging_first$patient_id, levels = ordering)

#annotate patients who develop new lesions 
plot_first_new<-ggplot(TP53_imaging_first, aes(patient_id, y=1, fill = new_lesion)) +
  geom_tile()+
  scale_fill_manual(values=c("white", "#E69F00"))+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), line = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))+
  scale_x_discrete(expand = c(.05, .05))

#annotate patients with progressing and responding lesions
plot_first_het<-ggplot(TP53_imaging_first, aes(patient_id, y=1, fill = het_response_first)) +
  geom_tile()+
  scale_fill_manual(values=c("#56B4E9", "white"))+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), line = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),plot.margin=unit(c(-0.1,0.1,0,0.1), "cm")
  )+
  scale_x_discrete(expand = c(.05, .05))

#annotate patients with heterogeneous response
plot_first_overall<-ggplot(TP53_imaging_first, aes(patient_id, y=1, fill = response_final_first)) +
  geom_tile()+
  scale_fill_manual(values=c("#009E73", "white"))+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), line = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),plot.margin=unit(c(-0.1,0.1,0.1,0.1), "cm")
  )+
  scale_x_discrete(expand = c(.05, .05))

figure1f<-ggarrange(plot_first, plot_first_new, plot_first_het,plot_first_overall,
                  ncol = 1, nrow = 4, heights = c(10, 0.5,0.3, 0.4) )

figure1f



###############
##
## Sup Figure 3a
##
###############


#Plot fold change from baseline to best response of individual lesions
tmp<-ggplot(TP53_lesion_specifc_max, aes(x = reorder(patient_id, foldchange_mean),Percent.change.at.MAX, fill=p53.MDM2))

plot_max<-tmp+stat_summary(fun.data=median_mad, geom= "errorbar", width =1, color="grey48")+
  stat_summary(fun = median,
               size = 0.1, colour ="grey48") +
  geom_dotplot(binaxis = "y", stackdir = 'center', stackratio = 0.1, dotsize=1.4,position = position_jitter(height = 2, width = 0) ,binwidth=5,stroke = 0.5)+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(
    axis.title.y = element_text( size=20),axis.text = element_text(size = 15) )+
  scale_fill_manual(values=c("royalblue", "red"))+
  theme(legend.position = "none", axis.title.x=element_blank())+
  theme(plot.margin = margin(0, 0.1, 0, 0.1, "cm"))+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank() )+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
  geom_hline(yintercept=c(10), linetype="dotted", colour = "gray18")+
  geom_hline(yintercept=c(-30), linetype="dotted", colour="gray18")+
  theme(
    strip.background = element_blank()
  )+
  scale_x_discrete(expand = c(.05, .05))

ordering<-as.vector(reorder(TP53_lesion_specifc_max$patient_id, TP53_lesion_specifc_max$foldchange_mean)%>%unique())

TP53_imaging_max$patient_id<-factor(TP53_imaging_max$patient_id, levels = ordering)

#annotate patients who develop new lesions 
plot_max_new<-ggplot(TP53_imaging_max, aes(x=reorder(patient_id,ordering), y=1, fill = new_lesion)) +
  geom_tile()+
  scale_fill_manual(values=c("white", "#E69F00"))+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), line = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),plot.margin=unit(c(0.1,0.1,0,0.1), "cm"))+
  scale_x_discrete(expand = c(.05, .05))

#annotate patients with progressing and responding lesions
plot_max_het<-ggplot(TP53_imaging_max, aes(patient_id, y=1, fill = het_response_max)) +
  geom_tile()+
  scale_fill_manual(values=c("#56B4E9", "white"))+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), line = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),plot.margin=unit(c(-0.1,0.1,0,0.1), "cm")
  )+
  scale_x_discrete(expand = c(.05, .05))

#annotate patients with heterogeneous response
plot_max_overall<-ggplot(TP53_imaging_max, aes(patient_id, y=1, fill = response_final_max)) +
  geom_tile()+
  scale_fill_manual(values=c("#009E73", "white"))+
  facet_grid( ~ p53.MDM2, scales = "free", space = "free")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), line = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),plot.margin=unit(c(-0.1,0.1,0.1,0.1), "cm")
  )+
  scale_x_discrete(expand = c(.05, .05))


Supplementary_figure_3a<-ggarrange(plot_max, plot_max_new, plot_max_het,plot_max_overall,
                   ncol = 1, nrow = 4, heights = c(10, 0.5,0.3, 0.4) )

Supplementary_figure_3a


options(warn = defaultW)
