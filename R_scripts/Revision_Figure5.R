
####
#### Final Figure 5 plots
####

### Also includes code for Supplemental figure S6

###Model for assortative mating for selfs only
Outself_fruits_combine_lmready_selfs<-Outself_fruits_combine_lmready %>% filter(selfcall=="self") #%>% filter(offspring_count>1)
summary(lm(fruitmom_sway~mom_nas_prop,data=Outself_fruits_combine_lmready_selfs))
summary(lm(fruitmom_sway~mom_nas_prop*Year,data=Outself_fruits_combine_lmready_selfs))
summary(lm(fruitmom_sway~mom_nas_prop*Year*stream,data=Outself_fruits_combine_lmready_selfs))

offspringxmaternal_outself_fruits_threeyears_selfs<-Outself_fruits_combine_lmready_selfs %>%
  ggplot(aes(x=mom_nas_prop,y=fruitmom_sway,color=Year,fill=Year,shape=selfcall,linetype=stream)) + 
  geom_smooth(method="lm",alpha=0.3) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  scale_fill_manual(values=colors_yr[c(2,3,4)]) + 
  scale_shape_manual(name="Mating history",breaks=c("self"),values=c(1),labels=c("selfed")) + 
  scale_linetype_manual(name="Stream",breaks=c("CAC_stream1","CAC_stream2"),values=c("solid","dotted"),labels=c("CAC_S1","CAC_S2")) + 
  guides(shape=guide_legend(override.aes=list(color="black"))) + 
  guides(linetype=guide_legend(override.aes=list(color="black"))) + 
  xlab("Maternal HI") + ylab("Offspring mean HI - Maternal HI") + 
  xlim(0,0.8) + ylim(-0.5,0.5) + 
  geom_hline(yintercept=0,linetype="dotted",color="black")
offspringxmaternal_outself_fruits_threeyears_selfs
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/offspringxmaternal_outself_fruits_threeyears_selfs.png",offspringxmaternal_outself_fruits_threeyears_selfs,device = 'png',height = 4,width=5,dpi=600,units='in')

###Model for assortative mating for oucrosses only
Outself_fruits_combine_lmready_outs<-Outself_fruits_combine_lmready %>% filter(selfcall=="out")
summary(lm(fruitmom_sway~mom_nas_prop,data=Outself_fruits_combine_lmready_outs))
summary(lm(fruitmom_sway~mom_nas_prop*Year,data=Outself_fruits_combine_lmready_outs))
summary(lm(fruitmom_sway~mom_nas_prop*Year*stream,data=Outself_fruits_combine_lmready_outs))

offspringxmaternal_outself_fruits_threeyears_outs<-Outself_fruits_combine_lmready_outs %>%
  ggplot(aes(x=mom_nas_prop,y=fruitmom_sway,color=Year,fill=Year,shape=selfcall,linetype=stream)) + 
  geom_smooth(method="lm",alpha=0.3) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  scale_fill_manual(values=colors_yr[c(2,3,4)]) + 
  scale_shape_manual(name="Mating history",breaks=c("out"),values=c(16),labels=c("outcrossed")) + 
  #scale_linetype_discrete(name="Stream") + 
  scale_linetype_manual(name="Stream",breaks=c("CAC_stream1","CAC_stream2"),values=c("solid","dotted"),labels=c("CAC_S1","CAC_S2")) + 
  guides(shape=guide_legend(override.aes=list(color="black"))) + 
  guides(linetype=guide_legend(override.aes=list(color="black"))) + 
  xlab("Maternal HI") + ylab("Offspring mean HI - Maternal HI") + 
  xlim(0,0.8) + ylim(-0.5,0.5) + 
  geom_hline(yintercept=0,linetype="dotted",color="black")
offspringxmaternal_outself_fruits_threeyears_outs
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/SEP_versions/offspringxmaternal_outself_fruits_threeyears_outs.png",offspringxmaternal_outself_fruits_threeyears_outs,device = 'png',height = 4,width=5,dpi=600,units='in')

###Both together
Rcomplete_Fig5<-plot_grid(offspringxmaternal_outself_fruits_threeyears_selfs,offspringxmaternal_outself_fruits_threeyears_outs,nrow=1)
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/Rcomplete_Fig5.pdf",Rcomplete_Fig5,device = 'pdf',height = 4,width=8,dpi=600,units='in')


###Supp Fig 6
###selfing rate with confidence intervals
markedfruits_strict_modelled<-markedfruits_strict %>% mutate(stream=factor(stream,levels=c("CAC_stream1","CAC_stream2")),
                                                             Year=factor(Year,levels=c(2019,2022))) %>%
  filter(mom_isSOOK == "no") %>%
  mutate(mom_nas_prop=ifelse(mom_nas_prop==1,0.999,mom_nas_prop))

mygam = gam(cbind(count_self,count_out)~nas_prop*Year*stream,data=Outself_moms_combine_glmready,family=binomial(link="logit"))
range.2019.s1<-filter(Outself_moms_combine_glmready,Year=="2019",stream=="CAC_stream1")$nas_prop
range.2021.s1<-filter(Outself_moms_combine_glmready,Year=="2021",stream=="CAC_stream1")$nas_prop
range.2021.s2<-filter(Outself_moms_combine_glmready,Year=="2021",stream=="CAC_stream2")$nas_prop
range.2022.s1<-filter(Outself_moms_combine_glmready,Year=="2022",stream=="CAC_stream1")$nas_prop
range.2022.s2<-filter(Outself_moms_combine_glmready,Year=="2022",stream=="CAC_stream2")$nas_prop
new.x.2019.stream1 <- data.frame(nas_prop = seq(min(range.2019.s1), max(range.2019.s1), length.out = 1000),Year="2019",stream="CAC_stream1")
new.x.2021.stream1 <- data.frame(nas_prop = seq(min(range.2021.s1), max(range.2021.s1), length.out = 1000),Year="2021",stream="CAC_stream1")
new.x.2021.stream2 <- data.frame(nas_prop = seq(min(range.2021.s2), max(range.2021.s2), length.out = 1000),Year="2021",stream="CAC_stream2")
new.x.2022.stream1 <- data.frame(nas_prop = seq(min(range.2022.s1), max(range.2022.s1), length.out = 1000),Year="2022",stream="CAC_stream1")
new.x.2022.stream2 <- data.frame(nas_prop = seq(min(range.2022.s2), max(range.2022.s2), length.out = 1000),Year="2022",stream="CAC_stream2")
new.x<-rbind(new.x.2019.stream1,new.x.2021.stream1,new.x.2021.stream2,new.x.2022.stream1,new.x.2022.stream2)
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, fraction_selfed = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = fraction_selfed - 1.96 * SE, upr = fraction_selfed + 1.96 * SE) # calculating the 95% confidence interval

Selfingrate_model_scatter_se<-Outself_moms_combine_glmready %>% 
  ggplot(aes(x=nas_prop,y=fraction_selfed,color=Year,shape=stream,linetype=stream)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values=colors_yr[c(2,3,4)]) + 
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity') + 
  scale_linetype_discrete(name="Stream",labels=c("CAC_S1","CAC_S2"),guide=guide_legend(override.aes=list(color="black"))) + 
  scale_shape_discrete(name="Stream",labels=c("CAC_S1","CAC_S2"),guide=guide_legend(override.aes=list(color="black"))) + 
  ylab("Proportion selfed offspring") + xlab("Maternal HI")
Selfingrate_model_scatter_se
ggsave("~/OneDrive - University of Georgia/General - Sweigart Lab/FarnitanoMatt/Ch2_CAC_popgen/Figures/DEC_revision/Selfingrate_model_scatter_se.pdf",Selfingrate_model_scatter_se,device = 'pdf',height = 4,width=5,dpi=600,units='in')


