### GAM MODELS

#Plotwise HI x median flower date
nrow(HI_plot_quartiles)
plot_quartiles_corr_null<-gam(median_HI~1,data=HI_plot_quartiles,family=betar(link="logit"))
plot_quartiles_corr_date<-gam(median_HI~median_flowerdate,data=HI_plot_quartiles,family=betar(link="logit"))
plot_quartiles_corr_stream<-gam(median_HI~median_flowerdate*stream,data=HI_plot_quartiles,family=betar(link="logit"))
plot_quartiles_corr_Year<-gam(median_HI~median_flowerdate*stream*Year,data=HI_plot_quartiles,family=betar(link="logit"))

lrtest(plot_quartiles_corr_null,plot_quartiles_corr_date,plot_quartiles_corr_stream,plot_quartiles_corr_Year)
summary(plot_quartiles_corr_date)

#Individual HI x flower date
nrow(markedfruits_strict_modelled)
fruit_timing_corr_null<-gam(mom_nas_prop~1,data=markedfruits_strict_modelled,family=betar(link="logit"))
fruit_timing_corr_date<-gam(mom_nas_prop~Days_start_Apr1,data=markedfruits_strict_modelled,family=betar(link="logit"))
fruit_timing_corr_stream<-gam(mom_nas_prop~Days_start_Apr1*stream,data=markedfruits_strict_modelled,family=betar(link="logit"))
fruit_timing_corr_Year<-gam(mom_nas_prop~Days_start_Apr1*stream*Year,data=markedfruits_strict_modelled,family=betar(link="logit"))

lrtest(fruit_timing_corr_null,fruit_timing_corr_date,fruit_timing_corr_stream,fruit_timing_corr_Year)
summary(fruit_timing_corr_Year)

#Selfing rate (with M. nasutus)
nrow(Outself_moms_combine_glmready_includenas)
selfing_moms_withnas_glm_null<-gam(cbind(count_self,count_out)~1,family = binomial(link = "logit"),data=Outself_moms_combine_glmready_includenas)
selfing_moms_withnas_glm_HI<-gam(cbind(count_self,count_out)~nas_prop,family = binomial(link = "logit"),data=Outself_moms_combine_glmready_includenas)
selfing_moms_withnas_glm_stream<-gam(cbind(count_self,count_out)~nas_prop*stream,family = binomial(link = "logit"),data=Outself_moms_combine_glmready_includenas)
selfing_moms_withnas_glm_Year<-gam(cbind(count_self,count_out)~nas_prop*stream*Year,family = binomial(link = "logit"),data=Outself_moms_combine_glmready_includenas)

lrtest(selfing_moms_withnas_glm_null,selfing_moms_withnas_glm_HI,selfing_moms_withnas_glm_stream,selfing_moms_withnas_glm_Year)
summary(selfing_moms_withnas_glm_Year)

#Selfing rate (no M. nasutus)
nrow(Outself_moms_combine_glmready)
selfing_moms_nonas_glm_null<-gam(cbind(count_self,count_out)~1,family = binomial(link = "logit"),data=Outself_moms_combine_glmready)
selfing_moms_nonas_glm_HI<-gam(cbind(count_self,count_out)~nas_prop,family = binomial(link = "logit"),data=Outself_moms_combine_glmready)
selfing_moms_nonas_glm_stream<-gam(cbind(count_self,count_out)~nas_prop*stream,family = binomial(link = "logit"),data=Outself_moms_combine_glmready)
selfing_moms_nonas_glm_Year<-gam(cbind(count_self,count_out)~nas_prop*Year,family = binomial(link = "logit"),data=Outself_moms_combine_glmready)

lrtest(selfing_moms_nonas_glm_null,selfing_moms_nonas_glm_HI,selfing_moms_nonas_glm_stream,selfing_moms_nonas_glm_Year)
lrtest(selfing_moms_nonas_glm_null,selfing_moms_nonas_glm_HI,selfing_moms_nonas_glm_Year)
summary(selfing_moms_nonas_glm_Year)

#Selfing rate (exclude HI>0.5)
nrow(filter(Outself_moms_combine_glmready,mom_nas_prop<0.5))
selfing_moms_nohigh_glm_null<-gam(cbind(count_self,count_out)~1,family = binomial(link = "logit"),data=filter(Outself_moms_combine_glmready,mom_nas_prop<0.5))
selfing_moms_nohigh_glm_HI<-gam(cbind(count_self,count_out)~nas_prop,family = binomial(link = "logit"),data=filter(Outself_moms_combine_glmready,mom_nas_prop<0.5))
selfing_moms_nohigh_glm_stream<-gam(cbind(count_self,count_out)~stream,family = binomial(link = "logit"),data=filter(Outself_moms_combine_glmready,mom_nas_prop<0.5))
selfing_moms_nohigh_glm_Year<-gam(cbind(count_self,count_out)~Year,family = binomial(link = "logit"),data=filter(Outself_moms_combine_glmready,mom_nas_prop<0.5))

lrtest(selfing_moms_nohigh_glm_null,selfing_moms_nohigh_glm_HI,selfing_moms_nohigh_glm_stream,selfing_moms_nohigh_glm_Year)
lrtest(selfing_moms_nohigh_glm_null,selfing_moms_nohigh_glm_stream,selfing_moms_nohigh_glm_Year)
lrtest(selfing_moms_nohigh_glm_null,selfing_moms_nohigh_glm_Year)
summary(selfing_moms_nohigh_glm_Year)

#Offspring HI (selfs only)
nrow(Outself_fruits_combine_lmready_selfs)
selfsonly_model_null<-lm(fruitmom_sway~1,data=Outself_fruits_combine_lmready_selfs)
selfsonly_model_HI<-lm(fruitmom_sway~mom_nas_prop,data=Outself_fruits_combine_lmready_selfs)
selfsonly_model_stream<-lm(fruitmom_sway~mom_nas_prop*stream,data=Outself_fruits_combine_lmready_selfs)
selfsonly_model_Year<-lm(fruitmom_sway~mom_nas_prop*stream*Year,data=Outself_fruits_combine_lmready_selfs)

lrtest(selfsonly_model_null,selfsonly_model_HI,selfsonly_model_stream,selfsonly_model_Year)



