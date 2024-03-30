
#Supplemental analysis: comparing AP Poisson vs AP Negative Binomial Models
#APCI Analysis: South Korea APC analysis#
#created: 6/20/2023#
#Updated: 1/20/2024 by YL#
################################################################################
#modeling
pacman::p_load(
  ggplot2,
  gnm,
  lme4,
  MASS, 
  Matrix,
  gridExtra
) 
#library(ggplot2);library(gnm);library(lme4)
#library(MASS);library(Matrix);library(gridExtra)
#library(gtable);library(grid);library(gridExtra);library(xlsx)


#table1.list=NULL #table with arrest counts (demonstration)
table2a.list=NULL
table2b.list=NULL
age.plot.l1=NULL; period.plot.l1=NULL
age.plot.l2=NULL; period.plot.l2=NULL


; table3a.list=NULL; table3b.list=NULL 
cohort.plot.l=NULL; age.plot.l=NULL; period.plot.l=NULL
figure2a.list=NULL; figure2b.list=NULL
figure3.list=NULL; model.list=NULL

#reorder offenses
dt.5$c3=factor(dt.5$c3, levels=c("Assault, Battery and Injury", 
                                 "Homicide","Ordinary Property", "Fraud"))
off.name=names(table(dt.5$c3))
#level: assault, battery & injury, battery and injury;
period=seq(1982,2017,by=5)
age=seq(17,52,by=5)
cohort=seq(1930,2000,5)
for (o in 1:4){
  a=off.name[o]
  dt=subset(dt.5, c3==a)
  dt$age.f=factor(dt$age.n)
  dt$period.f=factor(dt$yr.n)
  dt$pop=dt$tot_p
  
  a.n<-length(levels(dt$age.f)) #number of age groups
  p.n<-length(levels(dt$period.f))  #number of periods
  c.n<-a.n+p.n-1 #number of cohorts
  
  
  
  #sum to zero constraint
  options(contrasts=c("contr.sum","contr.sum"), 
          na.action = na.omit)
  #-------------------------------------------------------------------#
  
  #Step 2-------------------------------------------------------------
  #Step 2. Comparing poisson vs. negative binomial model
  
  ap1<-glm(count~age.f+period.f+offset(log(pop)), 
             family="poisson", data=dt)
  summary(ap1)
  
  ap2<-glm.nb(count~age.f+period.f+offset(log(pop)), 
                data=dt)
  summary(ap2)
  
  
  #Step 3------------------------------------------------------
  #Step 3b. GEt Age, Period main effect for ap2
  co2<-ap1$coefficients[1:c.n] #coefficients
  se2=summary(ap1)$coef[,2]#standard error
  a.co2=c(co2[2:a.n],-sum(co2[2:a.n])) #age coefficients, all ages sum to 0
  p.co2=c(co2[(a.n+1):c.n],-sum(co2[(a.n+1):c.n])) #period coefficients, all periods sum to 0
  
  #age effect SE#
  A=contr.sum(a.n,contrasts=TRUE)
  a.vcov=A%*%vcov(ap1)[2:a.n,2:a.n]%*%t(A)
  a.se=sqrt(diag(a.vcov))
  a.up=a.co2+2*a.se; a.lo=a.co2-2*a.se #95% cofidence interval
  
  #period effect SE#
  P=contr.sum(p.n,contrasts=TRUE)
  p.vcov=P%*%vcov(ap1)[(a.n+1):c.n,(a.n+1):c.n]%*%t(P)
  p.se=sqrt(diag(p.vcov))
  p.up=p.co2+2*p.se; p.lo=p.co2-2*p.se #95% cofidence interval
  
  #combine APC mixed model information
  api2.coef=c(co2[1],a.co2,p.co2) #age, period, main effects
  api2.se=c(se2[1],a.se,p.se) #age, period, main effect se
  cbind(api2.coef, api2.se)
  api2.aci=cbind(a.up, a.lo) #confidence interval of age
  api2.pci=cbind(p.up, p.lo) #confidence interval of period
  
  #for Table 2a. age period main effects
  #get age, period, cohort effects ready#
  age.c4=as.numeric(a.co2);period.c4=as.numeric(p.co2)
  age.se4=a.se; period.se4=p.se
  
  #Intercept#
  intercept=as.numeric(ap1$coefficients[1])
  in.se=summary(ap1)$coefficients[1,2]
  in.p=2*pnorm(-abs(intercept/in.se))#pvalues for constant
  in.sig=NA
  if (in.p<0.05&in.p>=0.01){in.sig="*"}else{
    if (in.p<0.01&in.p>=0.001){in.sig="**"}else{
      if (in.p<0.001){in.sig="***"}else{
        in.sig=""}}}
  #age effect#
  a.z=age.c4/age.se4;a.p=2*pnorm(-abs(as.numeric(a.z))) #pvalues for age effects
  a.sig=rep(NA,a.n)
  for (i in 1:a.n){
    pp=a.p[i]
    if (pp<0.05&pp>=0.01){a.sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){a.sig[i]="**"}else{
        if (pp<0.001){a.sig[i]="***"}else{
          a.sig[i]=""}}}}
  #period effect#
  p.z=period.c4/period.se4;p.p=2*pnorm(-abs(as.numeric(p.z))) #pvalues for period effects
  p.sig=rep(NA,p.n)
  for (i in 1:p.n){
    pp=p.p[i]
    if (pp<0.05&pp>=0.01){p.sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){p.sig[i]="**"}else{
        if (pp<0.001){p.sig[i]="***"}else{
          p.sig[i]=""}}}}
  #combine table: 
  table2a=data.frame(effect=c(intercept,age.c4, period.c4), se=c(in.se,age.se4,period.se4), sig=c(in.sig,a.sig,p.sig))
  table2a$name=c("Intercept",names(table(dt$age)), names(table(dt$period.f)))
  table2a.list[o]=list(table2a)
  
  
  #Figure 1-------------------------------------------------#
  #Figure 1 Age and Period main effects#GEt data ready
  
  #may need to add code to plot cohort plots
  #making age, period, cohort plots
  #period effects
  
  p.df<-data.frame(effect=period.c4,year=period,
                   ci_up=period.c4+2*period.se4,
                   ci_lo=period.c4-2*period.se4)
  #age effects
  a.df<-data.frame(effect=age.c4,age.n=age,
                   ci_up=age.c4+2*age.se4,
                   ci_lo=age.c4-2*age.se4)
  
  age.plot.l1[[o]]=list(a.df)
  period.plot.l1[[o]]=list(p.df)
  #----------------------------------------------------------------#
  #-----------------------------------------------------------------------#
  
  #Step 3a. GEt Age, Period main effect for ap2
  co2<-ap2$coefficients[1:c.n] #coefficients
  se2=summary(ap2)$coef[,2]#standard error
  a.co2=c(co2[2:a.n],-sum(co2[2:a.n])) #age coefficients, all ages sum to 0
  p.co2=c(co2[(a.n+1):c.n],-sum(co2[(a.n+1):c.n])) #period coefficients, all periods sum to 0
  
  #age effect SE#
  A=contr.sum(a.n,contrasts=TRUE)
  a.vcov=A%*%vcov(ap2)[2:a.n,2:a.n]%*%t(A)
  a.se=sqrt(diag(a.vcov))
  a.up=a.co2+2*a.se; a.lo=a.co2-2*a.se #95% cofidence interval
  
  #period effect SE#
  P=contr.sum(p.n,contrasts=TRUE)
  p.vcov=P%*%vcov(ap2)[(a.n+1):c.n,(a.n+1):c.n]%*%t(P)
  p.se=sqrt(diag(p.vcov))
  p.up=p.co2+2*p.se; p.lo=p.co2-2*p.se #95% cofidence interval
  
  #combine APC mixed model information
  api2.coef=c(co2[1],a.co2,p.co2) #age, period, main effects
  api2.se=c(se2[1],a.se,p.se) #age, period, main effect se
  cbind(api2.coef, api2.se)
  api2.aci=cbind(a.up, a.lo) #confidence interval of age
  api2.pci=cbind(p.up, p.lo) #confidence interval of period
  
  #for Table 2a. age period main effects
  #get age, period, cohort effects ready#
  age.c4=as.numeric(a.co2);period.c4=as.numeric(p.co2)
  age.se4=a.se; period.se4=p.se
  
  #Intercept#
  intercept=as.numeric(api$coefficients[1])
  in.se=summary(api)$coefficients[1,2]
  in.p=2*pnorm(-abs(intercept/in.se))#pvalues for constant
  in.sig=NA
  if (in.p<0.05&in.p>=0.01){in.sig="*"}else{
    if (in.p<0.01&in.p>=0.001){in.sig="**"}else{
      if (in.p<0.001){in.sig="***"}else{
        in.sig=""}}}
  #age effect#
  a.z=age.c4/age.se4;a.p=2*pnorm(-abs(as.numeric(a.z))) #pvalues for age effects
  a.sig=rep(NA,a.n)
  for (i in 1:a.n){
    pp=a.p[i]
    if (pp<0.05&pp>=0.01){a.sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){a.sig[i]="**"}else{
        if (pp<0.001){a.sig[i]="***"}else{
          a.sig[i]=""}}}}
  #period effect#
  p.z=period.c4/period.se4;p.p=2*pnorm(-abs(as.numeric(p.z))) #pvalues for period effects
  p.sig=rep(NA,p.n)
  for (i in 1:p.n){
    pp=p.p[i]
    if (pp<0.05&pp>=0.01){p.sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){p.sig[i]="**"}else{
        if (pp<0.001){p.sig[i]="***"}else{
          p.sig[i]=""}}}}
  #combine table: 
  table2b=data.frame(effect=c(intercept,age.c4, period.c4), se=c(in.se,age.se4,period.se4), sig=c(in.sig,a.sig,p.sig))
  table2b$name=c("Intercept",names(table(dt$age)), names(table(dt$period.f)))
  table2b.list[o]=list(table2b)
  #-----------------------------------------------------------------------#
  
  
  #Figure 1-------------------------------------------------#
  #Figure 1 Age and Period main effects#GEt data ready
  
  #may need to add code to plot cohort plots
  #making age, period, cohort plots
  #period effects
  
  p.df<-data.frame(effect=period.c4,year=period,
                   ci_up=period.c4+2*period.se4,
                   ci_lo=period.c4-2*period.se4)
  #age effects
  a.df<-data.frame(effect=age.c4,age.n=age,
                   ci_up=age.c4+2*age.se4,
                   ci_lo=age.c4-2*age.se4)
  
  age.plot.l2[[o]]=list(a.df)
  period.plot.l2[[o]]=list(p.df)
  #----------------------------------------------------------------#
  

  
}

################################################################################
#figure 2. estimated age effects, holding period and cohort constant

a1=data.frame(age.plot.l1[[1]]) #Assault
a2=data.frame(age.plot.l1[[2]]) #Homicide
a3=data.frame(age.plot.l1[[3]]) #Ordinary property
a4=data.frame(age.plot.l1[[4]]) #fraud


a.df1=data.frame(effect=c(a1$effect,a2$effect,
                          a3$effect,a4$effect),
                 age=rep(age,4),
                 ci_up=c(a1$ci_up,a2$ci_up,
                         a3$ci_up,a4$ci_up),
                 ci_lo=c(a1$ci_lo,a2$ci_lo,
                         a3$ci_lo,a4$ci_lo),
                 Offense=c(rep("Assault",a.n),
                           rep("Homicide",a.n),
                           rep("Ordinary theft",a.n),
                           rep("Fraud",a.n)))
a.df1$Offense=factor(a.df1$Offense,levels=c("Assault","Homicide","Ordinary theft","Fraud"))

age.plot=
  ggplot(a.df1, aes(x=age, y=effect))+
  facet_wrap(~Offense, ncol=2,nrow=2)+
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_up),width=1)+
  scale_x_continuous(breaks = seq(17,62, by=5),
                     labels=c("15-19","20-24","25-29","30-34",
                              "35-39","40-44","45-49","50-54",
                              "55-59","60-64"))+
  scale_y_continuous(limits = c(-1.5, 1.5))+
  geom_point() +
  geom_line()+
  geom_hline(yintercept = 0)+
  xlab("Age") + ylab("Age Effect")+ 
  #ggtitle("a.")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.position="bottom", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position=c(0.8, 0.8))

ggsave("Figure2_ap_poisson.png", path = figDir, age.plot, width = 7.5, height = 7.5, dpi = 300, bg = 'white')
ggsave("Figure2_ap_poisson.eps", path = figDir, age.plot, width = 7.5, height = 7.5, dpi = 300, bg = 'white')


################################################################################
#figure 2. estimated age effects, holding period and cohort constant

a1=data.frame(age.plot.l2[[1]]) #Assault
a2=data.frame(age.plot.l2[[2]]) #Homicide
a3=data.frame(age.plot.l2[[3]]) #Ordinary property
a4=data.frame(age.plot.l2[[4]]) #fraud


a.df1=data.frame(effect=c(a1$effect,a2$effect,
                          a3$effect,a4$effect),
                 age=rep(age,4),
                 ci_up=c(a1$ci_up,a2$ci_up,
                         a3$ci_up,a4$ci_up),
                 ci_lo=c(a1$ci_lo,a2$ci_lo,
                         a3$ci_lo,a4$ci_lo),
                 Offense=c(rep("Assault",a.n),
                           rep("Homicide",a.n),
                           rep("Ordinary theft",a.n),
                           rep("Fraud",a.n)))
a.df1$Offense=factor(a.df1$Offense,levels=c("Assault","Homicide","Ordinary theft","Fraud"))

age.plot=
  ggplot(a.df1, aes(x=age, y=effect))+
  facet_wrap(~Offense, ncol=2,nrow=2)+
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_up),width=1)+
  scale_x_continuous(breaks = seq(17,62, by=5),
                     labels=c("15-19","20-24","25-29","30-34",
                              "35-39","40-44","45-49","50-54",
                              "55-59","60-64"))+
  scale_y_continuous(limits = c(-1.5, 1.7))+
  geom_point() +
  geom_line()+
  geom_hline(yintercept = 0)+
  xlab("Age") + ylab("Age Effect")+ 
  #ggtitle("a.")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.position="bottom", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position=c(0.8, 0.8))

ggsave("Figure2_ap_nb.png", path = figDir, age.plot, width = 7.5, height = 7.5, dpi = 300, bg = 'white')
ggsave("Figure2_ap_nb.eps", path = figDir, age.plot, width = 7.5, height = 7.5, dpi = 300, bg = 'white')

