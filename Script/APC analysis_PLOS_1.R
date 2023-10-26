#APCI Analysis: South Korea APC analysis#
#created: 6/20/2023#
#Updated: 6/26/2023 by YL#
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


table1.list=NULL #table with arrest counts (demonstration)
table2a.list=NULL; table2b.list=NULL; table3a.list=NULL; table3b.list=NULL 
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
  
  #make a table with age and period arrest count summarize
  #Table 2b-Interaction terms
  #exporting cohort coefficients into a table#
  table1=data.frame(age.group=rep(c("15-19","20-24","25-29","30-34",
                                    "35-39","40-44","45-49","50-54"),each=2),
                    c1980=rep(NA,a.n*2),c1985=rep(NA,a.n*2),c1990=rep(NA,a.n*2),c1995=rep(NA,a.n*2),
                    c2000=rep(NA,a.n*2),c2005=rep(NA,a.n*2),c2010=rep(NA,a.n*2),c2015=rep(NA,a.n*2)
  )
  s=seq(1,p.n)
  for (i in s){
    odd=seq(1,15,2)
    even=seq(2,16,2)
    table1[odd,i+1]=dt$rate[(1+a.n*(i-1)):(a.n*i)]
    table1[even,i+1]=dt$cohort[(1+a.n*(i-1)):(a.n*i)]
  }
  
  table1.list[o]=list(table1)
  
  #sum to zero constraint
  options(contrasts=c("contr.sum","contr.sum"), 
          na.action = na.omit)
  #-------------------------------------------------------------------#
  
  #Step 2-------------------------------------------------------------
  #Step 2. Poison APC-I model, offset term to control for population: 
  
  api<-glm(count~age.f*period.f+offset(log(pop)), 
           family="poisson", data=dt)
  summary(api)
  model.list[o]=list(api)
  
  
  
  #cohort matrix
  cindex=array(rep(0,a.n*p.n),dim=c(a.n,p.n))
  for (j in 1:p.n){
    cindex[,j]=seq((a.n+j-1),j,-1)
  }
  
  
  #Step 3------------------------------------------------------
  #Step 3. GEt Age, Period main effect
  co2<-api$coefficients[1:c.n] #coefficients
  se2=summary(api)$coef[,2]#standard error
  a.co2=c(co2[2:a.n],-sum(co2[2:a.n])) #age coefficients, all ages sum to 0
  p.co2=c(co2[(a.n+1):c.n],-sum(co2[(a.n+1):c.n])) #period coefficients, all periods sum to 0
  
  #age effect SE#
  A=contr.sum(a.n,contrasts=TRUE)
  a.vcov=A%*%vcov(api)[2:a.n,2:a.n]%*%t(A)
  a.se=sqrt(diag(a.vcov))
  a.up=a.co2+2*a.se; a.lo=a.co2-2*a.se #95% cofidence interval
  
  #period effect SE#
  P=contr.sum(p.n,contrasts=TRUE)
  p.vcov=P%*%vcov(api)[(a.n+1):c.n,(a.n+1):c.n]%*%t(P)
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
  table2a=data.frame(effect=c(intercept,age.c4, period.c4), se=c(in.se,age.se4,period.se4), sig=c(in.sig,a.sig,p.sig))
  table2a$name=c("Intercept",names(table(dt$age)), names(table(dt$period.f)))
  table2a.list[o]=list(table2a)
  #-----------------------------------------------------------------------#
  
  #-----------------------------------------------------------------------#
  #Step 4. Cohort Interaction Terms#
  #computing full interaction index#
  T<-array(rep(0, a.n*p.n*(a.n-1)*(p.n-1)),dim=c(a.n*p.n, (a.n-1)*(p.n-1)))
  ind1=a.n*1:(p.n-1)
  ind2=(a.n*(p.n-1)+1):(a.n*p.n-1)
  ind3=a.n*p.n
  ind=c(ind1,ind2,ind3)       
  newind=1:(a.n*p.n)
  newind=newind[-ind]
  T[newind,]=diag((a.n-1)*(p.n-1))
  T[ind1,]=-diag(p.n-1)[,rep(1:(p.n-1),each=a.n-1)]
  T[ind2,]=-diag(a.n-1)[,rep(1:(a.n-1),p.n-1)]
  T[ind3,]=-rep(1,(a.n-1)*(p.n-1))
  
  #computing full interaction contrast#######
  api.c=api$coefficients
  iatemp=vcov(api)[(a.n+p.n):length(api.c),(a.n+p.n):length(api.c)]
  iavcov=T%*%iatemp%*%t(T)
  iaesti=as.vector(T%*%api.c[(a.n+p.n):length(api.c)]) #coefficient
  iase=sqrt(diag(iavcov)) 
  z=iaesti/iase
  pvalue=2*pnorm(-abs(z))
  
  cohortindex=as.vector(cindex)
  sig=rep(NA,a.n*p.n)
  for (i in 1:(a.n*p.n)){
    pp=pvalue[i]
    if (pp<0.05&pp>=0.01){sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){sig[i]="**"}else{
        if (pp<0.001){sig[i]="***"}else{
          sig[i]=""}
      }
    }
  }
  ia=data.frame(cohort=cohortindex, b=iaesti, se=iase,z, p=pvalue, sig=sig) 
  
  #Table 2b-Interaction terms
  #exporting cohort coefficients into a table#
  table2b=data.frame(age.group=c("15-19","20-24","25-29","30-34",
                                 "35-39","40-44","45-49","50-54"),
                     c1=rep(NA,a.n),p1=rep(NA,a.n),c2=rep(NA,a.n),p2=rep(NA,a.n),
                     c3=rep(NA,a.n),p3=rep(NA,a.n),c4=rep(NA,a.n),p4=rep(NA,a.n),
                     c5=rep(NA,a.n),p5=rep(NA,a.n),c6=rep(NA,a.n),p6=rep(NA,a.n),
                     c7=rep(NA,a.n),p7=rep(NA,a.n),c8=rep(NA,a.n),p8=rep(NA,a.n)
  )
  s=seq(2,p.n*2,by=2)
  for (i in s){
    table2b[,i]=iaesti[(1+(i/2-1)*a.n):(i/2*a.n)]
    table2b[,i+1]=ia$sig[(1+(i/2-1)*a.n):(i/2*a.n)]
  }
  
  table2b.list[o]=list(table2b)
  
  #-------------------------------------------------------------------------------#
  #4.1.inter-cohort deviations
  cint=rep(NA,c.n)
  cintse=rep(NA,c.n)
  cintz=rep(NA,c.n)
  cintp=rep(NA,c.n)
  
  for (k in 1:c.n){
    O=sum(cindex==k)
    k1=rep(1/O, O)
    k2=rep(0,a.n*p.n)
    k2[cindex==k]=k1
    
    contresti=k2%*%iaesti
    contrse=sqrt(t(k2)%*%iavcov%*%k2)
    z=contresti/contrse
    #p=2*pnorm(-abs(z))
    if(z>0){
      p=2*pnorm(z,lower.tail=F)
    }else{
      p=2*pnorm(z,lower.tail=T)
    }
    
    cint[k]=contresti
    cintse[k]=contrse
    cintz[k]=z
    cintp[k]=p
  }
  
  cgroup=seq(1,c.n)
  sig=rep(NA,c.n)
  for (i in 1:c.n){
    pp=cintp[i]
    if (pp<0.05&pp>=0.01){sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){sig[i]="**"}else{
        if (pp<0.001){sig[i]="***"}else{
          sig[i]=""}
      }
    }
  }
  #Table 2 a Intercohort deviation
  inter.cohort=data.frame(cohort=c(1:c.n), 
                          deviation=cint,s.e=cintse,
                          z=cintz,p=cintp,sig=sig)
  table3a.list[o]=list(inter.cohort)
  
  #---------------------------------------------------------------------#
  #4.2. intra-cohort: cohort slope variation
  cslope=rep(NA, c.n)
  cslopese=rep(NA, c.n)
  cslopez=rep(NA, c.n)
  cslopep=rep(NA,c.n)
  
  poly=1
  for(k in (poly+1):(c.n-poly)){
    O=sum(cindex==k)
    k1=contr.poly(O)
    k2=rep(0,a.n*p.n)
    k2[cindex==k]=k1[,poly]
    contresti=k2%*%iaesti
    contrse=sqrt(t(k2)%*%iavcov%*%k2)
    z=contresti/contrse
    
    if(z>0){
      p=2*pnorm(z,lower.tail=F)
      #pt(t,df,lower.tail=F)
    }else{
      p=2*pnorm(z,lower.tail=T)
      #pt(t,df,lower.tail=T)
    }
    
    cslope[k]=contresti
    cslopese[k]=contrse
    cslopez[k]=z
    cslopep[k]=p
  }
  
  sig=rep(NA,c.n)
  for (i in (poly+1):(c.n-poly)){
    pp=cslopep[i]
    if (pp<0.05&pp>=0.01){sig[i]="*"}else{
      if (pp<0.01&pp>=0.001){sig[i]="**"}else{
        if (pp<0.001){sig[i]="***"}else{
          sig[i]=""}
      }
    }
  }
  #table 3b. Intrac cohort slope
  cgroup=seq(1,c.n)
  cohortslope=data.frame(cohort=cgroup, slope=cslope, 
                         s.e=cslopese, z=cslopez, p=cslopep,sig.=sig)
  table3b.list[o]=list(cohortslope)
  
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
  
  #cohort effects
  cohortslope$slope=round(cohortslope$slope,digits=3)
  cohort.sig=NULL
  for (i in 1:c.n){
    cohort.sig[i]=paste(cohortslope$slope[i],cohortslope$sig.[i],sep="")
  }
  cohort.c4=as.numeric(cint);cohort.se4=as.numeric(cintse)
  #cohort=seq(1920,2000,by=5)
  c.df<-data.frame(effect=cohort.c4,cohort.n=cohort,
                   ci_up=cohort.c4+2*cohort.se4,
                   ci_lo=cohort.c4-2*cohort.se4,
                   intra=round(cohortslope$slope,digits=3),
                   intra.sig=cohort.sig)
  
  cohort.plot.l[[o]]=list(c.df)
  age.plot.l[[o]]=list(a.df)
  period.plot.l[[o]]=list(p.df)
  #----------------------------------------------------------------#
  
  #Figure 2--------------------------------------------------------
  # Figure 2a get data ready for figure 2a #heatmap
  f2=matrix(rep(NA,a.n*p.n),a.n, p.n)
  for (p in 1:p.n){
    f2[,p]=age.c4+period.c4[p]
  }
  #figure2a.list[[o]]=list(f2)
  
  # Figure 2b get data ready for figure 2b
  f2.d=as.matrix(subset(table2b,select=c(c1,c2,c3,c4,c5,c6,c7,c8)))
  #figure2b.list[[o]]=list(f2.d)
  #----------------------------------------------------------------------------#
  
  #Figure 3--------------------------------------------------------------------
  
  #Figure 3 cohort life course-------------------------------------------------#
  #assig period coefficients to each period
  #period1=seq(1978,2013,5)
  for (i in 1:p.n){
    j=period[i]
    dt$pco[dt$yr.n==j]=p.co2[i]
  } 
  #assign age coefficients to each age
  for (i in 1:a.n){
    j=age[i]
    dt$aco[dt$age.n==j]=a.co2[i]
  }
  #assign cohort interaction terms to each observation
  dt$cco=c(f2.d[,1],f2.d[,2],f2.d[,3],f2.d[,4],f2.d[,5],f2.d[,6],f2.d[,7],f2.d[,8])  
  #effects without interaction
  #dt$c1=exp(intercept+dt$aco+dt$pco)*100000
  #9/17: put effects together instead of predicted values
  #dt$c1=dt$aco+dt$pco
  #effects with interactionf
  #dt$c2=exp(intercept+dt$aco+dt$pco+dt$cco)*100000
  #9/17 replace predicted values with effects
  #dt$c2=dt$aco+dt$pco+dt$cco
  #dt$sig=as.character(ia$sig);dt$sig[dt$sig==""]="_"
  #set maximum values for y scale
  #ma=max(c(dt$c1,dt$c2))
  
  #Figure 4--------------------------------------------------------------------
  #Figure 4 cohort life course without period main effects---------------------------#
  #dt$c3=dt$aco+dt$cco
  
  #Figure 5--------------------------------------------------------------------
  #Figure 5 predicted arrest rates---------------------------#
  #age only
  dt$p1=exp(intercept+dt$aco)/(1+exp(intercept+dt$aco))*100000
  #with age+period
  dt$p2=exp(intercept+dt$aco+dt$pco)/(1+exp(intercept+dt$aco+dt$pco))*100000
  #with all effects prediction
  dt$p3=exp(intercept+dt$aco+dt$pco+dt$cco)/(1+exp(intercept+dt$aco+dt$pco+dt$cco))*100000
  
  #Figure 6--------------------------------------------------------------------
  #Figure 6 predicted arrest rates with ap only model---------------------------#
  
  #ap<-glm(count~age.f+period.f+offset(log(pop)), 
  #        family="poisson", data=dt)
  #summary(ap)
  #ap.in=coefficients(ap)[1]
  #co2<-ap$coefficients[1:c.n] #coefficients
  #se2=summary(api)$coef[,2]#standard error
  #a.co2=c(co2[2:a.n],-sum(co2[2:a.n])) #age coefficients, all ages sum to 0
  #p.co2=c(co2[(a.n+1):c.n],-sum(co2[(a.n+1):c.n])) #period coefficients, all periods sum to 0
  #period1=seq(1978,2013,5)
  #for (i in 1:p.n){
  #  j=period[i]
  #  dt$ap.pco[dt$yr.n==j]=p.co2[i]
  #} 
  #assign age coefficients to each age
  #for (i in 1:a.n){
  #  j=age[i]
  #  dt$ap.aco[dt$age.n==j]=a.co2[i]
  #}
  #dt$ap.p2=exp(ap.in+dt$ap.aco+dt$ap.pco)/(1+exp(ap.in+dt$ap.aco+dt$ap.pco))*100000
  ###############################################
  #age only model
  aonly<-glm(count~age.f+offset(log(pop)), 
             family="poisson", data=dt)
  summary(aonly)
  a.in=coefficients(aonly)[1]
  co2<-aonly$coefficients[1:a.n] #coefficients
  #se2=summary(api)$coef[,2]#standard error
  a.co2=c(co2[2:a.n],-sum(co2[2:a.n])) #age coefficients, all ages sum to 0
  #assign age coefficients to each age
  for (i in 1:a.n){
    j=age[i]
    dt$a.aco[dt$age==j]=a.co2[i]
  }
  dt$a.p2=exp(a.in+dt$a.aco)/(1+exp(a.in+dt$a.aco))*100000
  
  #dt$ac.p2=exp(intercept+dt$aco+dt$cco)/(1+exp(intercept+dt$aco+dt$cco))*100000
  
  ##############################################################################
  
  #get dataset ready
  figure3.list[[o]]=list(dt)
  
}

###############################################################################
#table 2 demonstrating the data structure 
tab1=data.frame(table1.list[3]) #homicide
write.csv(tab1,file=file.path(outDir,"table2_1005.csv"), row.names = FALSE)

################################################################################
#figure 2. estimated age effects, holding period and cohort constant

a1=data.frame(age.plot.l[[1]]) #Assault
a2=data.frame(age.plot.l[[2]]) #Homicide
a3=data.frame(age.plot.l[[3]]) #Ordinary property
a4=data.frame(age.plot.l[[4]]) #fraud


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
  #scale_y_continuous(limits = c(-1.3, 1))+
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

ggsave("Figure2_713.png", path = figDir, age.plot, width = 7.5, height = 7.5, dpi = 300, bg = 'white')
ggsave("Figure2_713.eps", path = figDir, age.plot, width = 7.5, height = 7.5, dpi = 300, bg = 'white')


###############################################################################
#now figure 3: predicted curve by age main effects & APC
plot3.list=NULL
off.name1=names(table(dt.5$c3));off.name1[1]="Assault"; off.name1[4]="Ordinary Theft"

for (i in 1:4){
  a=data.frame(figure3.list[[i]])
  a=subset(a,yr.n==1982|yr.n==2002|yr.n==2012)
  d1=data.frame(#rate=c(a$p1,a$rate),
    rate=c(a$p1,a$p3),
    age=rep(rep(age,3),2),
    yr=rep(rep(c("1980-84","2000-04","2010-14"),each=a.n),2),
    Model=rep(c("Age main effects","APC"),each=a.n*3))
  
  plot3=ggplot(d1, aes(x=age, y=rate, group=Model, shape=Model,colour=Model))+
    geom_line(aes(linetype=Model))+
    #geom_line()+
    geom_point()+
    facet_wrap(~yr, ncol=4)+
    #scale_y_continuous(limits = c(0, 12))+
    scale_x_continuous(limits = c(15, 55),breaks = seq(15,55, by=10))+
    xlab("Age") + ylab("Predicted Rate")+ 
    ggtitle(off.name1[i])+
    #theme_classic()+
    theme_bw()+
    scale_linetype_manual(values=c("dotted","solid"))+
    theme(panel.grid.minor = element_blank(),panel.grid.major= element_blank())+
    theme(legend.position="bottom")
  
  plot3.list[[i]]=plot3
}

plot3a=grid.arrange(plot3.list[[1]],plot3.list[[2]],plot3.list[[3]],plot3.list[[4]],
                    ncol=1,nrow=4)

ggsave("figure3_713.png",path = figDir, plot=plot3a,width=8.5,height=12, dpi = 300, bg = 'white')
ggsave("figure3_713.eps",path = figDir, plot=plot3a,width=8.5,height=12, dpi = 300, bg = 'white')

###############################################################################
#appendix. APC-I model results
#appendix 1. Main effects table
#Table 2a:Main effects
tb2.1=data.frame(table2a.list[1])
tb2.2=data.frame(table2a.list[2])
tb2.3=data.frame(table2a.list[3])
tb2.4=data.frame(table2a.list[4])

app.tb1=data.frame(var=tb2.1$name,ass.b=tb2.1$effect,ass.se=tb2.1$se, ass.sig=tb2.1$sig,
                   hom.b=tb2.2$effect,hom.se=tb2.2$se, hom.sig=tb2.2$sig,
                   pro.b=tb2.3$effect,pro.se=tb2.3$se, pro.sig=tb2.3$sig,
                   fra.b=tb2.4$effect,fra.se=tb2.4$se, fra.sig=tb2.4$sig
)
write.csv(app.tb1,file=file.path(outDir,"appendix_tb1_0713.csv"), row.names = FALSE)

#Table 2b:
tb2.1=data.frame(table2b.list[1])
tb2.2=data.frame(table2b.list[2])
tb2.3=data.frame(table2b.list[3])
tb2.4=data.frame(table2b.list[4])

write.csv(tb2.1,file=file.path(outDir,"appendix_tb2a.csv"), row.names = FALSE)
write.csv(tb2.2,file=file.path(outDir,"appendix_tb2b.csv"), row.names = FALSE)
write.csv(tb2.3,file=file.path(outDir,"appendix_tb2c.csv"), row.names = FALSE)
write.csv(tb2.4,file=file.path(outDir,"appendix_tb2d.csv"), row.names = FALSE)


#Table 3a: inter cohort deviaition
tb3.1=data.frame(table3a.list[1])
tb3.2=data.frame(table3a.list[2])
tb3.3=data.frame(table3a.list[3])
tb3.4=data.frame(table3a.list[4])

app.tb3=data.frame(cohort=cohort,
                   ass.b=tb3.1$deviation, ass.se=tb3.1$s.e,ass.sig=tb3.1$sig,
                   hom.b=tb3.2$deviation, hom.se=tb3.2$s.e,hom.sig=tb3.2$sig,
                   pro.b=tb3.3$deviation, pro.se=tb3.3$s.e,pro.sig=tb3.3$sig,
                   fra.b=tb3.4$deviation, fra.se=tb3.4$s.e,fra.sig=tb3.4$sig
)
write.csv(app.tb3,file=file.path(outDir,"appendix_tb3.csv"), row.names = FALSE)

###############################################################################
#Appendix Period plot
#period plot 1:
p1=data.frame(period.plot.l[[1]]) #assault
p2=data.frame(period.plot.l[[2]]) #homicide
p3=data.frame(period.plot.l[[3]]) #property
p4=data.frame(period.plot.l[[4]]) #fraud


p.df1=data.frame(effect=c(p1$effect,p2$effect,p3$effect,p4$effect),
                 year=rep(period,4),
                 ci_up=c(p1$ci_up,p2$ci_up,p3$ci_up,p4$ci_up),
                 ci_lo=c(p1$ci_lo,p2$ci_lo,p3$ci_lo,p4$ci_lo),
                 Offense=c(rep("1.Assault",8),
                           rep("2.Homicide",8),
                           rep("3.Ordinary Theft",8),
                           rep("4.Fraud",8)))

period.plot1=
  ggplot(p.df1, aes(x=year, y=effect))+
  facet_wrap(~Offense, ncol=4,nrow=1)+
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_up),width=1)+
  scale_x_continuous(breaks = seq(1982,2017,5),
                     labels=c("1980-84",
                              "1985-89","1990-94","1995-99","2000-04","2005-09",
                              "2010-15", "2015-19"))+
  #scale_y_continuous(limits = c(-1, 1))+
  geom_point() +
  geom_line( )+
  geom_hline(yintercept = 0)+
  xlab("Period") + 
  ylab("Period Effect")+  #ggtitle("1b")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.position="", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#cohort plot 1:
c1=data.frame(cohort.plot.l[[1]])
c2=data.frame(cohort.plot.l[[2]])
c3=data.frame(cohort.plot.l[[3]])
c4=data.frame(cohort.plot.l[[4]])


c.df1=data.frame(effect=c(c1$effect,c2$effect,c3$effect,c4$effect),
                 cohort=rep(cohort,4),
                 intra=c(c1$intra,c2$intra,c3$intra,c4$intra),
                 ci_up=c(c1$ci_up,c2$ci_up,c3$ci_up,c4$ci_up),
                 ci_lo=c(c1$ci_lo,c2$ci_lo,c3$ci_lo,c4$ci_lo),
                 Offense=c(rep("1.Assault",c.n),
                           rep("2.Homicide",c.n),
                           rep("3.Ordinary Theft",c.n),
                           rep("4.Fraud",c.n)))

#c.df1=subset(c.df1,cohort>=1925&cohort<=1990)


cohort.plot1=
  ggplot(c.df1, aes(x=cohort, y=effect))+
  facet_wrap(~Offense, ncol=4,nrow=1)+
  #geom_errorbar(aes(ymin=ci_lo, ymax=ci_up), width=1)+
  scale_x_continuous(breaks = seq(1930,2000,5))+
  #scale_y_continuous(limits = c(-1, 1))+
  geom_point() +
  geom_line()+
  geom_hline(yintercept = 0)+
  xlab("Cohort") + 
  ylab("Inter-cohort Deviation")+  #ggtitle("1c")+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(legend.position="", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file=file.path(figDir,"app_period_plot_0713.png"),period.plot1,width=11,height=3.5)
ggsave(file=file.path(figDir,"app_cohort_plot_0713.png"),cohort.plot1,width=11,height=3.5)
