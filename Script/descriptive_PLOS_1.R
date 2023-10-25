#Descriptive analysis: South Korea APC analysis#
#created: 6/20/2023#
#Updated: 6/26/2023 by YL#
#updated: 7/13/2023 Dropped the last age group--stopped at age 54#


library(extrafont)
font_import()
loadfonts(device = "win")


# Installs and loads packages automatically
library("pacman")                  # Load pacman package

# Install packages not yet installed & load them
pacman::p_load(
  here,       # relative file paths
  readr,      # for parsing out numbers from a variable
  foreign,    # read data
  plyr,       #
  dplyr,      # variable processing
  tidyr,      # reshaping data
  purrr,      # to use modify_at
  moments,     # calculate descriptive statistics such as skewness
  MESS,       # round prop & preserve sum to 100%
  data.table, #
  gtsummary,  # pretty weighted tables
  ggplot2,    # graphing  
  conflicted) # choose default packages


# Address any conflicts in the packages
conflict_scout() # identify the conflicts
conflict_prefer("here", "here")
conflict_prefer("mutate", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("select","dplyr")

# Set-up the Directories -------------------------------------------------------

## Set the project directory to the current working directory.
projDir <- here::here()  # File path to this project's directory
dataDir<-file.path(projDir,"Data")    # File path to where data will be downloaded
outDir  <- "output"                   # Name of the sub-folder where we will save results
figDir  <- file.path(outDir, "figs")  # Name of the sub-folder where we will save generated figures


## This will create sub-directory folders in the master project directory if doesn't exist
if (!dir.exists(here::here(outDir))){
  dir.create(outDir)
} else {
  print("Output directory already exists!")
}

if (!dir.exists(here::here(figDir))){
  dir.create(figDir)
} else {
  print("Figure directory already exists!")
}

###############################################################################
################################################################################
dt.f <- read.csv("Data/data.csv")
#6/23/2023: combine some of the offense categories
# assault=assault+battery+injury 
# fraud = fraud
# ordinary theft=larceny+robbery
# murder = murder

#combine assault
d1=subset(dt.f, c3=="assault")
d2=subset(dt.f, c3=="battery")
d3=subset(dt.f, c3=="willful infliction of bodily injury")

ass1=d1$count+d2$count+d3$count

vio1=d1; vio1$count=ass1;vio1$c3="Assault, Battery and Injury"
vio1$rate=vio1$count/vio1$tot_p*100000

#combine ordinary theft
d1=subset(dt.f, c3=="larceny")
d2=subset(dt.f, c3=="robbery")

prop.count=d1$count+d2$count
prop=d1; prop$count=prop.count; prop$c3="Ordinary Property"
prop$rate=prop$count/prop$tot_p*100000

fraud=subset(dt.f,c3=="frauds")
fraud$c3="Fraud"
murder=subset(dt.f,c3=="murder")
murder$c3="Homicide"

dt.5=rbind(vio1,prop,fraud, murder)

#now the dataframe is set up for analysis#
#--------------------------------------------------------------------------------
  #for 4/5 offenses only"dt.5", graphing for 1980, 2000, 2020
  #calculate PAI
  #sum rates by offense types and by years
  #calculate PAI
dt.5<-subset(dt.5, age.n<57)
dt.5<-dt.5%>% 
  group_by(c3, yr.n)%>%
  mutate(t.rate=sum(rate))
dt.5$pai<-dt.5$rate/dt.5$t.rate*100

#reorder offenses
dt.5$c3=factor(dt.5$c3, levels=c("Assault, Battery and Injury", 
                                 "Homicide","Ordinary Property", "Fraud"))

#select three periods only for presentation
a= subset(dt.5,yr.n==1982|yr.n==2002|yr.n==2012)
a$Period=NA
a$Period[a$year5=="[1980,1985)"]="1980 [1980-84]"
a$Period[a$year5=="[2000,2005)"]="2000 [2000-04]"
a$Period[a$year5=="[2010,2015)"]="2010 [2010-14]"


##Figure 1. Histogram##########################################################
off.name=names(table(dt.5$c3))
off.name1=names(table(dt.5$c3));off.name1[1]="Assault"; off.name1[4]="Ordinary Theft"
figure1a.list=NULL
for (i in 1:4){
  b=off.name[i]
  dt=subset(a, c3==b)
  c=subset(dt,yr.n==1982|yr.n==2002|yr.n==2012,select=c(yr.n,age.n,pai))
  c$Period=rep(c("1980-84","2000-04","2010-14"),each=8)
  
  plot=ggplot(c, aes(x=age.n, y=pai))+
    geom_bar(stat="identity") +
    facet_wrap(~Period, ncol=4)+
    scale_x_continuous(#limits=c(15,54),
      breaks = seq(17,55, by=5),
      labels=c("15-19","20-24","25-29","30-34",
               "35-39","40-44","45-49","50-54"))+
    #scale_y_continuous(limits = c(0, 35))+
    
    #geom_line(aes(linetype=Period))+
    xlab("Age") + ylab("PAI")+  
    ggtitle(off.name1[i])+ 
    theme_bw()+
    theme(panel.grid.minor = element_blank())+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(legend.position="")
  figure1a.list[i]=list(plot)
  
  
}

library(gridExtra)
figure1=grid.arrange(figure1a.list[[1]],figure1a.list[[2]],
                     figure1a.list[[3]], figure1a.list[[4]],
                     ncol=1,nrow=4)

#ggsave("Figure1_1024.png", path = figDir, figure1, width = 7.5, height = 10, dpi = 300, bg = 'white')
ggsave("Figure1_1024.eps", path = figDir, figure1, width = 7.5, height = 10, dpi = 300, bg = 'white')


################################################################################
#identify peak age, 1/2 peak, and skewness for the three periods
a= subset(dt.5,yr.n==1982|yr.n==2002|yr.n==2012)
#peak
a=a%>%
  group_by(c3,yr.n)%>%
  mutate(peak=max(pai))

peak.age<-a%>%
  group_by(c3,yr.n)%>%
  filter(pai==peak)

peak.age$peak.a=peak.age$age.n

#1/2 peak
a$half.peak=a$peak/2
peak.age<-peak.age%>%
  select(yr.n,c3,peak.a)
a=merge(a,peak.age,by=c("yr.n","c3"))
half.p<-a%>%
  filter(pai<=half.peak&age.n>=peak.a)
half.p<-half.p%>%
  group_by(c3,yr.n)%>%
  summarise(half.peak=min(age.n))
#mutate(half.peak=min(age.n))

#skewness
skewness=a%>%
  group_by(c3,yr.n)%>%
  summarise(skew=skewness(pai))

#merge the three
table1=merge(peak.age,  half.p, all=TRUE, by=c("c3","yr.n"))
table1=merge(table1, skewness, all=TRUE, by=c("c3","yr.n"))

dataDir<-file.path(projDir,"Data")    # File path to where data will be downloaded
outDir  <- "output" 
write.csv(table1, file=file.path(outDir,"table_1024.csv"), row.names=FALSE)



