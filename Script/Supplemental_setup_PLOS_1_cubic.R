################################################################################
################################################################################
#Supplemental analysis: setting up the dataset with cubic spline interpolation#
##Data set up: South Korea APC analysis#
#created: 5/18/2023#
#Updated: 10/25/2023 by YL#

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
  MESS,       # round prop & preserve sum to 100%
  data.table, #
  readxl,     # read excel file
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


#set up data

dt7913 <- read_excel("Data/SK_1979-2013.xlsx")
dt1421 <- read_excel("Data/SK_2014-2021.xlsx")
##setting up crime data (two periods)
table(dt7913$Classification4)
table(dt1421$Classification3)

#keep the following variables: 
#violence: assault, battery,murder, robbery, willful infliction of bodily injury
#property: frauds, intentional property damage, larceny, stolen goods
dt1=dt7913
dt1<-dt1%>%
  select(-Classification2)

dt1<-dt1%>%
  rename ("no"="No",                       "year"="Data Time",         "c1"="Classification",        
          "c2"="Classification3",          "c3"="Classification4",     "total"="Total",              "juv_tot"="Juvenile (subtotal)",
          "under14"="Under 14 years old",  "yr14"="14 years old",      "yr15"="15 years old",        "yr16"="16 years old",
          "yr17"="17 years old",           "yr18"="18 years old",      "adu_tot"="Adult (subtotal)", "yr19"="19 years old",      
          "yr20"="20 years old",           "yr21"="21 years old",      "yr22"="22 years old",        "yr23"="23 years old",
          "yr24"="24 years old",           "yr25"="25 years old",      "yr2630"="26-30 years old",   "yr3135"="31-35 years old",
          "yr3640"="36-40 years old",      "yr4150"="41-50 years old", "yr5160"="51-60 years old",   "yr6170"="61-70 years old",
          "yr71"="71 years old and over",  "unknown"="Unknown")
dt1$c3=tolower(dt1$c3)
dt1=subset(dt1, c3=="assault"  |c3=="battery"  |c3=="murder"  |c3== "robbery"|
             c3=="willful infliction of bodily injury"|
             c3=="frauds"   |c3=="intentional property damage"|
             c3== "larceny" |c3== "stolen goods"|
             c3== "punishment of violences, etc. act")
dt1=dt1%>%
  select(year, c3, yr14:yr18, yr19:yr6170)

dt1=dt1%>%
  mutate_at(vars(yr14:yr6170),readr::parse_number) #get rid of the "," in the numbers
#linear interpolation
#column 14: 20, age 15, 2630...6170
age2665=NULL
for (i in 1:350){
  a=as.numeric(dt1[i,14:20])
  a[2]=a[2]/5; a[3]=a[3]/5; a[4]=a[4]/5;a[5]=a[5]/10;a[6]=a[6]/10;a[7]=a[7]/10
  #b=approx(c(25,28,33,38,45,55,65),a,round(25:65))$y
  b=spline(c(25,28,33,38,45,55,65),a,length(25:65),method="natural")$y
  c=b[2:41]
  age2665=rbind(age2665,c)
}

new.var=paste0("yr",26:65)
colnames(age2665)=new.var
age2665=round(age2665,0)
dt1=cbind(dt1,age2665)
dt1=dt1%>%select(year, c3, yr15:yr25, yr26:yr64)


###############################################################################

#dt2 age categorization are different, so we have to do everything seperately
dt2=dt1421
dt2<-dt2%>%
  rename ("no"="No",                       "year"="Data Time",         "c1"="Classification1",        "c2"="Classification2", 
          "c3"="Classification3",          "total"="Total",            "juv_tot"="Juvenile (subtotal)",
          "under14"="Under 14 years old",  "yr14"="14 years old",      "yr15"="15 years old",        "yr16"="16 years old",
          "yr17"="17 years old",           "yr18"="18 years old",      "adu_tot"="Adult (subtotal)", "yr19"="19 years old",      
          "yr20"="20 years old",           "yr21"="21 years old",      "yr22"="22 years old",        "yr23"="23 years old",
          "yr24"="24 years old",           "yr25"="25 years old",      "yr2630"="26-30 years old",   "yr3135"="31-35 years old",
          "yr3640"="36~40 years old",      "yr4145"="41~45 years old", "yr4650"="46~50 years old",   "yr5155"="51~55 years old",       
          "yr5660"="56~60 years old",      "yr6170"="61~70 years old",
          "yr71"="71 years old and over",  "unknown"="Unknown")
dt2$c3=tolower(dt2$c3)
dt2=subset(dt2, c3=="assault"  |c3=="battery"  |c3=="murder"  |c3== "robbery"|
             c3=="willful infliction of bodily injury"|
             c3=="frauds"   |c3=="intentional property damage"|
             c3== "larceny" |c3== "stolen goods")

dt2=dt2%>%
  select(year, c3, yr14:yr18, yr19:yr6170)

dt2=dt2%>%
  mutate_at(vars(yr14:yr6170),readr::parse_number) #get rid of the "," in the numbers
#linear interpolation
#column 14: 22, age 15, 2630...6170
age2665=NULL
for (i in 1:72){
  a=as.numeric(dt2[i,14:22])
  a[2]=a[2]/5; a[3]=a[3]/5; a[4]=a[4]/5;a[5]=a[5]/5;a[6]=a[6]/5;a[7]=a[7]/5; a[8]=a[8]/5;a[9]=a[9]/10
  #b=approx(c(25,28,33,38,43,48,53,58,65),a,round(25:65))$y
  b=spline(c(25,28,33,38,43,48,53,58,65),a,length(25:65),method="natural")$y
  c=b[2:41]
  age2665=rbind(age2665,c)
}

new.var=paste0("yr",26:65)
colnames(age2665)=new.var
age2665=round(age2665,0)
dt2=cbind(dt2,age2665)
dt2=dt2%>%select(year, c3, yr15:yr25, yr26:yr64)

dt.full=rbind(dt1,dt2)
#after 2014 there is no punishment of violence#

###############################################
#combine indivdiual age to 5year age categories
dt.full= 
  dt.full %>% mutate(
    yr1519=rowSums(select(.,yr15:yr19)),
    yr2024=rowSums(select(.,yr20:yr24)),
    yr2529=rowSums(select(.,yr25:yr29)),
    yr3034=rowSums(select(.,yr30:yr34)),
    yr3539=rowSums(select(.,yr35:yr39)),
    yr4044=rowSums(select(.,yr40:yr44)),
    yr4549=rowSums(select(.,yr45:yr49)),
    yr5054=rowSums(select(.,yr50:yr54)),
    yr5559=rowSums(select(.,yr55:yr59)),
    yr6064=rowSums(select(.,yr60:yr64))
  ) 

dt.full<-dt.full%>%
  filter(year>1979&year<2020)
dt.full$year=as.numeric(dt.full$year)
dt.full<-dt.full%>%
  mutate(year5=cut(year, breaks=c(1980,1985,1990, 1995, 2000,2005,2010,2015,2020),right=F))
dt5=dt.full%>%
  select(year, c3, year5,yr1519, yr2024,yr2529, yr3034,yr3539, yr4044,yr4549,yr5054,yr5559,yr6064)

dt5a<-dt5%>%
  group_by(c3,year5)%>%
  summarise(across(c(yr1519, yr2024, yr2529, yr3034,yr3539,yr4044,yr4549,yr5054,yr5559,yr6064),mean))%>%
  as.data.frame

dt.l=dt5a%>% pivot_longer(cols=starts_with("yr"), names_to="age", values_to="count")  
dt.l<-dt.l%>%
  mutate(yr.n=recode(year5,
                     "[1980,1985)"=1982, "[1985,1990)"=1987, 
                     "[1990,1995)"=1992, "[1995,2000)"=1997,
                     "[2000,2005)"=2002, "[2005,2010)"=2007,
                     "[2010,2015)"=2012, "[2015,2020)"=2017),
         age.n=recode(age,
                      "yr1519"=17, "yr2024"=22, "yr2529"=27, "yr3034"=32, "yr3539"=37,
                      "yr4044"=42, "yr4549"=47, "yr5054"=52, "yr5559"=57, "yr6064"=62)
  )

#merging crime data with population data
pop <- read.csv("data/pop_data.csv")
dt.f=merge(dt.l,pop,by=c("yr.n","age.n"))
dt.f$count=round(dt.f$count,0)
dt.f$rate=dt.f$count/dt.f$tot_p*100000
dt.f$cohort=dt.f$yr.n-dt.f$age.n

write.csv(dt.f, file=file.path(dataDir,"data_spline.csv"), row.names=FALSE)

#end here: this is the full raw data set for SK age-specific arrests


