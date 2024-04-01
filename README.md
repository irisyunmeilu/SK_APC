# SK_APC
South Korea age-specific arrest data from 1980-2019 to examine the age, period, and cohort effects on crime. 
This folder includes scripts and data for the publication: "Examining the stability and change in age-crime relation in South Korea, 1980–2019: An age-period-cohort analysis" at PLOS One. 
The original dataset: "SK_1979-2013.xlsx" & "SK_2014-2021.xlsx" are downloaded from Analytical Statistics on Crime via the Crime and Criminal  Justice Statistics data portal: https://www.kicj.re.kr/crimestats/portal/stat/easyStatEngPage.do;jsessionid=Sx2vu5MPelzzlzr8YQBdejVQYipYpgZM6Dklh3A96mSbaWCwJ6qN8YlU21h4NS7Q.ccjs_web_servlet_engine1
"pop_data.csv" is the age-period-sepcific population data in South Korea
"data.csv" is the data set up for APC-I mdoel (see the script "setup_PLOS_1.R" for more details)
"data_spline.csv" & "data_postad.csv" are data set up for supplemental analysis. "data_spline.csv" is the data set up with cubic spline interpolation; and "data_postad.csv" is the data set up for supplemental analysis using post-adolescence age-crime data only. 
Script: 
- "setup_PLOS_1.R" sets up the data;
- "descriptive_PLOS_1.R" conducts the descriptive analysis Figure 1 & Table 1 in the paper;
- "APC analysis_PLOS_1.R" includes code for the APC analysis and Figure 2 & 3 in the paper (and other information in Appendix)
- the other files in the script folder are scripts for the supplemental analysis. 
