##########################################################
#
#  EXAMPLE usage of LWP-Trends functions
#
##########################################################

#This file provides three example applciations of the LWP-Trends code
#The three applications are for:
#   1. Monthly observations (a) Seasonal and (b) non-seasonal
#   2. Monthly observations with gaps
#   3. Annual Observations
#IN addition we include:
#   4. Trend aggregation examples


# Date: August 2024
# By Ton Snelder and Caroline Fraser
##########################################################

rm(list = ls())                                               # clear the memory
##########################################################
# 0.  Preparation
##########################################################
# Read in required packages
##########################################################
require(tidyverse)
#require(dplyr)
#require(ggpubr)   #this is just to plot ggplot objects in one picture, not actually required by LWPtrends

##########################################################
# Read in functions and sample data
##########################################################

#dir.code<-"C:\\LWP\\Harris Consulting\\Common Data and Code - Documents\\WQAnalysis\\R_TrendAnalysisFunctions\\"#UPDATE TO LOCAL DIRECTORY
#setwd(dir.code)

source("./LWPTrends_v2502/LWPTrends_v2502.r")

junk<-load("LWPTrends_ExampleData_v2502.rdata")


#########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         EXAMPLE 1: NonSeasonal Monthly Observations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#########################################################
#  Preliminary Setup
#########################################################
#Format myDate
WQData_Ex1$myDate <- as.Date(as.character(WQData_Ex1$sdate),"%Y-%m-%d")

#Add on time increment and extra date information - NOTE this dataset is based on Water Years
WQData_Ex1<-GetMoreDateInfo(WQData_Ex1,firstMonth = 7)

#Process censored values
WQData_Ex1 <- RemoveAlphaDetect(WQData_Ex1,ColToUse="Value")

summary(WQData_Ex1) # note there are no censored values (Censored == TRUE)

#########################################################
# Inspect Data
#########################################################

library(ggpubr)
#Make time series plot
Out1<-InspectTrendData(WQData_Ex1, EndYear = 2017, FlowtoUse="finalQ",
                        ReturnALLincr=TRUE,do.plot = TRUE,mymain="Example 1")
#Look at the numeric outputs:
Out1[[2]]
#look at the plots:
graphics.off();x11(width=30,height=25);ggarrange(Out1[[3]][[1]],ggarrange(Out1[[3]][[2]],Out1[[3]][[3]],Out1[[3]][[4]],nrow=1,align="h"),nrow=2)

WQData_Ex1<-Out1[[1]]
#If you didn't want the plots, can just get the data with teh selected time increment
#WQData_Ex1<-InspectTrendData(WQData_Ex1, EndYear = 2017)

##########################################################
# Seasonality Test on Raw Data
##########################################################
#Detemrine the appropriate season for the data (adds on column "Season") and conducts seasonality test
Seas1<-GetSeason(WQData_Ex1,ValuesToUse = "RawValue",printKW = TRUE,mymain="Example 1",do.plot = TRUE)

#Look at teh numeric outputs
Seas1[[2]] #The raw data is NOT seasonal
#look at the plot outputs
x11();Seas1[[3]]
WQData_Ex1<-Seas1[[1]]
##########################################################
# Perform Trend Tests on Raw Data
##########################################################
Trend1<-NonSeasonalTrendAnalysis(WQData_Ex1,mymain="Ex 1 Raw Trend",Year="CustomYear",do.plot=T)
Trend1[[1]]
x11();Trend1[[2]]



#########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         EXAMPLE 2: Seasonal Monthly Observations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Example dataset has already had censoring information added

WQData_Ex2<-GetMoreDateInfo(WQData_Ex2,firstMonth = 7)

#########################################################
# Inspect Data
#########################################################
#Make time series plot

#Make time series plot
Out2<-InspectTrendData(WQData_Ex2, EndYear = 2022, Year = "CustomYear",
                        ReturnALLincr=TRUE,do.plot = TRUE,mymain="Example 2")
lib#Look at the numeric outputs about data availablity for different time increments (note, using default of 90% of time increments)
Out2[[2]]
#look at the plots:
x11(width=30,height=25);ggarrange(Out2[[3]][[1]],ggarrange(Out2[[3]][[2]],Out2[[3]][[3]],Out2[[3]][[4]],nrow=1,align="h"),nrow=2)

WQData_Ex2<-Out2[[1]]
#If you didn't want the plots, can just get the data with the selected time increment
#WQData_Ex1a<-InspectTrendData(WQData_Ex1a, EndYear = 2017)

##########################################################
# get Season and perform Seasonality Test on Raw Data
##########################################################
#Detemrine the appropriate season for the data (adds on column "Season")
# and conduct seasonality test
Seas2<-GetSeason(WQData_Ex2,printKW = TRUE,mymain="Example 2: Raw Data",do.plot=TRUE)
Seas2[[2]];x11();Seas2[[3]]
#The raw data is strongly seasonal

#Get the data with the Season column added
WQData_Ex2<-Seas2[[1]]

##########################################################
# Perform Trend Tests on Raw Data
##########################################################
Trend_ex2<-SeasonalTrendAnalysis(WQData_Ex2,mymain="Ex 2 Raw Trend",do.plot=T)
Trend_ex2[[1]];x11();Trend_ex2[[2]]


#########################################################



#########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         EXAMPLE 2: Limited Monthly Observations (Quarterly analysis)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#########################################################
# Preliminary Setup
#########################################################
#Format myDate
WQData_Ex3$myDate <- as.Date(as.character(WQData_Ex3$sdate),"%Y-%m-%d")

#Add on seasons and extra date informaiton - NOTE, using a custom year from Jul-Jun
WQData_Ex3<-GetMoreDateInfo(WQData_Ex3,firstMonth = 7)

#Process censored values
WQData_Ex3 <- RemoveAlphaDetect(WQData_Ex3,ColToUse="Value")

#########################################################
# Inspect Data
#########################################################
Out3<-InspectTrendData(WQData_Ex3,  EndYear = 2017, FlowtoUse="finalQ",
                  ReturnALLincr=TRUE,do.plot = TRUE,Year="CustomYear",
                  mymain="Example 2")
#Look at the numeric outputs:
Out3[[2]]
#Plot summaries
graphics.off();x11(width=30,height=25);ggarrange(Out3[[3]][[1]],ggarrange(Out3[[3]][[2]],Out3[[3]][[3]],Out3[[3]][[4]],nrow=1,align="h"),nrow=2)

WQData_Ex3<-Out3[[1]]

##########################################################
# Seasonality Test on Raw Data
##########################################################
#Detemrine the appropriate season for the data (adds on column "Season")
#Conduct a seasonality test
Seas3<-GetSeason(WQData_Ex3,ValuesToUse = "RawValue",printKW = TRUE,mymain="Example 3",do.plot = TRUE)

Seas3[[2]]
#The raw data is seasonal - but biannually, not quarterly, like the time increment
x11();Seas3[[3]]
WQData_Ex3<-Seas3[[1]]
##########################################################
# Perform Trend Tests on Raw Data
##########################################################
Trend3<-SeasonalTrendAnalysis(WQData_Ex3,mymain="Ex 3 Raw Trend",
                         Year="CustomYear",do.plot=T)

Trend3[[1]]
x11();Trend3[[2]]

#########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         EXAMPLE 4: Annual Data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#########################################################
# Preliminary Setup
#########################################################
#Format myDate
WQData_Ex4$myDate <- as.Date(as.character(WQData_Ex4$sdate),"%Y-%m-%d")

#Add on seasons and extra date informaiton - NOTE, using a custom year from Jul-Jun
WQData_Ex4<-GetMoreDateInfo(WQData_Ex4,firstMonth = 7)

#Process censored values  (NOTE THERE AREN'T ACTUALLY ANY CENSORED VALUES IN THIS DATASET)
WQData_Ex4 <- RemoveAlphaDetect(WQData_Ex4,ColToUse="Value")

#########################################################
# Inspect Data
#########################################################
#Make time series plot
Out4<-InspectTrendData(WQData_Ex4, EndYear = 2017, ReturnALLincr=TRUE,do.plot = TRUE,Year="CustomYear",
                       mymain= "Time Series of Annual Data")
#Look at the numeric outputs:
Out4[[2]]
#look at the plots:
graphics.off();x11(width=30,height=25);ggarrange(Out4[[3]][[1]],ggarrange(Out4[[3]][[2]],Out4[[3]][[3]],Out4[[3]][[4]],nrow=1,align="h"),nrow=2)

WQData_Ex4<-Out4[[1]]
##########################################################
# Perform Trend Tests on Raw Data
##########################################################
Trend4<-NonSeasonalTrendAnalysis(WQData_Ex4,mymain="Ex 4 Raw Trend: Annual Data",
                               Year="CustomYear",do.plot=T)
Trend4[[1]];x11();Trend4[[2]]

##########################################################
##########################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Trend aggregation examples
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##########################################################
# Read in Example dataset
##########################################################
#This dataset includes observations from 3 variables over multiple sites
load("RegionExData_v2401.rdata")
#Add on the extra date infromation (removeAlphaDetect has already been done)
RegionExData<-GetMoreDateInfo(RegionExData,firstMonth = 7)

##########################################################
# Analyse trends by site x analyte
##########################################################
#Use wrapper function to undertake trend analysis across all sites and variables
My10yTrends <- ddply(RegionExData, c("siteID", "analyte"), function(x) {
  cat(as.character(x$siteID[1]), as.character(x$analyte[1]), "\n")
  doMyTrends_v2502(x, Year = "CustomYear", ValuesToUse="RawValue",
                   UseMidObs=TRUE,do.plot=FALSE)}, .progress = "win")

##########################################################
#Make stacked bar charts of simple improvement confidence categories
##########################################################
#First add on confidence category
My10yTrends$SimpImpCat<-AssignConfCat(My10yTrends,CatType="Improve",nCat="Simple")


#Make a table of proportions of sites in each confidence category - note, dropping "not analysed"
ConfCats<-ddply(My10yTrends[My10yTrends$SimpImpCat!="Not Analysed",], "analyte", function(x) 100*table(x$SimpImpCat)/sum(table(x$SimpImpCat)), .drop = F)[,-7]
rowSums(ConfCats[,2:6]) #Check to see that each analyte sums to 100%

# stack to make a plot
ConfCatsStack <- stack(ConfCats[, 2:6])
ConfCatsStack$ind <- factor(as.character(ConfCatsStack$ind), levels=names(ConfCats)[2:6])
ConfCatsStack$analyte <- ConfCats$analyte

#Assign a color ramp for the plot
require(RColorBrewer)
myCols <- c(brewer.pal(5,"RdYlBu"))

graphics.off(); x11(height = 7, width = 10)
ggplot(data=ConfCatsStack, aes(x=analyte, y=values, fill=ind)) + xlab("Analyte")+
  geom_bar(stat="identity") + scale_fill_manual(values=myCols,name="Confidence and\ndirection")+
  ylab("Proportion of sites (%)")+coord_cartesian(ylim= c(0,100),expand=FALSE) + theme_bw()



###########################################################
# Evaluate aggregate regional trend direction and strength
##########################################################
#Perform the regional assessment
TauRegional <-  dlply(My10yTrends,.(analyte), function(y) getTAU(x=y, obs=RegionExData[RegionExData$analyte==y$analyte[1],],SiteNameCol="siteID",VarNameCol="analyte",
                                                            ValuesToUse = "RawValue", Year = "CustomYear",
                                                            UseMidObs=T))  # analyse aggregate trends

#Extract summary table
TauRegionalSummary <- ldply(TauRegional, function(x) return(x$SummaryResults))
print(TauRegionalSummary)

#Make plot of correlations between sites, by variable
p<-PlotRegionalCor(TauRegional);x11();p

#make summary plot of regional aggregate trend confidence and direction
p1<-PlotRegionalConfDir(TauRegionalSummary,mymain="10-year regional aggregate trends");x11();p1

#Do the same, but use direction as improving or degrading - need to specify variables
#that should be excluded (i.e., temperature, pH, etc.) and variables where increasing trends
#indicate improvement (i.e., clarity, MCI etc. )
TauRegional_Improve <-  dlply(My10yTrends,.(analyte), function(y) getTAU(x=y, obs=RegionExData[RegionExData$analyte==y$analyte[1],],SiteNameCol="siteID",VarNameCol="analyte",
                                                              ValuesToUse = "RawValue", Year = "CustomYear",
                                                              UseMidObs=T,DirType="Improve") ) # analyse aggregate trends

#Extract summary table
TauRegionalSummary_Improve <- ldply(TauRegional_Improve, function(x) return(x$SummaryResults))
print(TauRegionalSummary_Improve)


p2<-PlotRegionalConfDir(TauRegionalSummary_Improve,mymain="10-year regional aggregate trends",DirType="Improve");x11();p2

##########################################################

