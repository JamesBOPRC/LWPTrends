library(tidyverse)

Raw_NIWA_Data <- read.csv("DF_for_Trends.csv")
Raw_NIWA_Data <- Raw_NIWA_Data %>%
    mutate(Time = parse_date_time(Time, orders = c("%Y-%m-%d %H:%M:%S","%Y-%m-%d"),tz="etc/GMT+12"))



Trend_DF <- Raw_NIWA_Data %>%
  filter(Site %in% c("IG265664","QJ471191")) %>%
  filter(Parameter %in% c("NNN (g/m^3)","VC - BD (m)")) %>%
  mutate(myDate = as.Date(Time, format = "%Y-%m-%d",tz="etc/GMT+12")) %>%
  select(Site, LocationName, myDate, Parameter, Value, Discharge)

Trend_DF <- GetMoreDateInfo(Trend_DF)
Trend_DF <- RemoveAlphaDetect(Trend_DF,ColToUse="Value")


My10yTrends <- ddply(Trend_DF, c("Site", "Parameter"), function(x) {doMyTrends_v2502(x, Year = "Year", propIncrTol=0.9, propYearTol=0.9,TrendPeriod=10, EndYear=2020, ValuesToUse="RawValue", UseMidObs=TRUE, do.plot=FALSE)}, .progress = "win")





debug(doMyTrends_v2502)
Trend_DF %>%
  group_by(Site, Parameter) %>%
  summarise(doMyTrends_v2502(.,Year = "Year", propIncrTol=0.9, propYearTol=0.9,TrendPeriod=10, EndYear=2020, ValuesToUse="RawValue", UseMidObs=TRUE, do.plot=FALSE))


###### Lets stick with our functions.


Trend_DF <- Raw_NIWA_Data %>%
  filter(Site %in% c("IG265664","QJ471191")) %>%
  filter(Parameter %in% c("NNN (g/m^3)","VC - BD (m)")) %>%
 # mutate(myDate = as.Date(Time, format = "%Y-%m-%d",tz="etc/GMT+12")) %>%
  select(Site, LocationName, Time, Parameter, Value, Discharge) %>%
  distinct() %>%
  pivot_wider(id_cols = c(Site,LocationName, Time, Discharge),values_from = Value, names_from = Parameter) %>%
  rename("FLOW" = Discharge, 'NNN' = `NNN (g/m^3)`, "VC" = `VC - BD (m)`)


# Data,
# param_colnames,
# siteID_col,
# time_col = "Time",
# Timeframe="All",
# Seasonal = F,
# Co_variate_Adjust=F,
# Co_variate_col=NA,
# log_parameters = c("ECOLI","TSS","FC","ENT")

source('./R/LWPFunctions.R')

undebug(Fixed_Trend_Analysis)
undebug(ContinuousTS)


Trend_DF %>%
  ggplot()+
  geom_point(aes(x=Time,y=NNN))+
  facet_wrap(~Site)




Output_NN <- Fixed_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),siteID_col = "Site",
                     Seasonal=F,Co_variate_Adjust=F,time_col = "Time")






debug(ImprovementConfCat_LAWA)

ImprovementConfCat_LAWA(Output_NN)


undebug(AssignConfCat)
Output_NN$Direction <- AssignConfCat(Output_NN,CatType="Improve",Reverse=c("VC"))


Output_YN <- Fixed_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),siteID_col = "Site",
                                  Seasonal=T,Co_variate_Adjust=F,time_col = "Time")


Output_NY <- Fixed_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),siteID_col = "Site",
                                  Seasonal=F,Co_variate_Adjust=T,Co_variate_col = "FLOW",time_col = "Time")

Output_YY <- Fixed_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),siteID_col = "Site",
                                  Seasonal=T,Co_variate_Adjust=T,Co_variate_col = "FLOW",time_col = "Time")

#what if trend DF has no covariate?

Trend_DF2 <- Trend_DF %>% select(-FLOW)

Output_NN2 <- Fixed_Trend_Analysis(Data = Trend_DF2, param_colnames = c("VC","NNN"),siteID_col = "Site",
                                  Seasonal=F,Co_variate_Adjust=F,time_col = "Time")

Output_YN2 <- Fixed_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),siteID_col = "Site",
                                  Seasonal=T,Co_variate_Adjust=F,time_col = "Time")


Output_NY2 <- Fixed_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),siteID_col = "Site",
                                  Seasonal=F,Co_variate_Adjust=T,Co_variate_col = "FLOW",time_col = "Time",
                                  Timeframe = 10)

# Data,
# param_colnames,
# siteID_col,
# time_col = "Time",
# Timeframe="All", #or a number that will be converted to years
# Co_variate_col=NA,
# Data_gap = 3,
# min_values = 20,
# log_parameters = c("ECOLI","TSS","FC","ENT")


debug(Automated_Trend_Analysis)

Automated_Trend_Analysis(Data = Trend_DF, param_colnames = c("VC","NNN"),Co_variate_col="FLOW",
                         siteID_col = "Site",Timeframe = 10)



# Suggested procedure for trend analysis.
# 1. Evaluate the data for the time-period of interest. Assess if there is sufficient data using MakeTimeSeries, PlotData
# 2. determine the appropriate time-period increment  e.g. months or quarters and define Year and time increment in the dataframe.
# 3.Based on  experience, select likely candidates for covariates (e.g., Flow).
# 4.Examine the covariate variation with time. If there is a significant trend in the covariate then it's use may cause a trend in the variable.
# 5.Determine if variable is correlated with covariate using Spearman non-parametric correlation. Confirm by plotting variable against covariate and deciding the best form of the relationship (Linear, log-log, or GAM).
# If GAM select the appropriate degrees of freedom to get a good fit and a monotonic relationship.
# If the relationship between the variable and covariate is a straight line with close to zero slope then applyng the covariate correction will have no effect on the result.
# 6.Carry out Seasonal Kendall test with the number of seasons per year equal to the sampling interval and covariate adjustment, if appropriate.

# implement  (NOTE set all values below the highest censored to be censored values. )
# The option to set all values (censored or otherwise) less than highest censoring limit as censored at the
# highest censor limit eliminates the possibility of that these data cause a spurious trend. For example, if
# a series of values were 6, 5 , <4, 3, 4, 2, <2, 5, 2 , 1, a Kendall trend analysis without setting values to
# the maximum would show a trend (Mann-Kendall P=0.06). Increasing all values less than the highest limit of
# 4 to a value of <4, results in a less significant trend (Mann-Kendall P=0.15)

