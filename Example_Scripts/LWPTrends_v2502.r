#######################################################################
#######################################################################
#
#              ********  LWP-TRENDS  *********    
#
#######################################################################
#######################################################################
# Ton Snelder and Caroline Fraser
# LWP Ltd
# Created: September 2017
# Most Recent update: 
#September 2018; changes to calcualation of Sen slopes 
# and probability trend was decreasing 
# November 2019: Changes to the implementation of the "hi-censor" rules
#   allowing ggplot objects to be returned on request
#   including bimonthly as built in a time increment option
#   error trap on seasonality test
# April 2021: Changes in language from Proability to Confidence.  Updates
# to HICensor filter, and use of flow adjustment function with HiCensor filter

#######################################################################e#######################################################################
# Please read accompanying documentation: LWPTrendsHelp_c2503.pdf and
# note the disclaimer therein for the use of these scripts
#
# See also companion example files:
#       RunLWPTrendsExample_v2503.R
#       LWPTrends_ExampleData.rdata
#
#######################################################################
#######################################################################
# Functions for analysis of water quality trends
# These function replicate, as closely as possible, the Mann-Kendall and Seasonal Kendall tests in the TimeTrends software.
# Differences between this code and TT occur due to differences in the handling of censored values when calculating S and its variance.
# This is due to the key function is "GetKendal", which computes the Kendall test for censored data.
# The function is adapted from the "cenken" function in the NADA package
# The function returns the Kendall statistic and variance as described by Helsel (2012) "Statistics for censored environmental data using MINITAB and R", John Wiley & Sons. 


# the main analysis functions are as follows:
    # Impute.upper - provides imputed values for left Censored data
    # Impute.lower - provides imputed values for right censored data
    # RemoveAlphaDetect - removes non detect ("<" and ">") and produces a dataframe with a numeric value and a boolean defining if value is censored. 
    # GetFlowAdjustment - flow adjusts values
    # InspectTrenData - options ot view data as timeseries or a heat plot showing WQ data or censored occasions
    # NonSeasonalTrendAnalysis - wrapper function that calls MannKendall and SenSlope
    # SeasonalTrendAnalysis - wrapper function that calls SeasonalKendall and SeasonalSenSlope
        # MannKendall - Mann-Kendall test
        # SenSlope - Calculates the annual sen slope 
        # SeasonalKendall - Seasonal Kendall test.
        # SeasonalSenSlope - Calculates the seasonal sen slope
#Additonal Trend Aggregation functions  (applied to trend outputs across one variable at many sites)
    # AnalysePd - provides proportion of improving trends (Pd Statistic)
    # FaceValueCounter - counts numbers of increasing and decreasing sites
#Direction
    # ImprovementConfCat - Assigns a cateorical probability of improvment to each trend (based on IPCC categories)
    # ImprovementConfCatLAWA - Assigns a cateorical probability of improvment to each trend (based on LAWA categories)


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

#######################################################################
#Utility Functions
#######################################################################
# if(is.na(match("plyr", installed.packages()[,1]))) stop("these functions require that the R package plyr is installed. Install plyr and try again.")
require(plyr); require(tidyr); require(viridis)
#require(gridExtra)
# if(is.na(match("NADA", installed.packages()[,1]))) stop("these functions require that the R package NADA is installed. Install NADA and try again.")
require(NADA) 
require(lubridate)

# count unique elements in a given vector
countUnique <- function(myVec) {length(unique(myVec))}

# recover numerical value from number stored as factor
unfactor <- function(f) {return(as.numeric(levels(f))[f])} 

unpaste <- function (str, sep = "/", fixed = T) {
  w <- strsplit(str, sep, fixed = fixed)
  w <- matrix(unlist(w), ncol = length(str))
  nr <- nrow(w)
  ans <- vector("list", nr)
  for (j in 1:nr) ans[[j]] <- w[j, ]
  ans
}

# split string (similar to unpaste)
DougSplit  <- function(mytext, myPart = 1, ...) {  # use: sapply(MyNames, DougSplit, split = "\\.", myPart = 2)
  myOut <- strsplit(mytext, ...)[[1]][myPart]      # note space is DougSplit("1996-02-07 00:00:00", split = "\\s")
  return(myOut)
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

#######################################################################
#Imputation Functions
######################################################################
Impute.lower <- function(x, ValuesToUse = "RawValue", forwardT="log", reverseT="exp", do.plot=F,minObsforfit=10) {     
  # this uses ROS to impute values for obs below multiple detection limits. There is no longer any randomisation of the imputed values 
  
  # x is a data-frame containing (at a minimum)
  # x$RawValue : the raw data 
  # x$CenType            :a flag indicating whether the raw data are left censored ("lt"), right censored ("gt"), 
  #                     or within the lower and upper detection limits ("ok")
  # x$myDate: the date of the samples (to ensure correct ordering)
  # forwardT: forward transformation to be used in the NADA:ros() function
  # reverseT: reverse transformation to be used in the NADA:ros() function
  
  myData <- x
  myData$CenType <- as.character(myData$CenType)
  
  if(sum(myData$CenType == "lt") == 0 | sum(myData$CenType =="lt") == nrow(myData)) {       # if no left censored obs, OR all are censored, return the values.
    
    if(sum("lt"==myData$CenType) == 0 ) { # first deal with no censored obs
      print("no observations are left censored")
      myData$i1Values <- myData[, ValuesToUse]
      myData$LeftImpute <- "No censored - no imputation required"
    }
    if(sum(myData$CenType =="lt") == nrow(myData)) {
      print("all observations are left censored - cannot fit model")
      myData$i1Values <- myData[, ValuesToUse]/2
      myData$LeftImpute <- "All censored - cannot impute"
    }
    
  } else {
    
    r.na2 <- is.na(myData[, ValuesToUse])  # Treat NAs as censored values, then sub NA back  at end
    # Since they are marked as censored, they will not influence the ROS regression.  
    # They will receive imputed values, but those will later be overwritten with NAs
    CenType.original <- myData$CenType
    myData$CenType[r.na2] <- "lt"
    myData <- myData[order(myData$myDate),]
    rownames(myData) <- 1:nrow(myData) # to keep order of values in time
    print("some below limits")
    
    ######################   
    # catch some bad myData if censored values exceed max of uncensored values
    #cat("Catch_", length(myData[myData$CenType != "lt", "RawValue"]), "\n")
    if(sum(myData$CenType =="lt") != length(myData$CenType)) MaxUncen <- max(myData[myData$CenType != "lt", "RawValue"]) 
    
    
    if( sum(myData[myData$CenType == "lt", ValuesToUse] >=  MaxUncen, na.rm=T)>0 ) { # if any censored values exceed the observed uncensored values
      print("some bad data")
      BadDates <- (myData$CenType == "lt" & myData[, ValuesToUse] >=  MaxUncen) 
      myData_RawValue_orig <- myData[, ValuesToUse]
      myData$RawValue[BadDates] <- 0.5 * MaxUncen
      myData$CenType[BadDates] <- "lt" # this is already the case  if condition above is True
      
      data2 <- myData  
      #################
      
    } else {
      print("no bad dates")
      data2 <- myData
    }
    rownames(data2) <- 1:nrow(data2) # to keep order of values in time
    usedValues <- data2[, ValuesToUse]
    names(usedValues) <- rownames(data2)    
    
    
    usedValues[usedValues <= 0] <- min(data2[data2$CenType == "lt", ValuesToUse], na.rm=T) # replace any zero values with the minimum DL
    
    # added a catch as ros can fail occasionally if the data are not conducive to a regression
    myros <- try(ros(obs = usedValues, censored = data2$CenType == "lt",forwardT=forwardT,reverseT=reverseT))
    
    if (any(class(myros) == c("ros", "lm"))&sum(data2$Censored==FALSE)>minObsforfit) {  # if a model is fitted then contiune 
      
      SortValues <- sort(usedValues, na.last = T)  # sort order of values # na.last to preserve the NA values on sorting
      
      NewROS <- cbind.data.frame(SortValues,as.data.frame(myros))
      
      # the plot helps to show why imputed replacements for censored values are often greater than
      # non-censored values. 
      if(do.plot)  {
        #graphics.off();x11(); 
        plot(myros)
        with(NewROS[!NewROS$censored,],  points(x=qnorm(pp), y=obs, pch=16, col="red"))
        with(NewROS[NewROS$censored,],  points(x=qnorm(pp), y=modeled, cex=2, pch=16, col="blue"))
        legend("topleft",  legend=c("Uncensored","Censored"), col = c("red","blue"), text.col = "black",  pch = c(16, 16), bg = 'white', inset = .05)  
      }
      
      NewROS$OrigOrder <- as.numeric(names(SortValues)) # this is the original time order
      SortROS <- NewROS[order(NewROS$OrigOrder), ]  # reorder ROS to original time order
      
      data2$i1Values[SortROS$OrigOrder] <- SortROS[,"modeled"] # retain the order.
      data2$i1Values[r.na2] <- NA
      data2$CenType[r.na2] <- CenType.original[r.na2]
      data2$LeftImpute <- "Imputed"
      # note that any data where CenType=="lt" but (original) converted_value > MaxUncen now has an imputed value from ROS
    } else {
      data2$i1Values <- myData$RawValue
      data2$i1Values[data2$CenType=="lt"] <- myData$RawValue[data2$CenType=="lt"]/2
      data2$LeftImpute <- "Not Imputed - model fit failed"
      print("Warning: ROS model was not fitted. Imputed values are original values")
    }
    myData <- data2
  }
  return(myData)
}

Impute.upper <- function(x, ValuesToUse = "i1Values") {   
  # this function uses survreg (package::survival) to impute values right censored observations - i.e., ABOVE (multiple) detection limits

  myData <-  x
  myData<-myData[order(myData$myDate),]
  
  if(sum(myData$CenType == "gt") == 0 | sum(myData$CenType =="gt") == nrow(myData)) { # if no right censored obs, OR all are censored, return the values.
    
    if(sum("gt"==myData$CenType) == 0 ) { # first deal with no censored obs
      print("no observations are right censored")
      myData$i2Values <- myData[, ValuesToUse]
      myData$RightImpute <- "No censored - no imputation required"
    }
      
      if(sum(myData$CenType =="gt") == nrow(myData)) {
        print("all observations are right censored - cannot fit model")
        myData$i2Values <- myData[, ValuesToUse]
        myData$RightImpute <- "All censored - cannot impute"
      }
      
    } else { # otherwise impute
      
  
# a survreg model cannot be fit with n<24; in this case, simply set all right censored to 1.1*RawValue  
 if(nrow(myData) < 24 | sum(myData$CenType == "gt")==0 ) { 
    myData$i2Values <- myData[, ValuesToUse]
    myData$i2Values[myData$CenType == "gt"] <- 1.1*myData$i1Values[myData$CenType == "gt"] # add 10%
    myData$RightImpute <- "Increased by 10%"
  } else {
    myData$Temp <- myData[, ValuesToUse]
    myData$Temp[myData[, ValuesToUse]==0]<-min(myData[, ValuesToUse][myData[, ValuesToUse]!=0])  #This is  a catch for bad data that makes the weibull model fall over
    # fit distribution
    myData$status <- myData$CenType != "gt"   # note well if myData flag is "gt" (meaning censored) this is equivalent to "not dead" status in a survival model and therefore is
    myMod <- survreg(Surv(Temp, status) ~ 1, data = myData, dist="weibull") # using original observed non-censored
    RScale <- as.numeric(exp(myMod$coefficients)) # this is the rweibull scale (see "survreg" function notes)
    RShape <- 1/myMod$scale   # this is the weibull shape
    
    SynthDist <- sort(rweibull(10000, shape=RShape, scale=RScale))   # synthetic distribution
    # NB large sample taken to ensure that there are several values in SynthDist > the non detects
    # otherwise the function below falls over. 
    # Include a catch here in the cases that the non-detects are much larger than suggested by the distrbution,
    # replace RawValues with 99th percentil of the distribution
    
    myData$Temp[myData$CenType=="gt" & myData$Temp > quantile(SynthDist,0.995)]<-quantile(SynthDist,0.995)
    
    myData$i2Values<-myData$Temp
    myData$i2Values[myData$CenType=="gt"]<-sapply(myData$Temp[myData$CenType=="gt"],function(x) resample(SynthDist[SynthDist > x], size=1))
    myData$status=NULL
    
    myData <- myData[, -which(names(myData) == "Temp")] # drop the Temp column
    myData$RightImpute <- "Imputed"
       }
    }
  return(myData)
}


#######################################################################
# Pre-Process censored data
#######################################################################
RemoveAlphaDetect <-  function(Data,ColToUse="Value") {   #  removes non detect ("<" and ">") and produces a dataframe with a numeric value and a boolean defining if value is censored.
  x<-Data[,ColToUse]
  if(is.numeric(x)) xNumeric <- x   # check values may be already numeric
  if(is.factor(x)) xNumeric <- unfactor(x) # converts to numeric (but note the values < and > are returned as NA)
  if(is.character(x)) xNumeric <- as.numeric(x) # converts to numeric (but note the values < and > are returned as NA)
  isNonDetLT <- grep("<", x)  # which values are below detection
  isNonDetGT <- grep(">", x)  # which values are above detection
  DetectionLimLT <- as.numeric(sapply(as.character(x[isNonDetLT]), function(x) DougSplit(x, myPart=2, split= "<", fixed = T ))) # returns the numeric values
  DetectionLimGT <- as.numeric(sapply(as.character(x[isNonDetGT]), function(x) DougSplit(x, myPart=2, split= ">", fixed = T ))) # returns the numeric values
  xNumeric[isNonDetLT] <- DetectionLimLT 
  xNumeric[isNonDetGT] <- DetectionLimGT 
  Censored <- rep(FALSE, times = length(x))
  Censored[isNonDetLT] <- TRUE
  Censored[isNonDetGT] <- TRUE
  CenType <-  rep("not", times = length(x))  # classification of the type of censoring 
  CenType[isNonDetLT] <- "lt" # less than censored
  CenType[isNonDetGT] <- "gt" # greater than censored
  CenType <- factor(CenType, levels = c("gt", "lt", "not"))
  # cat("*** It is safe to ignore warning message ***")
  
  Censored <- ifelse(Censored == "TRUE", T, F) # ensure this a binary
  Data$RawValue<-xNumeric
  Data$Censored<-Censored
  Data$CenType<-CenType
  return(Data) # raw values means NOT Flow-Adjusted and no < or > values 
}

#######################################################################
# Add extra date information onto original dataset
#######################################################################
#This function adds on some extra date formats onto the dataframe for subsequent analysis
#Including a default time increment
GetMoreDateInfo<-function(Data,firstMonth=1,FindDateShifts=TRUE){
  #First month is the first month of the custom year.  1 is defaul (Jan-Dec year)
 
  
  #First, assign a temporary date to observations, if FindDateShifts=TRUE that shifts observations to teh next month in the case tat they are at the end of the month
  # and there was another observation at the start of the month, and there is no observation in the next month.
  #If false, we just use myDate
  ifelse(FindDateShifts,Data<-ShiftTempDate(Data),Data$tempDate<-Data$myDate)

  if(firstMonth>1){
    yearshift<-365.25-as.numeric(strftime(as.Date(paste("2000",firstMonth,"1",sep="-")), format = "%j"))
    MyNumMonthString<-c(NumMonthString[firstMonth:12],NumMonthString[1:firstMonth-1])
  }else{
    MyNumMonthString<-NumMonthString
    yearshift<-0
  }  
  
  Data$Year <- as.numeric(format(Data$tempDate, "%Y"))
  if(firstMonth!=1){# Customised  year
    Data$CustomYear<-Data$Year
    Data$CustomYear[as.numeric(format(Data$tempDate,"%m"))>=firstMonth]<-Data$Year[as.numeric(format(Data$tempDate,"%m"))>=firstMonth]+1
  }

  Data$Month <- format(Data$tempDate, "%b")    # abbreviated  months
  Data$Month <- factor(Data$Month, levels = MyNumMonthString) 
 
  #Make a data frame of all the extra time increment option s
  A<-data.frame(Month=factor(MyNumMonthString,levels=MyNumMonthString))
  A$Num<-1:12

  A$BiMonth<-paste(A$Month[ceiling((A$Num)/2)*2-1],A$Month[ceiling((A$Num)/2)*2],sep=".to.")
  A$Qtr<-paste(A$Month[ceiling((A$Num)/3)*3-2],A$Month[ceiling((A$Num)/3)*3],sep=".to.")
  A$BiAnn<-paste(A$Month[ceiling((A$Num)/6)*6-5],A$Month[ceiling((A$Num)/6)*6],sep=".to.")

  A$BiMonth<-factor(A$BiMonth,levels=unique(A$BiMonth))
  A$Qtr<-factor(A$Qtr,levels=unique(A$Qtr))
  A$BiAnn<-factor(A$BiAnn,levels=unique(A$BiAnn))
  
  A$Num<-NULL
  
  #Determine the middle day for each time increment
 MidTimeIncrList<-list()
  for(i in 1:12){MidTimeIncrList[[i]] <- seq(365/i/2,365-365/i/2,by=365/i)-yearshift-1}
 #Commit to the gloabl environment for later use
 assign("MidTimeIncrList",MidTimeIncrList,envir = .GlobalEnv)

  #Glue all time increment options onto original dataframe
Data<-merge(Data,A)
Data$tempDate<-NULL

return(Data)}



###############################################################################
# This function moves dates that are really belonging to the sampling of the next month
###############################################################################

ShiftTempDate<-function(x){
  mynames<-names(x)
  x<-x[order(x$myDate),]
  x$dT_bwd<-as.numeric(c(0,diff(x$myDate)))
  x$dT_fwd<-as.numeric(c(diff(x$myDate),0))
  x$Month<-as.numeric(format(x$myDate, "%m"))
  x$Year<-as.numeric(format(x$myDate, "%Y"))
  x$Day<-as.numeric(format(x$myDate, "%d"))
  x$dMonB<-c(0,diff(x$Month))
  x$dMonB[x$dMonB<0]<-x$dMonB[x$dMonB<0]+12
  x$dMonF<-c(diff(x$Month),0)
  x$dMonF[x$dMonF<0]<-x$dMonF[x$dMonF<0]+12
  
  #Index of observations that should be moved forward a month
  indF<-x$dMonB==0&x$dMonF>1&x$dT_bwd>20&x$dT_fwd>20
  #Index of observations that should be moved backward a month
  indB<-x$dMonF==0&x$dMonB>1&x$dT_bwd>20&x$dT_fwd>20
  
  x$NextMonth<-x$Month+1
  x$NextMonth[x$NextMonth==13]<-1
  x$PrevMonth<-x$Month-1
  x$PrevMonth[x$NextMonth==1]<-12
  
  x$tempDateF<-make_date(year=x$Year,month=x$NextMonth,day=1)

  x$tempDateB<-make_date(year=x$Year,month=x$PrevMonth,day=1)

  x$tempDate<-x$myDate
  x$tempDate[indF==TRUE]<-x$tempDateF[indF==T]
  x$tempDate[indB==TRUE]<-x$tempDateB[indB==T]
  
  x<-x[,c(mynames,"tempDate")]
  return(x)
}



#######################################################################
# Define time increment strings
#######################################################################
NumMonthString<-format(as.Date(sapply(seq(1,12),function(x) paste0("2000-",x,"-01"))),"%b")

WhatistheIncr<-function(x,val="TimeIncr"){
  #This is a function to determine what is the time increment type, just to make nice titles in plots
  # NtheTimeIncr<-length(levels(x[,val]))
  ind<-setdiff(1:ncol(x),which(names(x) %in% c("TimeIncr","Season")))
  mycol<-names(x)[ind[sapply(ind,function(y) identical(as.character(x[,val]),as.character(x[,y])))]]
  
  if(mycol %in% c("Year","CustomYear")){
    myTimeIncr="Annual"
  }else{
    myTimeIncr=SeasonLabs$mySeas[SeasonLabs$SeasName==mycol]
  }
  return(myTimeIncr)
}
#######################################################################
# Covariate adjustment 
#######################################################################
# fit a lowess model to log-log and where Ci and Qi are concentration and Covariate time series
# method options are "Gam" or "Loess"
# the argument Span < 1, the neighbourhood includes proportion "span" of the points, and these have tricubic weighting (proportional to (1 - (dist/maxdist)^3)^3)
# this is evidently NOT the same as the "% of points to fit" parameter in TimeTrends.
if(is.na(match("gam", installed.packages()[,1]))) stop("these functions require that the R package gam is installed. Install gam and try again.")
require(gam)

AdjustValues <- function(myFrame, method = c("Gam", "LogLog", "LOESS"), ValuesToAdjust = "RawValue", 
                         Covariate = "Flow", Span = c(0.5, 0.75), 
                         do.plot = F, plotpval=T,plotloglog=F,
                         Covariate_Lab="Flow",mymain=NA,HiCensor=FALSE,...) {  # myFrame <- WQData[WQData$siteID == "AK1" & WQData$analyte == "TN", ]

  # make sure the data are in date order!
  myFrame <- myFrame[order(myFrame$myDate), ]
  row.names(myFrame) <- 1:nrow(myFrame)        # deliberately naming the rows in ascending order so that adjusted values can be matched to the orginal values at the end
  
  OriginalValues <- myFrame[,ValuesToAdjust]  #
  
  if (HiCensor!=FALSE){  #This makes a hicenosr adjustment just for the outputs, we still use all data for the fitting
    x<-ApplyHiCensor(myFrame,ValuesToAdjust,HiCensor)
    OriginalValues[x[,'CenType']=="lt"]<-OriginalValues[x[,'CenType']=="lt"]/2
    OriginalValues[x[,'CenType']=="gt"]<-OriginalValues[x[,'CenType']=="gt"]*1.1
    
  }else{
    OriginalValues[myFrame[,'CenType']=="lt"]<-OriginalValues[myFrame[,'CenType']=="lt"]/2
    OriginalValues[myFrame[,'CenType']=="gt"]<-OriginalValues[myFrame[,'CenType']=="gt"]*1.1
  }
  
  # number of analyses and different adjustments output  
  
  if (any(method == "LOESS")) {
    NoAdjusts <- (length(method)-1) + length(Span)
    AdjustNames <- c(method[method != "LOESS"], paste0("LOESS", Span))
  } else {
    NoAdjusts <- length(method)
    AdjustNames <-method
  }
  
  AdjustList <- vector("list", NoAdjusts) # maintain a list of the adjustments made using different models
  names(AdjustList) <- AdjustNames
  
  # dataframe for output adjusted values
  OUT <- data.frame(matrix(nrow=length(OriginalValues), ncol=NoAdjusts*3))
  names(OUT) <- c(AdjustNames,sapply(AdjustNames,function(x) paste(x,"R",sep="_")),sapply(AdjustNames,function(x) paste(x,"p",sep="_")))
  OUT$myDate <- myFrame$myDate
  OUT$OriginalRawValues <- myFrame[,ValuesToAdjust]
  
  
  #For censored values, use half of <DL and 1.1 >AL
  myFrame[myFrame[,'CenType']=="lt",ValuesToAdjust]<-  myFrame[myFrame[,'CenType']=="lt",ValuesToAdjust]/2
  myFrame[myFrame[,'CenType']=="gt",ValuesToAdjust]<-  myFrame[myFrame[,'CenType']=="gt",ValuesToAdjust]*1.1
  myFrame$ValueswithHiCen<-OriginalValues
  #OUT$analyte<-myFrame[,"analyte"]
  #OUT$siteID<-myFrame[,"siteID"]
  myplot<-NA
  # what percentage of data has Covariates?
  myFrame$HasFlowAndData <- !is.na(myFrame[,Covariate]) & !is.na(OriginalValues) # rows with Covariate and values. 
  PropWithFlow <- sum(myFrame$HasFlowAndData)/nrow(myFrame)
  #
  #
  if(all(is.na(myFrame[,Covariate])) |  (nrow(myFrame) < 10) | sum(myFrame$HasFlowAndData)<10 |  # no Covariate data, not much data - PropWithFlow < 0.8 |
     length(unique(myFrame[, ValuesToAdjust]))/nrow(myFrame) < 0.05 |   # many values are the same 
     max(unlist(rle(as.vector(myFrame[,Covariate]))[1]))/nrow(myFrame) > 0.3) # or long sequence of indentical Covariate
  {     # dont even try Covariate adjustment
  } else {
    myFrame2 <- myFrame[myFrame$HasFlowAndData, ]    # makes a frame with no NA values
    myFrame2$XXValues <- myFrame2[,ValuesToAdjust]
    myFrame2$XXFlow <- myFrame2[,Covariate]
    
    
    hasZeroValues <- ifelse(sum(myFrame2$XXValues <= 0) >0 | sum(myFrame2$XXFlow <= 0)>0, T, F)  # are there any zero or negative Covariate or variable values
    ZerosHere <- myFrame2$XXValues == 0 | myFrame2$XXFlow == 0                                                       # keep an index of zero values
    myFrame2$XXValues[myFrame2$XXValues == 0] <-  min(myFrame2$XXValues[myFrame2$XXValues > 0])/2                                         # make zeros equal to small values 
    myFrame2$XXFlow[myFrame2$XXFlow <= 0] <- min(myFrame2$XXFlow[myFrame2$XXFlow > 0])/2                                         # make zeros equal to small values 
    myx<-seq(from=min(myFrame2$XXFlow),to=max(myFrame2$XXFlow),length.out=100) #Xvalues for plotting
    
    # fit GAM
    if(any(method == "Gam")) {
      myGam <- gam::gam(XXValues ~ s(XXFlow), data = myFrame2,contrasts = list(a = "contr.sum"))
      concHat <- predict(myGam, newdata = myFrame2)
      Residuals <- myFrame2$XXValues - concHat
      FlowAdjvalues <- median(OriginalValues, na.rm=T) + myFrame2$ValueswithHiCen - concHat 
      R2 <- round(100* (1 - (sum(residuals(myGam)^2)/sum((mean((myFrame2$XXValues)) - (myFrame2$XXValues))^2))))
      R<-cor.test(concHat,myFrame2$XXValues,method="pearson")$estimate
      p<-cor.test(concHat,myFrame2$XXValues,method="pearson")$p.value
      myy<-predict(myGam, newdata = data.frame(XXFlow=myx))
      AdjustList[["Gam"]] <- list(FlowAdjvalues = FlowAdjvalues, Residuals = Residuals, concHat = concHat, Fitted = fitted(myGam), R2 = R2,R=R,p=p,method="Gam",myy=myy,myx=myx)
      # print("**** Safe to ignore MODEL.MATRIX warning *****")
      #Apparenlty this is a long term issue with gam, that now throughs a warning
      # in R v 3.6.0+ but does not affect our evaluations https://stackoverflow.com/questions/57664927/warning-in-gam-with-release-of-r-3-6-1
    }
    
    # fit log-log
    if(any(method == "LogLog")) {
      myLogLog <- lm(log10(XXValues) ~ log10(XXFlow), data = myFrame2)
      concHat <- 10^predict(myLogLog, newdata = myFrame2)
      Residuals <- myFrame2$XXValues - concHat
      FlowAdjvalues <- median(OriginalValues) + myFrame2$ValueswithHiCen - concHat 
      R2 <- round(100 * (1 - (sum(residuals(myLogLog)^2)/sum((mean(log10(myFrame2$XXValues)) - log10(myFrame2$XXValues))^2))))  #Note =- this R2 is in the log-log space
      R<-cor.test(concHat,myFrame2$XXValues,method="pearson")$estimate
      p<-cor.test(concHat,myFrame2$XXValues,method="pearson")$p.value
      myy<-10^predict(myLogLog, newdata =  data.frame(XXFlow=myx))
      AdjustList[["LogLog"]] <- list(FlowAdjvalues = FlowAdjvalues, Residuals = Residuals, concHat = concHat, Fitted = 10^fitted(myLogLog), R2 = R2,R=R,p=p,method="LogLog",myy=myy,myx=myx)  # NB back transformation here (should I cotrrect  for bias???) 
    }
    
    #loess 
    if(any(method == "LOESS")) {
      LoessMods <- vector("list", length=length(Span))          # list to store each loess model
      for(i in 1:length(Span)) {                             # for each value of Span. 
        thisName <- paste0("LOESS", Span[i])
        op <- options(warn=2)
        tt <- try(loess((XXValues)~(XXFlow), span=Span[i], data = myFrame2))
        
        if(is(tt,"try-error"))  {
          LoessMods[[i]] <- list(FlowAdjvalues = rep(NA, nrow(myFrame2)), Residuals = rep(NA, nrow(myFrame2)), concHat = rep(NA, nrow(myFrame2)), R2 = NA,R=NA,p=NA,method=thisName,myy=NA,myx=myx)
        } else {
          op <- options(warn=0) # set back to default warnings.
          Lowess <- loess((XXValues)~(XXFlow), span=Span[i], data = myFrame2)
          concHat <- predict(Lowess, newdata = myFrame2)  # the estimated concentrations for each days Covariate
          Residuals <- myFrame2$XXValues - concHat
          FlowAdjvalues <- median(OriginalValues) + myFrame2$ValueswithHiCen - concHat
          R2 <- round(100 * (1 - (sum(residuals(Lowess)^2)/sum((mean((myFrame2$XXValues)) - (myFrame2$XXValues))^2))))
          R<-cor.test(concHat,myFrame2$XXValues,method="pearson")$estimate
          p<-cor.test(concHat,myFrame2$XXValues,method="pearson")$p.value
          myy<-predict(Lowess, newdata = data.frame(XXFlow=myx))
          AdjustList[[thisName]] <- list(FlowAdjvalues = FlowAdjvalues, Residuals = Residuals, concHat = concHat, Fitted = fitted(Lowess), R2 = R2,R=R,p=p,method=thisName,myy=myy,myx=myx)
        }
      }
    }    
    
    
    if(do.plot!=FALSE) {
      
      
      myFits<-data.frame(do.call("rbind", lapply(AdjustList, function(df) cbind.data.frame(df$myx,df$myy,df$method))))
      names(myFits)<-c("myx","myy","Method")
      myFits$myx<-as.numeric(as.character(myFits$myx))
      myFits$myy<-as.numeric(as.character(myFits$myy))
      myFits$Method<-factor(myFits$Method,levels=AdjustNames)
      if(plotpval) { 
        ModR2_p <- ldply(AdjustList, function(x) {
          PVal <- "NOT SIGNIFICANT p>0.05"
          if(x$p < 0.05) PVal <- "SIGNIFICANT p<0.05"
          return(paste0("R2 =", round((x$R)^2*100), "% ", PVal))  #Note - R2 here is calcualted in teh original unit space
        })
        
        ModR2_p$lab <- paste(ModR2_p[,1], ModR2_p[,2])}else{ModR2_p=NA}
      
      if(do.plot==TRUE){
        myplot<-AdjustValues_GGPlots(myFits=myFits,myFrame2=myFrame2,
                                     ZerosHere=ZerosHere,plotpval=plotpval,
                                     ModR2_p=ModR2_p,
                                     Covariate=Covariate_Lab,mymain=mymain,...)
      }}
    
    for(i in 1:NoAdjusts) OUT[,i][match(names(AdjustList[[i]]$Residuals), row.names(myFrame))] <- AdjustList[[i]]$Residuals#[order(myFrame$myDate)[myFrame$HasFlowAndData]] # put in original date order only for oginal data with value and flow. 
    for(i in 1:NoAdjusts) OUT[,NoAdjusts+i][myFrame$HasFlowAndData] <- AdjustList[[i]]$R
    for(i in 1:NoAdjusts) OUT[,NoAdjusts*2+i][myFrame$HasFlowAndData] <- AdjustList[[i]]$p 

    op <- options(warn=0) # set back to default warnings.
  }  
  if(do.plot==TRUE){OUT<-list(OUT,myplot)}
  #browser()
  return(OUT)  # pass back the frame with all data including the Covariate adjustment
}

AdjustValues_GGPlots=function(myFits=NA,myFrame2=NA,ZerosHere=NA,
                              plotpval=FALSE,ModR2_p=NA,
                              Covariate='Flow',mymain=NA){

  myplot<-ggplot(myFits,aes(x=myx,y=myy,colour=Method))+geom_line(linewidth=1.3)+
                 theme_bw()+geom_point(data=myFrame2,aes(y=XXValues,x=XXFlow),colour="black")+
                theme(legend.position = "inside",
                  legend.position.inside = c(0.05,0.95),
                      legend.background = element_rect(colour="black",fill="white"),
                      legend.justification = c(0, 1))+
    xlab(Covariate)+ylab("Value")
  
  if(is.na(mymain)){myplot<-myplot+ggtitle("Covariation plot") }else{
    myplot<-myplot+ggtitle(paste0(mymain),"\n Covariation plot")}
  
  if(plotpval==TRUE){
    myplot<-myplot+scale_color_discrete(labels=ModR2_p$lab)
  }
  return(myplot)
}
#######################################################################
# Inspect the data. 
#######################################################################

InspectTrendData <- function(x,  
                        Year = "Year", 
                        TrendPeriod = NA, EndYear = NA, propIncrTol=0.9,propYearTol=0.9,
                        ReturnALLincr=FALSE,AnnualOK=c("MCI"),TimeIncrOpts=NA,
                        do.plot = FALSE, UseMidObs=TRUE,FlowtoUse=NA,
                        ...) { # ... arguments passed to the plots 
  
  #The purpose of this function is to (1) prepare the data for subsequent trend analysis, but cutting down to the
  # specified trend period, and determining the time increment for analysis 
  # And (2) optionally to provide metadata about data availability for ALL possible time increments and 
  # (3) optionally to make plots to demonstrate the data availability

  #When ReturnALLincr is FALSE, only the original data set is returned with a new column "TimeIncr" specifying the 
  # highest frequency time increment that meets the minimum data requirements given by propIncrtol (the proportion of
  # time increments that must have at least one observation) and propYearTol (the proportion of time "years" in the trend
  # period that must have at least one observation)
  
  
  
  if(!is.na(FlowtoUse)){
    x$Flow<-x[,FlowtoUse]
  }
  
  #If the start and end year are not provided, automatically extract from dataset
  if (is.na(EndYear)) EndYear<-max(x[,Year])
  if (is.na(TrendPeriod)) TrendPeriod<-EndYear-min(x[,Year])+1
  #Get Date range:
  if(levels(x$Month)[1]=="Jan"){
    EndDate<-as.Date(paste0(EndYear+1,"-01-01"))-1
    StartDate<-as.Date(paste0(EndYear-TrendPeriod+1,"-01-01"))
  }else{

    EndDate<-as.Date(paste0(EndYear,"-",levels(x$Month)[1],"-01"),format="%Y-%b-%d")-1
    StartDate<-as.Date(paste0(EndYear-TrendPeriod,"-",levels(x$Month)[1],"-01"),format="%Y-%b-%d")
  }

  x <- x[x[["myDate"]] >= StartDate & x[["myDate"]] <= EndDate, ] # trims down the data to the year range
  x<-x[!is.na(x$RawValue),] #Get rid of any rows where the observation is NA
  if(nrow(x)>0){
  
#Make data availablity summary table
    if(is.character(TimeIncrOpts)){checkSeas<-SeasonLabs$SeasName[SeasonLabs$mySeas %in% TimeIncrOpts]
    }else{
    checkSeas<-SeasonLabs$SeasName}
    if("analyte" %in% colnames(x)){if(x$analyte[1] %in% AnnualOK){ checkSeas<-c(checkSeas,Year)}}
    
    
  DataAvailability<-ldply(checkSeas,function(y,Data=x){
    Data$IncYear<-paste(Data[[y]], Data[[Year]], sep = "-")
    out<-data.frame(Incr=y,
                    TrendPeriodL=TrendPeriod,
                    nobs=nrow(Data),
                    nYear=length(unique(Data[[Year]])),
                    propYear=length(unique(Data[[Year]]))/TrendPeriod,
      nIncrYear=length(unique(Data$IncYear)),
    propIncrYear=length(unique(Data$IncYear))/(TrendPeriod*length(levels(Data[,y]))))
    if(y==Year){out$propIncrYear=out$propYear}
    return(out)
  })
  
  DataAvailability$propCen <- sum(x$Censored[!is.na(x$RawValue)])/nrow(x)
  DataAvailability$nCenLevelsLT<-length(unique(x$RawValue[x$CenType=="lt"]))
  DataAvailability$nCenLevelsGT<-length(unique(x$RawValue[x$CenType=="gt"]))
  if(!is.null(x$Flow)) {
    DataAvailability$nFlow <- sum(!is.na(x$Flow[!is.na(x[,"RawValue"])]))
  } else {  DataAvailability$nFlow = 0 }
  
  DataAvailability$DataOK<-FALSE
  DataAvailability$DataOK[DataAvailability$propYear>=propYearTol&
                            DataAvailability$propIncrYear>=propIncrTol]<-TRUE
  
  #Select the highest frequency increment that meets the data requirements
  if(any(DataAvailability$DataOK)){
  mytimeincr<-DataAvailability$Incr[which(DataAvailability$DataOK==TRUE)[1]]
  x$TimeIncr<-x[,mytimeincr] #Assign as time increment for the analysis
  }else{ #If none are ok, then we leave TimeIncr as a note ot remove from analysis
    x$TimeIncr<-"none"
  }
  if (ReturnALLincr==TRUE){
    output<-list(x,DataAvailability)
   
  if(do.plot==TRUE){
  myYears <-  min(x[,"Year"]):max(x[,"Year"]) 
  #
  YearsInPeriod <- length(unique(x[,Year]))

  myTimeIncrYearFrame <- expand.grid(list(TimeIncr = NumMonthString, Year = myYears))
  myTimeIncrYearFrame$TimeIncrYear <- paste(myTimeIncrYearFrame$TimeIncr, myTimeIncrYearFrame$Year, sep = "-")
  myTimeIncrYearFrame$TimeIncrYear <- factor(myTimeIncrYearFrame$TimeIncrYear, levels = myTimeIncrYearFrame$TimeIncrYear)
  
  x$TimeIncrYear <- paste(x[["Month"]], x[["Year"]], sep = "-")
  x$TimeIncrYear <- factor(x$TimeIncrYear, levels = myTimeIncrYearFrame$TimeIncrYear)
  x$TimeIncr<-x$Month
  # this takes the median or mid-point of observations in mont (so a month only has one value)
  # if the median of the raw values is less than OR EQUAL to the max censored value, call the value censored, otherwise NOT censored.
  nValperTimeIncr<-ddply(x,.(TimeIncrYear),function(x) data.frame(nVal=nrow(x)))
  
  Data<-ValueForTimeIncr(x, ValuesToUse = "RawValue",Year="Year",UseMidObs)
  Data<-merge(Data,nValperTimeIncr)
  
  # this determines censor type for Data  based on CenType that has highest prevalence in the time increment
  Data$CenType <- ddply(x, "TimeIncrYear", function(y) names(which.max(table(y$CenType))))[,2]
  myTimeIncrYearFrame$Row <- 1:nrow(myTimeIncrYearFrame)
  PlotData <- merge(myTimeIncrYearFrame, Data, all=TRUE) # this over-writes to include all years and TimeIncrs specified by StartYear & EndYear
  PlotData<-PlotData[order(PlotData$Row),]

    output[[3]]<-InspectData_GGPlots(PlotData=PlotData,DateRange=c(StartDate,EndDate),Year="Year",...)#mymain=mymain,Year="Year",
  }
  }else{
    output<-x
  }
  return(output)
  }else{return(NULL)}
  
}




InspectData_GGPlots<-function(PlotData=NA,mymain=NA,Year='Year',DateRange=NA){
  #This produves 3 ggplot plots showing (1) a timeseries of data (2) a matri of values
  # and (3) a matrix of censoring
  #PlotData$theYear<-PlotData[,Year]
  #Only make plots if ggplot is loaded
  if("ggplot2" %in% (.packages())){
    PlotData$Censored<-factor(PlotData$Censored,levels=c("FALSE","TRUE"))
    plotlist=list()
    myTimeIncr<-WhatistheIncr(PlotData[!is.na(PlotData$V1),],val="TimeIncr")
    
    #1. Make plot of timeseries ##########
      mytitle_1<-paste0("Time Series of data")
    if(!is.na(mymain)){mytitle_1<-paste0(mymain,": ",mytitle_1)}

    myplot_1<-ggplot(data=PlotData[!is.na(PlotData$V1),], aes(x=NewDate,y=V1,colour=Censored,shape=Censored))+
      geom_point(show.legend = T)+theme_bw()+ylab('Value')+xlab('Date')+
      scale_colour_manual(values=c("FALSE"="black","TRUE"="red"),
                          labels=c("Observations","Censored"),
                          limits=c("FALSE","TRUE"),na.value="white", drop = FALSE)+
      scale_shape_manual(values=c("FALSE"=4,"TRUE"=16),
                         labels=c("Observations","Censored"),
                         limits=c("FALSE","TRUE"), drop = FALSE)+
      theme(legend.position.inside = c(0.2, 0.8),
            legend.background = element_rect(fill = "white",colour="black"),
            legend.title = element_blank(),
            plot.title = element_text(hjust=0.5))+
      ggtitle(mytitle_1)
    if(!is.na(DateRange[1])){
      myplot_1<-myplot_1+labs(subtitle = paste0(DateRange[1]," to ",DateRange[2]," (",
                                                round(as.numeric(DateRange[2]-DateRange[1])/365,0)," years)"))+
        theme(plot.subtitle = element_text(hjust=0.5))}
    
    plotlist[[1]]<-myplot_1
    #2. Make plot of matrix of values ##########
    # if(myTimeIncr!="Annual"){ #can't make matrix plots for annual data
      mytitle_2<-paste0("Matrix of values")
    if(!is.na(mymain)){mytitle_2<-paste0(mymain,": ",mytitle_2)}
    PlotData$YearFac<-factor(PlotData$Year,levels=rev(seq(min(PlotData$Year),max(PlotData$Year))))
    
    myplot_2<-ggplot(data=PlotData, aes(x=TimeIncr,y=YearFac,fill=V1))+
      geom_tile(colour="black",linewidth=0.1,show.legend = T)+theme_bw()+ylab(Year)+xlab("Month")+
      scale_fill_gradient(low="lemonchiffon",high="red2",na.value="white")+
      labs(fill="Value")+scale_x_discrete(expand = c(0, 0)) +
     scale_y_discrete(expand = c(0, 0)) + 
      theme(axis.ticks = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 0),
            plot.title = element_text(hjust=0.5),
            legend.position="bottom")+
      ggtitle(mytitle_2)
    
    plotlist[[2]]<-myplot_2
    
    
    #3. Make plot of matrix of cesnored values ##########
      mytitle_3<-paste0("Matrix of censoring")

    if(!is.na(mymain)){mytitle_3<-paste0(mymain,": ",mytitle_3)}
    
    myplot_3<-ggplot(data=PlotData, aes(x=TimeIncr,y=YearFac,fill=Censored))+
      geom_tile(colour="black",linewidth=0.1,show.legend = T)+theme_bw()+ylab(Year)+xlab("Month")+
      scale_fill_manual(values=c("FALSE"="skyblue2","TRUE"="red"),
                          labels=c(" NOT Censored"," Censored"),
                          limits=c("FALSE","TRUE"),na.value = "white", drop = FALSE)+
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) + 
      theme(axis.ticks = element_blank(), legend.title = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            plot.title = element_text(hjust=0.5),legend.position="bottom")+
      ggtitle(mytitle_3)
    plotlist[[3]]<-myplot_3
    
    #4. Make matrix of number of samples per month
    
    mytitle_4<-paste0("Matrix of obs. per month")
    if(!is.na(mymain)){mytitle_4<-paste0(mymain,": ",mytitle_4)}
    
    myplot_4<-ggplot(data=PlotData, aes(x=TimeIncr,y=YearFac,fill=nVal))+
      geom_tile(colour="black",linewidth=0.1,show.legend = T)+theme_bw()+ylab(Year)+xlab("Month")+
      scale_fill_gradient(low="pink",high="maroon4",na.value="white")+
      labs(fill="Count")+scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) + 
      theme(axis.ticks = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 0),
            plot.title = element_text(hjust=0.5),legend.position = "bottom")+
      ggtitle(mytitle_4)
    
    plotlist[[4]]<-myplot_4
  }else{
    print("*** WARNING: ggplot list not produced as ggplot2 not loaded ***")
  }
  return(plotlist)
}

#######################################################################
# Trend analysis utility functions
#######################################################################
#These pieces of code were previously imbedded within the Kendall and Sen
# functions, but were repeated.  They have been extracted in order to ensure 
# that there is consistency in their application if/when updates are required


ApplyHiCensor <-function(x, ValuesToUse = "RawValue",HiCensor=TRUE,RawValues=TRUE){  
# HiCensor can either be logical, or numeric.  When ==TRUE all values below the highest left censored value are set as censored.
  #This code is only called if HiCensor is numeric OR HiCensor
# when numeric, only values below the value are reassigne as censored if they are not already so identified.
  
  CenMax <- max(x[, ValuesToUse][x$CenType == "lt"], na.rm=T)
  
  if (!is.logical(HiCensor)){ # if there is a numeric high censor value, check to see whether the largest censored value reaches this - if not, hi censor level to largest censored value
  CenMax=min(CenMax,HiCensor)}
  
  x$Censored[x[, ValuesToUse] < CenMax] <- TRUE
  x$CenType[x[, ValuesToUse] < CenMax] <- "lt"  # and treat these as all censored. 
  if(RawValues==TRUE){  #If the data are Flow adjusted, we just want to get the new censoring labels and not change the values
  x[, ValuesToUse][x[, ValuesToUse] < CenMax] <- CenMax # set all values below the max censored value, to the max censored
  }
  return(x)
}

GetTimeIncrYear<-function(x,ValuesToUse = "RawValue",Year="Year"){
  
  x <- x[order(x$myDate), ]  # make sure the data are in date order!
  if(length(unique(x$TimeIncr))==1){ 
    x$TimeIncrYear<-x[[Year]]
    x$myFreq<-1
  }else{
    x$TimeIncrYear <- paste(x[["TimeIncr"]], x[[Year]], sep = "-")
    x$myFreq<-length(unique(x$TimeIncr))
  }
  take <- !is.na(x[, ValuesToUse]) & !is.na(x$Censored)  # remove values and censored values that are NA.
  x <- x[take, ] # drop the NA values
  return(x)
}

ValueForTimeIncr <-function(x, ValuesToUse = "RawValue",Year="Year",UseMidObs=TRUE,TimeIncrMed=TRUE){
  if(is.null(x$Season)){x$Season<-x$TimeIncr} #Just a catch in case season has not been assigned yet
  # this takes the median of observations and DATES in time increment (so a time increment only has one value)
  # if the median of the raw values is less than OR EQUAL to the max censored value, call the value censored, otherwise NOT censored.
  
  #If the "UseMidObs" variable is set to "TRUE" - the observation that is closest in time to the middle of the time increment is used (rather
  # than takin the median of all observations over time increments).  This variant is appropriate to use when there has a been a systematic change
  # in sampling interval over the trend period
 if(TimeIncrMed){
   if (UseMidObs){
    #Determine mid-dates for time increments.  NOTE, this will only work if seasons are based on Months and all time increments are the same
    #nubmer of months long
    nTimeIncr=length(levels(x$TimeIncr));nTimeIncr=max(nTimeIncr,1)
    midTimeIncr=data.frame(theMidDate=MidTimeIncrList[[nTimeIncr]])
    midTimeIncr$TimeIncr<-levels(x$TimeIncr)
    
    x<-merge(x,midTimeIncr)
    
    Data <-  ddply(x, "TimeIncrYear", function(y) {
      if(length(y$TimeIncrYear)>1){
      midTimeIncrdate<-as.Date(paste0(y$Year[1],"-01-01"))+y$theMidDate[1]
      #Find date closest to the middle of the time increment
      theid<- which(abs(y$myDate-midTimeIncrdate)==min(abs(y$myDate-midTimeIncrdate)))[1]   
      }else{
      theid<-1
      }

    return(data.frame(V1 = y[theid,ValuesToUse], NewDate = y$myDate[theid], Censored = y$Censored[theid], Year=y[theid,Year], 
                      Month=y[theid,"Month"], BiMonth=y[theid,"BiMonth"],Qtr=y[theid,"Qtr"],BiAnn=y[theid,"BiAnn"],
                      TimeIncr=y$TimeIncr[theid],Season=y$Season[theid],CenType=as.character(y$CenType[theid])))
    })
    
  }else{
  Data <-  ddply(x, "TimeIncrYear", function(y) {
    Median <- median(y[, ValuesToUse]) # get the median of the values
    NewDate <- median(y[, "myDate"]) # get the median of the dates
    if(sum(y$Censored) == 0) {
      Censored = FALSE  
    } else { 
      MaxCensored <- max(y[, ValuesToUse][y$Censored],na.rm = T)
      Censored <- ifelse(Median <= MaxCensored, TRUE, FALSE) 
    }
    return(data.frame(V1 = Median, NewDate = NewDate, Censored = Censored, Year=y[1,Year], 
                      Month=y[1,"Month"], BiMonth=y[1,"BiMonth"],Qtr=y[1,"Qtr"],BiAnn=y[1,"BiAnn"],
                      TimeIncr=y$TimeIncr[1],Season=y$Season[1]))
  })
  # this determines censor type for Data  based on CenType that has highest prevalence in the time increment
  Data$CenType <- ddply(x, "TimeIncrYear", function(y) names(which.max(table(y$CenType))))[,2]# 
  }
 
 }else{
   Data=(data.frame(TimeIncrYear=as.character(x[,"TimeIncrYear"]),V1 = x[,ValuesToUse], NewDate = x$myDate, 
                    Censored = as.character(x$Censored), 
                     Year=x[,Year], TimeIncr=as.character(x$TimeIncr),Season=as.character(x$Season),CenType=as.character(x$CenType)))
 }
  return(Data)
  }

GetAnalysisNote<-function(Data, ValuesToUse = "RawValue",IsSeasonal=FALSE,SecondTierTest=FALSE,...){
  #This function is used to check that there are sufficient data at several points in the analysis
  #to contintue with the calculations, and if there are not to provide an error message describing the issue
  AnalysisNote="ok"
  #Set generic levels here about numbers of unique values and non-cenosred values required
  #allows user to adjust and apply consistently across all functions
  noUnique<-3
  noNonCen<-5
  noUniqueSeas<-2
   #First round of filtering
  if(all(is.na(Data[, ValuesToUse]))){
    AnalysisNote="Data all NA values"
    
    }else if(length(unique((Data[, ValuesToUse][!Data$Censored&!is.na(Data[,ValuesToUse])]))) < noUnique){
    AnalysisNote=paste0("< ",noUnique," unique values")
    
    }else if(length(which(Data[,'Censored']==FALSE & !is.na(Data[,ValuesToUse]))) < noNonCen){
      AnalysisNote=paste0("< ",noNonCen," Non-censored values")
  }
     
  #Check to see whether we have failed at the first step, AND that we are requested to do next
  #set of tests (these are for data post averagein over seasons.  Continue on to check more details 
  if(AnalysisNote=="ok"&SecondTierTest){
    #Different tests whether Seasonal or not
    if(IsSeasonal){ #Seasonal checks
      #Check to see that there are sufficient non-NA noncensored values in each season
      EnoughSeason <- min(table(Data[,"Season"][!is.na(Data[,ValuesToUse])]))< noUnique
      #Check to see that there are sufficient unique values in each season
      EnoughSeason_2 <- min(ddply(Data,c("Season"),function(x)length(unique(x[, ValuesToUse])))[,"V1"]) <  noUniqueSeas #There must be at least 3 unique values per season

      if (EnoughSeason==TRUE){
        AnalysisNote=paste0("< ",noUnique,"non-NA values in Season")
      }else if (EnoughSeason_2==TRUE){
        AnalysisNote=paste0("< ",noUnique," unique values in Season")
      }else{
        #Then check the run length
        # if TRUE data has long sequence of indentical values within one or more season and will crash 
        RunLengthSeason <- ddply(Data, "Season", function(y) {  # for each season y = x[x$Season == "Q1",]
          take <- !is.na(y[, ValuesToUse])  & !is.na(y$Censored)
          theValues <- y[take, ValuesToUse]
          if(length(theValues)>1) { 
            longRun <- max(unlist(rle(diff(theValues))[1]))/length(theValues) > 0.75  # this catches seasons with one value to avoid warning. 
          } else {
            longRun <- FALSE
          }
          return(longRun)
        })
        RunLength <- sum(RunLengthSeason$V1)>0 
        
        if(RunLength==TRUE) {AnalysisNote="Long run of single value in a Season"}
      }
    }else{#Then not seasonal 
      # if TRUE data has long sequence of indentical values and will crash 
      RunLength <-  max(unlist(rle(diff(Data[, ValuesToUse]))[1]))/length(Data[, ValuesToUse]) > 0.5 
      if(RunLength==TRUE) {AnalysisNote="Long run of single value"}
    }
  }
  
  return(AnalysisNote)
}

GetTrendDirectionandClass<-function(A3){
  TrendClass <- "Insufficient Data" # default assumption is  that there is insufficient data to define trend direction
  if(A3$Cd >= 0.95)  TrendClass <- "Decreasing"
  if(A3$Cd <= 0.05)  TrendClass <- "Increasing"
  if(is.na(A3$C))  TrendClass <- "Not Analysed"
  
  TrendDirection<- "Indeterminate" # rare - but dooes occur!
  if(A3$Cd > 0.5 )  TrendDirection <- "Decreasing"
  if(A3$Cd < 0.5)  TrendDirection <- "Increasing"
  if(is.na(A3$S))  TrendDirection <- "Not Analysed"
  return(data.frame(TrendDirection=TrendDirection,TrendCategory=TrendClass)) 
}

GetInterObservationSlopes<-function(Data,RawValues=TRUE){
  # number of Time Increments between observations, this excludes comparing observations in the same year+tiem increment
  if(is.null(Data$Season[1])){x$Season<-x$TimeIncr} #Catch in case season is not specified
  
  Data$SeasYr<-apply(Data,1,function(y) paste0(y['Season'],y['Year']))
  AllTimeIncs <- outer(as.numeric(Data$NewDate), as.numeric(Data$NewDate), `-`)/365.25 # the time increment in years
  SameYearSeason <- outer(Data$SeasYr, Data$SeasYr, "==")
  AllTimeIncs[SameYearSeason] <- NA # remove any increments withing the same year + season
 
  # take each observation and compute differences for all other observations
  #Using the RAW data - this will be used for the checks later to account for changing censor levels
  AllDifferences1 <- outer(Data$V1, Data$V1, `-`)
  AllDifferences1[SameYearSeason] <- NA # remove any distances withing the same year + season

  CenLab <- outer(Data$CenType, Data$CenType, 'paste')
  Slopes1 <- AllDifferences1/AllTimeIncs
  
  # For the purposes of of computing the sen slope, replace "lt" censored values with half 
  # the censored values and "gt" values with 1.1 the censored value (to aoiv ties with non-cenosred obs)
  if(RawValues==TRUE){ # We don't do this step if the data is FA - this is already done as part of the flow adjustement
    Data[Data$CenType=="lt","V1"]<-Data[Data$CenType=="lt","V1"]*0.5
    Data[Data$CenType=="gt","V1"]<-Data[Data$CenType=="gt","V1"]*1.1
    
    AllDifferences <- outer(Data$V1, Data$V1, `-`)
    AllDifferences[SameYearSeason] <- NA # remove any distances withing the same year + season
    Slopes <- AllDifferences/AllTimeIncs
    
  }else{
    AllDifferences<-AllDifferences1
    Slopes<-Slope1
  }

  OUTPUTS <- data.frame(Slopes=as.vector(Slopes [lower.tri(Slopes, diag = FALSE)]),
                      CensorLabel=as.vector(CenLab [lower.tri(CenLab, diag = FALSE)]),
                      Slopes1=as.vector(Slopes1 [lower.tri(Slopes1, diag = FALSE)]))
  
  #NOW - tidy this up to account for censoring
  #1.  There can be NO slope between any pairs of censored values:
  OUTPUTS$Slopes[OUTPUTS$CensorLabel %in% c("gt gt","lt lt")]<-0
  #2. There cannot be a slope between a non-censored value and a censored higher value (i.e., a positive slope)
  OUTPUTS$Slopes[OUTPUTS$Slopes1>0 & OUTPUTS$CensorLabel %in% c("lt not")]<-0
  #3. There cannot be a slope between a censored value and a non-censored lower value (i.e., a negative slope)
  OUTPUTS$Slopes[OUTPUTS$Slopes1<0 & OUTPUTS$CensorLabel %in% c("not lt")]<-0
  #4. There cannot be a slope between a right censored value and a non-censored higher value (i.e., a positive slope)
  OUTPUTS$Slopes[OUTPUTS$Slopes1>0 & OUTPUTS$CensorLabel %in% c("not gt")]<-0
  #5. There cannot be a slope between a non-censored value and a right censored lower value (i.e., a negative slope)
  OUTPUTS$Slopes[OUTPUTS$Slopes1<0 & OUTPUTS$CensorLabel %in% c("gt not")]<-0
  
  OUTPUTS <- OUTPUTS[!is.na(OUTPUTS$Slopes),]
  OUTPUTS$Slopes1<-NULL
  return(OUTPUTS)
}

GGPlotTrend<-function(x,x1,mymain=NULL,Intercept=NA,ValuesToUse = "RawValue",AnnualSenSlope=NA,
                    Lci=NA,Uci=NA,IsSeasonal=FALSE,myTimeIncr=NA,mySeason=NA,Ymed=NA,Tmed=NA,Percent.annual.change=NA,Cd=NA,legend.pos="top",
                    the_ylab="Values",CropYaxis=FALSE,UseMidObs=TRUE){
  if(is.null(x$myDate)) x$myDate<-x$NewDate
  
  Ymed<-median(x[,ValuesToUse])
  Interceptu <- Ymed - (Uci*Tmed)
  Interceptl <-Ymed - (Lci*Tmed)
  
  # value at end of time period
  T1 <- Intercept + as.numeric(diff(range(x$myDate)))/365.25 * AnnualSenSlope
  T1l <- Interceptl + as.numeric(diff(range(x$myDate)))/365.25 * Lci
  T1u <- Interceptu + as.numeric(diff(range(x$myDate)))/365.25 * Uci
  
  Trends<-data.frame(t=as.Date(c(rep.int(as.character(min(x$myDate)),3),rep.int(as.character(max(x$myDate)),3)),origin="1-1-1970"),
                     y=c(Intercept,Interceptl,Interceptu,T1,T1l,T1u),
                     type=c("Trend","90% C.I.","90% C.I.","Trend","90% C.I.","90% C.I."),
                     gp=c("a","b","c","a","b","c"))
  if(UseMidObs==TRUE){
  x$DataType<-"Observations\n(mid of time increment)"}else{
    x$DataType<-"Observations\n(median of time increment)"
  }
  names(x)[names(x)=="NewDate"]<-"myDate"
  x1$DataType<-"Raw Observations"
  x1$CensoredOrig<-x1$Censored
  
  xall<-rbind(x[,c("myDate",ValuesToUse,"DataType","CensoredOrig")],
              x1[,c("myDate",ValuesToUse,"DataType","CensoredOrig")])
  names(xall)[c(2,4)]<-c("Values","Censored")
  xall$Censoring<-"Non-censored"
  xall$Censoring[xall$Censored=="TRUE"]<-"Censored"
  xall$Censoring<-factor(xall$Censoring,levels=c("Censored","Non-censored"))
  
  cencol<-c("red","black")
  names(cencol)<-c("Censored","Non-censored")
  
  myplot<-ggplot()+#xall,aes(x=myDate,y=Values))+
    geom_point(data=xall,aes(x=myDate,y=Values,colour=Censoring,shape=DataType,alpha=DataType),show.legend = T)+
    scale_color_manual(values=cencol, drop = FALSE)+
    scale_shape_manual(values=c(16,21), drop = FALSE)+
    scale_alpha_manual(values=c(1,0.5), drop = FALSE)+
    theme_bw()+ylab(the_ylab)+
    geom_line(data=Trends,aes(x=t,y=y,group=gp,linetype=type),colour="blue",linewidth=.9)+
    scale_linetype_manual(values=c(2,1))+
    labs(linetype="Trends",colour="Censoring",alpha=NULL,shape="Data Type")+
    xlab("Time")+ylab("Values")+guides(alpha="none")+
    theme(legend.position = "bottom",legend.box="horizontal",
          legend.direction = "vertical",
          legend.spacing=unit(0.1,"lines"),
          legend.box.just = "top",
          legend.title=element_text(hjust=0.5))+
    guides(colour = guide_legend(override.aes = list(shape = 15,size=3)),
           linetype = guide_legend(override.aes = list(shape = NA)))

  SenLabel<-ifelse(IsSeasonal,"Annual Seasonal Sen Slope","Annual Sen Slope")
  mysub<-ifelse(IsSeasonal,paste0(myTimeIncr," time increments; ",mySeason," seasons"),paste0(myTimeIncr," time increments; non-seasonal"))
  if(is.null(mymain)==FALSE){mymain<-paste0(mymain,"\n",mysub)}else{
    mymain<-mysub}
  myplot<-myplot+ggtitle(mymain)+theme(plot.title = element_text(hjust = 0.5))
  
  
mylegend=c(paste("% ",SenLabel," = ", round(Percent.annual.change,1),"%\n",
                 SenLabel," = ", signif(AnnualSenSlope,3),
                 "\nConfidence trend is decreasing = ",round(Cd,3))) 
xp=diff(range(xall$myDate))*0.03+min(xall$myDate)

if(CropYaxis !=FALSE){
gb=ggplot_build(myplot)
ymin = gb$layout$panel_params[[1]]$y.range[1]
ymax1=gb$layout$panel_params[[1]]$y.range[2]
ymax=max(min(ymax1,CropYaxis*sqrt(var(xall$Values))+median(xall$Values)),max(Trends$y))

myplot<-myplot+ylim(c(ymin,ymax))
if(ymax<ymax1){
mynote=c(paste0("Note: y-axis cropped at ",CropYaxis," s.d.\n",
                 "Max obs. value  = ", signif(max(xall$Values),3))) 

myplot<-myplot + 
  geom_label(data=data.frame(label=mynote,x=as.Date(Inf), y=as.Date(Inf)),
             aes(label=label,x=x, y=y),fill="white",label.size=0,vjust="inward",
             hjust="inward",label.padding=unit(0.7,"lines"))
}
}

myplot<-myplot + 
  geom_label(data=data.frame(x=as.Date(-Inf), y=as.Date(Inf),label=mylegend),aes(
             label=label,x=x, y=y),fill="grey95",label.size=0,vjust="inward",
             hjust=0,label.padding=unit(0.7,"lines"))

# 
# myplot<-myplot + 
#   geom_label(label=mylegend,x=-Inf, y=Inf,fill="grey95",label.size=0,vjust="inward",
             # hjust=0,label.padding=unit(0.7,"lines"))
return(myplot)

}


#######################################################################
#Wrapper functions
#######################################################################

# Wrapper function to incorporate seasonality assessment and then selecting appropriate method
doMyTrends_v2502<-function(x,TimeIncrOpts=NA,do.plot=F,
                           propIncrTol=0.9,propYearTol=0.9,AnnualOK=c("MCI"),
                           TrendPeriod=NA,EndYear=NA,Year="Year",
                           ...){ # elipsis { ... } passes arguments to the internal functions. (e.g., Year="Year",HiCensor=F,ValuesToUse = "RawValues", do.plot=F)
  
  #Process data to cut down to the trend period, and add an appropriate time increment
  x<-InspectTrendData(x,Year=Year,TrendPeriod = TrendPeriod,EndYear=EndYear,propIncrTol=propIncrTol,
                      propYearTol=propYearTol,AnnualOK = AnnualOK,TimeIncrOpts=TimeIncrOpts)
  
  #Does the site meet minimum data requirements?
  if(is.null(dim(x))){ #Then there are no observations in the trend period
    return(NULL)
  }else if(x$TimeIncr[1]=="none"){ #then there are not enough observations in the trend period
    return(NULL)
  }else{
  
  x<-GetSeason(x,do.plot=FALSE,...)[[1]]
  # print(MySeasonality)
  if(x$Seasonal[1]==FALSE){ #Then we will perform non-seasonal trend tests
    OUTPUTS<-NonSeasonalTrendAnalysis(x,do.plot=do.plot,...)
  }else{ #Otherwise it's seasonal
    OUTPUTS<-SeasonalTrendAnalysis(x,do.plot=do.plot,,...)
  }
  if(is.data.frame(OUTPUTS)&do.plot){ #This is when you have asked for plots, but the trend assessment is not analysed
    OUTPUTS<-list(OUTPUTS,NA)
  }
  return(OUTPUTS)
}}


NonSeasonalTrendAnalysis<-function(x,do.plot=F,...){
  
  A1<-MannKendall(x,...)
  if(A1$AnalysisNote =="ok"){#Then carry on and do the Sen Slope test
    A1$AnalysisNote <- NULL
    A2<-SenSlope(x,Cd=A1$Cd,do.plot=do.plot,...)
  }else{
    A2<-data.frame(Median=NA, Sen_VarS=NA, AnnualSenSlope=NA, Intercept=NA, 
                   Sen_Lci=NA, Sen_Uci=NA,Sen_Probability=NA, Sen_Probabilitymax=NA, Sen_Probabilitymin=NA,
                   Percent.annual.change=NA)
  }
  if(class(A2)=="list"){
    A<-cbind.data.frame(A1,A2[[1]])
  }else{
  A<-cbind.data.frame(A1,A2)}
  
  #Add on trend category and trend direction
  if(is.na(A$C)){
    A$TrendCategory<-"Not Analysed"
    A$TrendDirection<-"Not Analysed"
  }else{
  A<-cbind.data.frame(A,GetTrendDirectionandClass(A))}
 

  if(class(A2)=="list"){
    A<-list(A,A2[[2]])
  }
  
  #Choose to no longer return some of the older outputs.  They are retained in case anyone wants to find them
  #for comparative purposes (back compatability to older versions), but are no longer part of the core outputs
 
  if(class(A2)=="list"){
    A[[1]][,c("Sen_Probability","Sen_VarS","Intercept","Sen_Probabilitymax","Sen_Probabilitymin","TrendCategory")]<-NULL
  }else{

   A[,c("Sen_Probability","Sen_VarS","Intercept","Sen_Probabilitymax","Sen_Probabilitymin","TrendCategory")]<-NULL}
  
  return(A)
  
}

SeasonalTrendAnalysis<-function(x,do.plot=F,...){
  #Do Seasonal Kendall Test
  A1<-SeasonalKendall(x,...)
  
  if(A1$AnalysisNote =="ok"){
    #Continue on to do Seasonal Sen Slope if Seasonal Kendall Test was sucessful
    A1$AnalysisNote <- NULL
    A2<-SeasonalSenSlope(x,Cd=A1$Cd,do.plot=do.plot,...)
    
  }else{
    A2<-data.frame(Median=NA, Sen_VarS=NA, AnnualSenSlope=NA, Intercept=NA, 
                   Sen_Lci=NA, Sen_Uci=NA,Sen_Probability=NA, Sen_Probabilitymax=NA, Sen_Probabilitymin=NA,
                   Percent.annual.change=NA)
    # if(do.plot=="doGGPlot"){A2<-list(A2,NA)}
  }
  
  if(class(A2)=="list"){
    A<-cbind.data.frame(A1,A2[[1]])
  }else{
    A<-cbind.data.frame(A1,A2)}
  
  if(is.na(A$Cd)){
    A$TrendCategory<-"Not Analysed"
    A$TrendDirection<-"Not Analysed"
  }else{
    A<-cbind.data.frame(A,GetTrendDirectionandClass(A))}
  
  if(class(A2)=="list"){
    A<-list(A,A2[[2]])
  }
  
  #Choose to no longer return some of the older outputs.  They are retained in case anyone wants to find them
  #for comparative purposes (back compatability to older versions), but are no longer part of the core outputs
  if(class(A2)=="list"){
    A[[1]][,c("Sen_Probability","Sen_VarS","Intercept","Sen_Probabilitymax","Sen_Probabilitymin","TrendCategory")]<-NULL
  }else{
    
    A[,c("Sen_Probability","Sen_VarS","Intercept","Sen_Probabilitymax","Sen_Probabilitymin","TrendCategory")]<-NULL}
  return(A)
}

######################################################################
# Determine what season to use
######################################################################
#Make a list of the seasons that can be used for each of the time increments
TimeIncrToSeas<-list()
TimeIncrToSeas[["Monthly"]]<-c("Monthly","Bi-monthly","Quarterly","Bi-annual")
TimeIncrToSeas[["Bi-monthly"]]<-c("Bi-monthly","Bi-annual")
TimeIncrToSeas[["Quarterly"]]<-c("Quarterly","Bi-annual")
TimeIncrToSeas[["Bi-annual"]]<-c("Bi-annual")

#DefaultSeason labels
SeasonLabs<-data.frame(mySeas=c("Monthly","Bi-monthly","Quarterly","Bi-annual"),
                       SeasName=c("Month","BiMonth","Qtr","BiAnn"))

#This function determines which time increment to use as a season.  The season can 
# be as small as the analysis time increment (TimeIncr), or any other time increment
# that is a multiple of TimeIncr.  It runs the seasonality test 
GetSeason<-function(x,ValuesToUse = "RawValue",RawValues=TRUE,printKW=FALSE,do.plot=FALSE,
                    UseMidObs=TRUE,TimeIncrMed=TRUE,...){
  
x<-x[!is.na(x[,ValuesToUse]),]
  
  minSeas<-WhatistheIncr(x,val="TimeIncr") #Just get the name of the minimum time increment

if(minSeas!="Annual"){
  StartSeas<-which(SeasonLabs$mySeas==minSeas)
  mylist<-list()
for(i in TimeIncrToSeas[[minSeas]]){
  x$Season<-x[,SeasonLabs$SeasName[SeasonLabs$mySeas==i]]
  mylist[[SeasonLabs$SeasName[SeasonLabs$mySeas==i]]]<-SeasonalityTest(x,do.plot=FALSE,ValuesToUse = ValuesToUse,RawValues = RawValues,
                                                                       UseMidObs=UseMidObs,TimeIncrMed=TimeIncrMed,...)
}
  
  MyTab<-ldply(mylist)
  if(printKW){
    print(MyTab)
  }

  #Look for the smallest possible time increment that is seasonal and set this to season.
  # if none are seasonal, set to the TimeIncr
  if(any(MyTab$p<0.05&MyTab$SeasNote=="ok")){
    MyTab<-MyTab[MyTab$SeasNote=="ok"&MyTab$p<0.05,]
    MyTab<-MyTab[order(MyTab$KWstat,decreasing = TRUE),]
    x$Season<-x[,MyTab$.id[1]]
    
    temp1<-SeasonalityTest(x,do.plot=do.plot,ValuesToUse = ValuesToUse,RawValues = RawValues,
                           UseMidObs=UseMidObs,TimeIncrMed=TimeIncrMed,...)
    x$Seasonal=TRUE
  }else{
    x$Season<-x$TimeIncr
    x$Seasonal=FALSE
    temp1<-SeasonalityTest(x,do.plot=do.plot,ValuesToUse = ValuesToUse,RawValues = RawValues,
                           UseMidObs=UseMidObs,TimeIncrMed=TimeIncrMed,...)
  }

}else{ #Then we have an annual time step and we can't do seasonality test
  x$Season<-x$TimeIncr
  x$Seasonal<-FALSE
  temp1<-"Annual time increment"
}
  out<-list()
  out[[1]]<-x
  if(class(temp1)=="list"){
  out[[2]]<-temp1[[1]]
  out[[3]]<-temp1[[2]]
  }else{out[[2]]<-temp1}
  return(out)
}


#######################################################################
#Seasonality Test
#######################################################################
SeasonalityTest <-  function(x, ValuesToUse = "RawValue", HiCensor=FALSE,main=NULL,Year="Year",
                             do.plot=FALSE,mymain=NA,RawValues=TRUE,UseMidObs=TRUE,TimeIncrMed=TRUE,...) { # x = DataForTest
  

  #Stop if season is not defined
  if(is.null(x$Season[1])){stop("Season must be defined", call. = FALSE)}
  # 
  #Stop if season is not a multiple of time increment
  if(!is.null(x$Season[1])){
  check<-ddply(x,.(TimeIncr),function(x) data.frame(myCount=length(unique(x$Season))))
  if(max(check$myCount>1)){stop("Season must be a multiple of the time increment", call. = FALSE)}}
  
  Observations <- sum(!is.na(x[, ValuesToUse]))
  mySeas<-WhatistheIncr(x,"Season")
  myTimeIncr<-WhatistheIncr(x,"TimeIncr")
  x<-GetTimeIncrYear(x,ValuesToUse,Year)
  
 x<-ValueForTimeIncr(x,ValuesToUse=ValuesToUse,Year=Year,UseMidObs=UseMidObs,TimeIncrMed = TimeIncrMed)
  # if(HiCensor==TRUE|is.numeric(HiCensor)) x<-ApplyHiCensor(x,"V1",HiCensor,RawValues)
  # CHECK TO SEE IF DATA OK TO CONTINUE
  AnalysisNote <- GetAnalysisNote(x,ValuesToUse="V1",IsSeasonal = TRUE,SecondTierTest = TRUE)
  
  if(Observations > 5 & length(unique(x[, "V1"]))>1&AnalysisNote=="ok") { # check there is minimum data and it has some variation
    
    #If required, apply the hi-censor

      
    if(do.plot==TRUE){
      x$Value<-x[, "V1"]
      thetitle<-paste0(ifelse(is.na(mymain),"Seasonality plot",paste0(mymain,": Seasonality plot")),"\n",
                       myTimeIncr," time increments, ",mySeas," seasons")
      
        myplot<-ggplot(x,aes(y=Value,x=Season))+geom_boxplot(fill="lightblue")+
          theme_bw()+ggtitle(thetitle)+theme(plot.title = element_text(hjust = 0.5))
      }
    
    if(length(table(x[, "Season"][!is.na(x[, "V1"])]))>1)
    { # check there is more than one season represented 
      #if(!is.factor(x$Season))  x$Season<-factor(as.character(x[, "Season"]),levels=unique(x[,"Season"]))
      
      KW <- kruskal.test(x[, "V1"] ~ factor((x[, "Season"])))
      OUT <- data.frame(Observations=Observations, KWstat = KW$statistic, 
                        pvalue = KW$p.value,SeasNote=AnalysisNote,
                        TimeIncr=myTimeIncr,Season=mySeas)
      
      if(do.plot==TRUE){
        mybox<-paste("Observations =", Observations,
                       "\nKruskal Wallis statistic = ", round(OUT$KWstat,3),
                     "\nP value =", round(OUT$pvalue,3))
        myplot<-myplot + 
        geom_label(label=mybox,x=-Inf, y=Inf,fill="grey95",label.size=0,vjust="inward",
                   hjust=0,label.padding=unit(0.7,"lines"),label.r=unit(0,"lines"))

        OUT<-list(OUT,myplot)
        }
    } else {
      OUT <-data.frame(Observations=Observations, KWstat = NA, pvalue = NA,SeasNote=AnalysisNote,TimeIncr=myTimeIncr,Season=mySeas)
    }
    
  } else {
    OUT <- data.frame(Observations=Observations, KWstat = NA, pvalue = NA,SeasNote=AnalysisNote,TimeIncr=myTimeIncr,Season=mySeas)
  }
  OUT$TimeIncr<-myTimeIncr
  OUT$Season<-mySeas
  return(OUT)
}

#######################################################################
# Mann-Kendall test.
#######################################################################
MannKendall <- function(x, ValuesToUse = "RawValue", Year = "Year", 
                        RawValues=TRUE,UseMidObs=TRUE,TimeIncrMed=TRUE,...) { # 
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have Observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the time increment (name =  "TimeIncr") and the year (name = "Year" )
  
  # CHECK TO SEE IF DATA IS OK TO CONTINUE
  AnalysisNote <- GetAnalysisNote(x,ValuesToUse)
  
  if(AnalysisNote !="ok") { # if any TRUE dont do the test 
    KendalTest <- data.frame(nObs = nrow(x), VarS=NA, S = NA, D = NA, tau = NA, Z=NA, p=NA,C=NA,Cd=NA)
  } else {
  
  #..............................................................................#
  #Tidy data and add on TimeInceYear column
  x<-GetTimeIncrYear(x,ValuesToUse,Year)
  x$Season<-x$TimeIncr #This is just a catch so the a different season is not sent through the analysis
  x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])
  
  #Take medians over each time increment if required
  Data<-ValueForTimeIncr(x,ValuesToUse,Year,UseMidObs,TimeIncrMed)

  #..............................................................................#
  
  # CHECK TO SEE IF DATA OK TO CONTINUE
    AnalysisNote <- GetAnalysisNote(Data,ValuesToUse="V1",SecondTierTest = TRUE)
    
  if(AnalysisNote != "ok" ) { # if any TRUE dont do the test 
    KendalTest <- data.frame(nObs = nrow(x),nTimeIncr = nrow(Data), VarS=NA, S = NA, D = NA,tau = NA, Z=NA, p=NA,C=NA,Cd=NA)
  } else {
    
    # organise the ordering of data
    Data <- Data[order(Data$NewDate), ] # restores the time ordering

    KendalTest <- GetKendal(x = Data)
    KendalTest$C<-1-KendalTest$p/2
    KendalTest$Cd<-KendalTest$C
    KendalTest$Cd[which(KendalTest$S>0)]<-KendalTest$p[which(KendalTest$S>0)]/2
    names(KendalTest)[names(KendalTest)=="vars"]<-"VarS"
    KendalTest<-cbind(data.frame(nObs=nrow(x)),as.data.frame(KendalTest))
    if(is.na(KendalTest$C)==T) AnalysisNote = "Not Analysed"
  }
}  # end of first if statement

  KendalTest$AnalysisNote<-AnalysisNote
  KendalTest$prop.censored <- sum(x$Censored)/nrow(x) # proportion of values that were censored
  KendalTest$prop.unique <- length(unique(x[, ValuesToUse]))/nrow(x) # proportion of values that were censored
  KendalTest$no.censorlevels<-length(unique(x[x$Censored=="TRUE",ValuesToUse]))
  KendalTest$TimeIncr<-WhatistheIncr(x,"TimeIncr")
  KendalTest$SeasIncr<-"NonSeasonal"
  # KendalTest$HiCensorNote<-HiCensor
  
  return(KendalTest)
}

#######################################################################
# Sen slope test
#######################################################################
SenSlope <- function(x, ValuesToUse = "RawValue", ValuesToUseforMedian=NA,Year = "Year", 
                     do.plot=FALSE,mymain=NULL,RawValues=TRUE,UseMidObs=TRUE,TimeIncrMed=TRUE, ...) {  # x=DataForTest
  # Calculates the annual Sen slope and its 90% CIs - for 95% confidence in trend direction. 
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the time increment (name =  "TimeIncr") and the year (name = "Year" )
  
  # CHECK TO SEE IF DATA OK TO CONTINUE
  AnalysisNote <- GetAnalysisNote(x, ValuesToUse = ValuesToUse,...)
  
  if(AnalysisNote != "ok") { # if any TRUE dont do the test and return NA values
    Median<-NA; VarS<-NA; AnnualSenSlope<-NA; Intercept<-NA; Lci<-NA; Uci<-NA;  
    Sen_Probability<-NA; Sen_Probabilitymin<-NA;Sen_Probabilitymax<-NA; Percent.annual.change<-NA
  } else {  
  #..............................................................................#
    #Tidy data and add on TimeIncrYear column
    x<-GetTimeIncrYear(x,ValuesToUse,Year)
    myFreq<-x$myFreq[1]; x$myFreq<-NULL
    x1<-x;x1$V1<-x1[,ValuesToUse]# Save raw data here for plotting later

     x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])
    


    #Take medians or mid-values over each time increment
    Data<-ValueForTimeIncr(x,ValuesToUse,Year,UseMidObs,TimeIncrMed)

    if(is.na(ValuesToUseforMedian)){
      Median <- median(Data$V1, na.rm=T) #The median of the data after summarising to time increments
     }else{
      Median <- median(x[,  ValuesToUseforMedian], na.rm=T) #A catch to get the median of raw data when doing covariate adjustment
    }
    #..............................................................................#
    # else {

      Data <- Data[order(Data$NewDate), ] # restores the time ordering
      
      TheSlopes <- GetInterObservationSlopes(Data,RawValues) # calculated slopes
      AnnualSenSlope <- median(TheSlopes$Slopes, na.rm=T) # the slopes are between individual observations that are incremenets X frequency apart in time.
      indexforMedian<-which(abs(TheSlopes$Slopes-AnnualSenSlope)==min(abs(TheSlopes$Slopes-AnnualSenSlope)))
      MySenCen<-as.character(unique(TheSlopes$CensorLabel[indexforMedian]))

      #Provide some warnings about Censored values used in the derivation of the Sen Slope
      if(!all(MySenCen == "not not")){
        if(all(MySenCen %in% c("lt lt","gt gt","lt gt","gt lt"))){
          AnalysisNote<-"WARNING: Sen slope based on two censored values"
        }else{
          AnalysisNote<-"WARNING: Sen slope influenced by censored values"}
      }else{ 
        if (AnnualSenSlope==0) AnalysisNote<- "WARNING: Sen slope based on tied non-censored values"
      }
      
      Data$CensoredOrig<-Data$Censored #Save teh original censoring for the plotting phase
      Data$CenTypeOrig<-Data$CenType
      Data$Censored <- FALSE # 
      Data$CenType <- "not"
      
      # estimate confidence intervals for 1-2*Alpha (where alpha = 0.05)
      VarS <- GetKendal(x = Data[, ])[["vars"]] # get the variance using the Kendall test function 
      Z <- 1-(0.05) # NB, (2*0.05/2) 2 X alpha but alpha/2 as per http://vsp.pnnl.gov/help/Vsample/Nonparametric_Estimate_of_Trend.htm 

      ########
      nC2 <-  length(TheSlopes$Slopes)
      RL <- (nC2 - qnorm(Z)*sqrt(VarS))/2 # Rank of lower  confidence limit 
      RU <- (nC2 + qnorm(Z)*sqrt(VarS))/2 # Rank of upper confidence limit  
      RankOfSlopes <- 1:nC2
      
      ConfInts <- approx(x=RankOfSlopes, y=sort(TheSlopes$Slopes), xout = c(RL, RU))$y
      Lci <- ifelse(is.na(ConfInts[1]),min(TheSlopes$Slopes),ConfInts[1])
      Uci <- ifelse(is.na(ConfInts[2]),max(TheSlopes$Slopes),ConfInts[2])

      # calculate the probability that the slope was truly below zero
      # rank of slope zero by interpolation
      R0 <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0,ties=median)$y
      
      #BUT if all slopes are either negative or positive, this fails, so catch
      if (sum(TheSlopes$Slopes<0)==length(TheSlopes$Slopes)){R0<-max(RankOfSlopes)}
      if (sum(TheSlopes$Slopes>0)==length(TheSlopes$Slopes)){R0<-min(RankOfSlopes)}
      
      Z1minusAlpha <- (2*R0 - nC2)/sqrt(VarS) 
      R0max <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0,ties=max)$y 
      Z1minusAlphamax <- (2*R0max - nC2)/sqrt(VarS) 
      R0min <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0,ties=min)$y 
      Z1minusAlphamin <- (2*R0min - nC2)/sqrt(VarS) 
      
      # The probability of slope being less than zero is
      Sen_Probability <- pnorm(Z1minusAlpha)
      Sen_Probabilitymax <- pnorm(Z1minusAlphamax)
      Sen_Probabilitymin <- pnorm(Z1minusAlphamin)
      
      # get intercept.
      Time <- (1:length(Data$V1))-1
      Ymed <-  median(Data$V1, na.rm=T) # the median of the measurements that are used to compute slope. 
      Tmed <- as.numeric(max(Data$NewDate)-min(Data$NewDate))/365.25/2# the median of the time
      Intercept <- Ymed - (AnnualSenSlope*Tmed)
      
      Percent.annual.change=AnnualSenSlope/abs(Median) * 100 #have to use abs, as sometimes median -ve after flow adjustment
      
  
      if(do.plot==TRUE) {
        if(RawValues==TRUE){ # We don't do this step if the data is FA - this is already done as part of the flow adjustement
          Data[Data$CenTypeOrig=="lt","V1"]<-Data[Data$CenTypeOrig=="lt","V1"]*0.5
          Data[Data$CenTypeOrig=="gt","V1"]<-Data[Data$CenTypeOrig=="gt","V1"]*1.1
        }
        
        myTimeIncr<-WhatistheIncr(x,"TimeIncr")
        mySeason<-ifelse(myTimeIncr=="Annual",NA,WhatistheIncr(x,"Season"))
        myplot<-GGPlotTrend(Data,x1,Intercept=Intercept,AnnualSenSlope=AnnualSenSlope,Lci=Lci,Uci=Uci,mymain=mymain,ValuesToUse = "V1",
                            IsSeasonal=FALSE,myTimeIncr=myTimeIncr,mySeason=mySeason,Ymed=Ymed,Tmed=Tmed,
                            Percent.annual.change=Percent.annual.change,UseMidObs=UseMidObs,...)
        
      }
      
    # } # end else
  } # end of first ifelse statement 
  
  output<-data.frame(Median=Median, Sen_VarS=VarS, AnnualSenSlope=AnnualSenSlope, Intercept=Intercept, 
                     Sen_Lci=Lci, Sen_Uci=Uci, AnalysisNote=AnalysisNote, Sen_Probability=Sen_Probability, Sen_Probabilitymax, Sen_Probabilitymin,
                     Percent.annual.change=Percent.annual.change)
  if(do.plot==TRUE&!is.na(output$AnnualSenSlope)){output<-list(output,myplot)}
  return(output)
} # endSenSlope function


#######################################################################
# Seasonal Kendall test.
#######################################################################
SeasonalKendall <- function(x, ValuesToUse = "RawValue", Year = "Year", RawValues=TRUE,
                            UseMidObs=TRUE,TimeIncrMed=TRUE, ...) {  
  # Calculates the seasonal Kendall test 
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the season (name =  "Season") and the year (name = "Year" )

  AnalysisNote <- GetAnalysisNote(x,ValuesToUse)   # CHECK TO SEE IF DATA OK TO CONTINUE
  
  if(AnalysisNote != "ok") { # if any TRUE dont do the test and return NA values
    VarS<-NA; S<-NA; D<-NA; tau<-NA;  Z<-NA; p<-NA; n<-length(nrow(x)) ;C=NA; Cd=NA;SeasIncr=NA
    } else { 

  #Tidy data and add on tine increment-year column
  x<-GetTimeIncrYear(x,ValuesToUse,Year)
  
  myTimeIncr<-WhatistheIncr(x,val="TimeIncr")
  mySeas<-WhatistheIncr(x,val="Season")
  x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])
  
  #Take medians over each time increment
  Data<-ValueForTimeIncr(x,ValuesToUse,Year,UseMidObs,TimeIncrMed)
  # CHECK TO SEE IF DATA OK TO CONTINUE    
  AnalysisNote <- GetAnalysisNote(Data,ValuesToUse = "V1",IsSeasonal = TRUE,SecondTierTest = TRUE)
 
  if(AnalysisNote != "ok") { # dont do the test if not sufficient variation, non censored values or long runs
    VarS<-NA; S<-NA;  D<-NA;  tau<-NA;  Z<-NA;  p<-NA;  n<-nrow(Data);C=NA;Cd=NA;
  } else {
  

    Data <- Data[order(Data$NewDate), ] # restores the time ordering

    # take each season, and compute kendall statistics 
    thisSeason <- by(Data, Data$Season, function(y) {  # for each season y = Data[Data$Season == "Jan",]
      GetKendal(x = y)
    })
    
    # sum kendall statistics over seasons
    S <- sum(sapply(thisSeason, function(r) return(r[["S"]])), na.rm=T) # nb put na remove here - some seasons do not have enough replicates to estimate S and return NaN from GetKendall
    VarS <- sum(sapply(thisSeason, function(r) return(r[["vars"]])))
    D <- sum(sapply(thisSeason, function(r) return(r[["D"]])))
    tau <- S/D
    n <- nrow(Data) # total observations
    
    if (n >= 10 & !is.na(VarS)&VarS>0) {
      SigmaS <- VarS^0.5        
      if(S > 0)  Z <- (S-1)/SigmaS
      if(S == 0) Z <- 0
      if(S < 0)  Z <- (S+1)/SigmaS
      if(Z > 0)  p <- pnorm(Z, lower.tail = F)*2
      if(Z == 0) p <- 1
      if(Z < 0)  p <- pnorm(Z, lower.tail = T)*2
      
      C<-1-p/2
      Cd<-C
      Cd[which(S>0)]<-p[which(S>0)]/2
    } else {
      AnalysisNote="Insufficent data to complete Seasonal Mann Kendall"
      Z <- NA
      p <- NA
      C<-NA
      Cd<-NA
    }
  } # end of second if-else
} # end of first if-else  
  prop.censored <- sum(x$Censored)/nrow(x) # proportion of values that were censored
  prop.unique <- length(unique(x[, ValuesToUse]))/nrow(x) # proportion of values that were censored
 no.censorlevels<-length(unique(x[x$Censored=="TRUE",ValuesToUse]))
  return(data.frame(nObs=nrow(x),nTimeIncr=n, S=S, VarS=VarS, D=D,tau=tau, Z=Z, p=p,C=C, Cd=Cd,
                    AnalysisNote=AnalysisNote,prop.censored=prop.censored,prop.unique=prop.unique,no.censorlevels=no.censorlevels,
                    TimeIncr=myTimeIncr,SeasIncr=mySeas))
}

#######################################################################
#Seasonal Sen slope.
#######################################################################
SeasonalSenSlope <- function(x, ValuesToUse = "RawValue",ValuesToUseforMedian=NA, Year = "Year",
                             do.plot=FALSE, mymain=NULL,RawValues=TRUE,UseMidObs=TRUE,
                             TimeIncrMed=TRUE,...) {  
  # Calculates the seasonal sen slope and its 90% CIs - for 95% confidence in trend direction. 
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the season (name "Season") and the year  (name "Year")
  
  AnalysisNote <- GetAnalysisNote(x, ValuesToUse = ValuesToUse,...)   # CHECK TO SEE IF DATA OK TO CONTINUE
  if(AnalysisNote != "ok") { # if not "ok" dont do the test 
    Median<-NA; VarS<-NA;   AnnualSenSlope<-NA;  Intercept<-NA;  Lci<-NA;  Uci<-NA;  
    Sen_Probability<-NA;  Sen_Probabilitymax <- NA;   Sen_Probabilitymin <- NA;  Percent.annual.change<-NA
    myplot<-NA; 
  } else {
    
    #Tidy data and add on timeincrement - Year column
    x<-GetTimeIncrYear(x,ValuesToUse,Year)
    myFreq<-x$myFreq[1]; x$myFreq<-NULL
    x1<-x ; x1$V1<-x1[,ValuesToUse]#Saved here to use for plotting later
    
    #Take medians or mid-values over each time increment
    Data<-ValueForTimeIncr(x,ValuesToUse,Year,UseMidObs,TimeIncrMed)
    
    if(is.na(ValuesToUseforMedian)){
      Median <- median(Data$V1, na.rm=T) #The median of the data after summarising to time increments
    }else{
      Median <- median(x[,  ValuesToUseforMedian], na.rm=T) #A catch to get the median of raw data when doing covariate adjustment
    }
    
    
    Data$Season<-factor(as.character(Data[, "Season"]),levels=unique(Data[,"Season"])) #Adjust the factors for the seasons for the analysis to allow completely missing seasons

      Data <- Data[order(Data$NewDate), ] # restores the time ordering
      Data$myDate<-Data$NewDate
      
      TheSlopes <- ddply(Data, c("Season"), function(y) GetInterObservationSlopes(y,RawValues))
      AnnualSenSlope <- median(TheSlopes$Slopes, na.rm=T) # the slopes are between individual observations that are incremenets X frequency apart in time.
      indexforMedian<-which(abs(TheSlopes$Slopes-AnnualSenSlope)==min(abs(TheSlopes$Slopes-AnnualSenSlope)))
      MySenCen<-as.character(unique(TheSlopes$CensorLabel[indexforMedian]))

      #Provide some warnings about Censored values used in the derivation of the Sen Slope
      if(!all(MySenCen == "not not")){
        if(all(MySenCen  %in% c("lt lt","gt gt","gt lt","lt gt"))){
          AnalysisNote<-"WARNING: Sen slope based on two censored values"
        }else{
          AnalysisNote<-"WARNING: Sen slope influenced by censored values"}
      }else{ 
        if (AnnualSenSlope==0) AnalysisNote<- "WARNING: Sen slope based on tied non-censored values"
        }
      
      Data$CensoredOrig<-Data$Censored
      Data$CenTypeOrig<-Data$CenType
      Data$Censored <- FALSE # Note the calculation of VarS and Senslope does not include censored values, but these fields are required by GetKendall
      Data$CenType <- "not"
      
      # estimate confidence intervals for 1-2*Alpha (where alpha = 0.05)
      VarS <- SeasonalKendall(x=Data, ValuesToUse="V1")$VarS # get the variance using the Seasonal Kendall test function 
      Z <- 1-(0.05) # NB, (2*0.05/2) 2 X alpha but alpha/2 as per http://vsp.pnnl.gov/help/Vsample/Nonparametric_Estimate_of_Trend.htm 
      
      ########
      nC2 <-  length(TheSlopes$Slopes)
      RL <- (nC2 - qnorm(Z)*sqrt(VarS))/2       # Rank of lower  confidence limit 
      RU <- (nC2 + qnorm(Z)*sqrt(VarS))/2      # Rank of upper confidence limit  
      
      RankOfSlopes <- 1:nC2
      
      ConfInts <- approx(x=RankOfSlopes, y=sort(TheSlopes$Slopes), xout = c(RL, RU))$y
      Lci <- ifelse(is.na(ConfInts[1]),min(TheSlopes$Slopes),ConfInts[1])
      Uci <- ifelse(is.na(ConfInts[2]),max(TheSlopes$Slopes),ConfInts[2])
      
      # calculate the probability that the slope was truly below zero
      # rank of slope zero by interpolation
      R0 <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0,ties=median)$y 
      Z1minusAlpha <- (2*R0 - nC2)/sqrt(VarS) 
      #BUT if all slopes are either negative or positive, this fails, so put in a catch
      if (sum(TheSlopes$Slopes<0)==length(TheSlopes$Slopes)){R0<-max(RankOfSlopes)}
      if (sum(TheSlopes$Slopes>0)==length(TheSlopes$Slopes)){R0<-min(RankOfSlopes)}
      
      R0max <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0,ties=max)$y 
      Z1minusAlphamax <- (2*R0max - nC2)/sqrt(VarS) 
      R0min <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0,ties=min)$y 
      Z1minusAlphamin <- (2*R0min - nC2)/sqrt(VarS) 
      
      # The probability of slope being less than zero is
      Sen_Probability <- pnorm(Z1minusAlpha)
      Sen_Probabilitymax <- pnorm(Z1minusAlphamax)
      Sen_Probabilitymin <- pnorm(Z1minusAlphamin)
      
      # get intercept.
      Time<-(1:length(Data$V1))-1 # need all dates despite some being censored. 
      Ymed <-  median(Data$V1, na.rm=T) # the median of the measurements that are used to compute slope. 
      Tmed <- as.numeric(max(Data$NewDate)-min(Data$NewDate))/365.25/2  # the median of the time
      Intercept <- Ymed - (AnnualSenSlope*Tmed)
      
      Percent.annual.change = AnnualSenSlope/abs(Median)*100 #using abs of median as sometimes median is -ve after flow adjustment
      
 
      if(do.plot==TRUE) { 
        if(RawValues==TRUE){ # We don't do this step if the data is FA - this is already done as part of the flow adjustement
          Data[Data$CenTypeOrig=="lt","V1"]<-Data[Data$CenTypeOrig=="lt","V1"]*0.5
          Data[Data$CenTypeOrig=="gt","V1"]<-Data[Data$CenTypeOrig=="gt","V1"]*1.1
        }
        
        myTimeIncr<-WhatistheIncr(x,"TimeIncr")
        mySeason<-ifelse(myTimeIncr=="Annual",NA,WhatistheIncr(x,"Season"))
        myplot<-GGPlotTrend(Data,x1,Intercept=Intercept,AnnualSenSlope=AnnualSenSlope,Lci=Lci,Uci=Uci,mymain=mymain,ValuesToUse = "V1",
                            IsSeasonal=TRUE,myTimeIncr=myTimeIncr,mySeason=mySeason,Ymed=Ymed,Tmed=Tmed,Percent.annual.change=Percent.annual.change,UseMidObs=UseMidObs,...)
        
      }
    # } # end of second if-else statement 
  }  # end of first if-else statement
  
output<-data.frame(Median=Median, Sen_VarS=VarS, AnnualSenSlope=AnnualSenSlope, Intercept=Intercept, 
                   Sen_Lci=Lci, Sen_Uci=Uci, AnalysisNote=AnalysisNote,
                   Sen_Probability=Sen_Probability, Sen_Probabilitymax = Sen_Probabilitymax,  Sen_Probabilitymin = Sen_Probabilitymin, 
                   Percent.annual.change=Percent.annual.change)
if(do.plot==TRUE&!is.na(output$AnnualSenSlope)){output<-list(output,myplot)}
  return(output)
} # end of function

#######################################################################
# Kendall test for censored data.
#######################################################################
# this is the heart of the procedures - calculates Kendall statistic and variance 
# accounting for censored values.
# adapted from the cenken function in the NADA package
# The S statistic is calculated as described on page 228 of Helsel's (2012) book
# "Statistics for censored environmental data using MINITAB and R", John Wiley & Sons. 

GetKendal <- function (x) { # x=Data
  
  # make some adjustments to data to get Helsel to agree (more closely) with TimeTrends.
  # TT sets all greater thans to a value slightly higher than their face value (and makes them all the same), 
  # then switches off all censored values (all=0). This treats all greater thans as ties.
  if(sum(x$CenType == "gt") > 0) { 
    MaxGTvalue <- max(x[x$CenType == "gt", "V1"])
    MaxGTvalue <- MaxGTvalue + 0.1
    x[x$CenType == "gt", "V1"] <- MaxGTvalue
    x[x$CenType == "gt", "Censored"] <- FALSE 
  }
  
  
  # For less thans, TimeTrends checks to see if there are multiple less thans with measured values  
  # between them. If so TimeTrends uses just uses Helsel code. If there are no real values between 
  # them TT sets them all equal and slightly less than their face value (i.e., a < 0.005s the same value as <0.001s) 
  # so that they form 1 group of ties. Helsel would treat them as two groups of ties.
  # if(sum(x$CenType == "lt") > 0) {
  #  TabLT <- as.data.frame.matrix(table(x$CenType, x$V1))
  #  ind <- TabLT["not", ] <= TabLT["lt", ] # what censored value has real values that are less?
  #  ind <- as.vector(ind)
  #  NoRealLess <- rle(ind)$lengths[1] # run length of index for which there are no real values less than;
  #  if(NoRealLess > 1) { # if first run>1, then tied less thans if true then set equal
  #    MaxCenVal <- as.numeric(names(TabLT)[NoRealLess])
  #    x$V1[x$CenType == "lt" & x$V1 <= MaxCenVal] <- MaxCenVal - 0.1*MaxCenVal
  #  }}
  
  x$SeasYr<-apply(x,1,function(y) paste0(y['Season'],y['Year']))
  SeasYrRank<-data.frame(SeasYr=unique(x$SeasYr))
  SeasYrRank$RankTime<-1:nrow(SeasYrRank)
  x$TimeRank<-SeasYrRank$RankTime[match(x$SeasYr,SeasYrRank$SeasYr)]
  
  xx <- x$V1
  cx <- x$Censored
  yy <- x$TimeRank # y is time and only needs to be  a sequentially increasing value
  cy <- rep(F, length.out = length(xx)) # no censored values for time
  n  <- length(xx)
  
  delx <- min(diff(sort(unique(xx))))/1000
  dely <- min(diff(sort(unique(yy))))/1000
  dupx <- xx - delx * cx
  diffx <- outer(dupx, dupx, "-")
  diffcx <- outer(cx, cx, "-")
  xplus <- outer(cx, -cx, "-")
  dupy <- yy - dely * cy
  diffy <- outer(dupy, dupy, "-")
  diffcy <- outer(cy, cy, "-")
  yplus <- outer(cy, -cy, "-")
  signyx <- sign(diffy * diffx)
  tt <- (sum(1 - abs(sign(diffx))) - n)/2
  uu <- (sum(1 - abs(sign(diffy))) - n)/2
  cix <- sign(diffcx) * sign(diffx)
  cix <- ifelse(cix <= 0, 0, 1)
  tt <- tt + sum(cix)/2
  signyx <- signyx * (1 - cix)
  ciy <- sign(diffcy) * sign(diffy)
  ciy <- ifelse(ciy <= 0, 0, 1)
  uu <- uu + sum(ciy)/2
  signyx <- signyx * (1 - ciy)
  xplus <- ifelse(xplus <= 1, 0, 1)
  yplus <- ifelse(yplus <= 1, 0, 1)
  diffx <- abs(sign(diffx))
  diffy <- abs(sign(diffy))
  tplus <- xplus * diffx
  uplus <- yplus * diffy
  tt <- tt + sum(tplus)/2
  uu <- uu + sum(uplus)/2
  
  
  test<-signyx * (1 - xplus) * (1 - yplus)

  itot <- sum(signyx * (1 - xplus) * (1 - yplus))
  kenS <- itot/2
  tau <- (itot)/(n * (n - 1))
  
  J <- n * (n - 1)/2
  taub <- kenS/(sqrt(J - tt) * sqrt(J - uu))
  
  varS <- n * (n - 1) * (2 * n + 5)/18
  
  intg <- 1:n
  dupx <- xx - delx * cx
  dupy <- yy - dely * cy
  dorder <- order(dupx)
  dxx <- dupx[dorder]
  dcx <- cx[dorder]
  dorder <- order(dupy)
  dyy <- dupy[dorder]
  dcy <- cy[dorder]
  tmpx <- dxx - intg * (1 - dcx) * delx
  tmpy <- dyy - intg * (1 - dcy) * dely
  rxlng <- rle(rank(tmpx))$lengths
  nrxlng <- table(rxlng)
  rxlng <- as.integer(names(nrxlng))
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  rylng <- rle(rank(tmpy))$lengths
  nrylng <- table(rylng)
  rylng <- as.integer(names(nrylng))
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  delc <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * 
                                                        (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n-1)) 
  
  x4 <- nrxlng * (rxlng - 1)
  y4 <- nrylng * (rylng - 1)
  tmpx <- intg * dcx - 1
  tmpx <- ifelse(tmpx < 0, 0, tmpx)
  nrxlng <- sum(tmpx)
  rxlng <- 2
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  tmpy <- intg * dcy - 1
  tmpy <- ifelse(tmpy < 0, 0, tmpy)
  nrylng <- sum(tmpy)
  rylng <- 2
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  deluc <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * 
                                                         n * (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n - 1)) - (sum(x4) + sum(y4))
  
  dxx <- dxx - intg * dcx * delx
  dyy <- dyy - intg * dcy * dely
  rxlng <- rle(rank(dxx))$lengths
  nrxlng <- table(rxlng)
  rxlng <- as.integer(names(nrxlng))
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  rylng <- rle(rank(dyy))$lengths
  nrylng <- table(rylng)
  rylng <- as.integer(names(nrylng))
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  delu <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * 
                                                        (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n -  1))
  
  varS <- varS - delc - deluc - delu
  
  if (n >= 3 & !is.na(varS)&varS>0) {
    SigmaS <- varS^0.5        
    if(kenS > 0)  Z <- (kenS-1)/SigmaS
    if(kenS == 0) Z <- 0
    if(kenS < 0)  Z <- (kenS+1)/SigmaS
    
    if(Z > 0)  p <- pnorm(Z, lower.tail = F)*2
    if(Z == 0) p <- 1
    if(Z < 0)  p <- pnorm(Z, lower.tail = T)*2
  } else {
    Z <- NA
    p <- NA
  }
  
  return(list(nTimeIncr = n, S = kenS, vars=varS,  D = J, tau = kenS/J, Z=Z, p=p))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#######################################################################
#                   TREND AGGREGATION FUNCTIONS
#######################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#This function returns the mode 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getTAU <- function(x, obs,SiteNameCol="siteID",VarNameCol="analyte",ValuesToUse = "RawValue", Year = "Year", 
                   UseMidObs=T,DirType="Decrease",Reverse=c("CLAR","MCI"),excl=c("pH")) { # 
  x$siteID<-x[,SiteNameCol];obs$siteID<-obs[,SiteNameCol]
  x$analyte<-x[,VarNameCol];obs$analyte<-obs[,VarNameCol]
  
  x <- x[!is.na(x$C), ] # remove NA values from trend evaluations that result from NotAnalysed
  
  myObs <- obs[obs$analyte == x$analyte[1], ] # observation data pertaining to this variable
  sitesInCommon <- intersect(unique(myObs$siteID), x$siteID) # sites in common to the orginal data and trend evaluations
  
  x<-x[x$siteID %in% sitesInCommon,]
  myObs <- myObs[myObs$siteID %in% sitesInCommon, ]  # ensure sites are consistent and get observations for the common sites
  
  M <- nrow(x)                # the total number of sites
  
  #Determine the smallest frequency used in the trend assessements
  if(all(x$TimeIncr=="Annual")){
    Frequency<-"Annual"
  } else{Frequency <- SeasonLabs$SeasName[min(match(unique(x$TimeIncr),SeasonLabs$mySeas))]    }  #
  
  sites.all <- x$siteID # all the sites
  names(sites.all) <- sites.all
  
  ModalDirection <- sign(sum(sign(x$S)) + sum(x$S ==0)/2) 
  if(ModalDirection==0){ModalDirection=1}#Just a catch for when it is exaclty equal
  DT <- ifelse(ModalDirection==1, "Increasing", "Decreasing") # aggregate trend direction
  
  TAU <- (sum(sign(x$S) == ModalDirection) + sum(x$S ==0)/2)/M # aggregate trend strength (the proportion of trends in the modal direction)
  
  ProbModalDirection <- ddply(x, "siteID", function(y) ifelse(sign(y$S) == ModalDirection, y$C, 1-y$C))
  x$P_Modal  <- ProbModalDirection[match(x$siteID, ProbModalDirection$siteID),2]
  
  VarTAU <- 1/M^2 * sum(x$P_Modal * (1-x$P_Modal))          # the non-adjusted variance
  
  SumVarTAU <- sum(x$P_Modal * (1-x$P_Modal))               # this is the first term From Yue and Wang 2002 Equation 17a (sum of individual variances).This is the first term in manuscript Eqn 9.
  
  # Following From Yue and Wang 2002 Equation 17a. Calculate pairwise variance (note this is square matrix we need only lower triangle). This is manuscript Eqn 10 without the cross correlation term
  Cov <- outer(x$P_Modal*(1 - x$P_Modal), x$P_Modal*(1-x$P_Modal)) # nb. variance per site is Cd(1-Cd)
  dimnames(Cov)  <- list(x$siteID, x$siteID)
  
  
  ###########################################################################
  # the adjustment to the variance to account for covariance
  
  if(Frequency == "Annual") {  # only one season = annual
    cat("Annual variable", as.character(x$analyte[1]), "\n")
    # get correlation of the observations (this is the pck,k+l term in  Yue and Wang 2002 Equation 17a.)
    
    # Add on time increment and time incr- year for all sites
    myObs$TimeIncr<-myObs[,Year]
    myObs.TimeIncr<-ddply(myObs,.(siteID,analyte),function(d){
      d<-GetTimeIncrYear(x=d, ValuesToUse = ValuesToUse,Year=Year)
      d <- ValueForTimeIncr(x=d, ValuesToUse=ValuesToUse, Year=Year, UseMidObs=UseMidObs) 
      return(d)
    })
    
    AllDates<-myObs.TimeIncr[,c("TimeIncr","TimeIncrYear")]
    AllDates<-AllDates[!duplicated(AllDates),]
    
    AllDates<-merge(AllDates,myObs.TimeIncr,all=T)
    AA<- as.data.frame(pivot_wider(AllDates[,c("TimeIncrYear","siteID","V1")],names_from="siteID",values_from="V1"))
    rownames(AA)<-AA$TimeIncrYear
    AA$TimeIncrYear<-NULL
    
    CorMatrix <- cor(AA, use="pairwise.complete.obs")  # get correlation of the observations (this is the pck,k+l term in  Yue and Wang 2002 Equation 17a.). This is the first term in manuscript Eqn 10.
  } 
  #######################################################################
  
  if (Frequency != "Annual") { # calculate the pairwise-site observation correlations for data that has frequency higher than annual
    # assume that the data can have frequency as high as "Frequency" - get a time series at the Frequency time-step with missing months represented by NA values 
    
    cat(SeasonLabs$mySeas[SeasonLabs$SeasName==Frequency]," variable", as.character(x$analyte[1]), "\n")
    
    # Add on time increment and time incr- year for all sites
    myObs$TimeIncr<-myObs[,Frequency]
    myObs.TimeIncr<-ddply(myObs,.(siteID,analyte),function(d){
      d<-GetTimeIncrYear(x=d, ValuesToUse = ValuesToUse,Year=Year)
      d <- ValueForTimeIncr(x=d, ValuesToUse=ValuesToUse, Year=Year, UseMidObs=UseMidObs) 
      return(d)
    })
    
    AllDates<-myObs.TimeIncr[,c("TimeIncr","TimeIncrYear")]
    AllDates<-AllDates[!duplicated(AllDates),]
   
    AllDates<-merge(AllDates,myObs.TimeIncr,all=T)
    AA<- as.data.frame(pivot_wider(AllDates[,c("TimeIncrYear","siteID","V1")],names_from="siteID",values_from="V1"))
    rownames(AA)<-AA$TimeIncrYear
    AA$TimeIncrYear<-NULL
    
    CorMatrix <- cor(AA, use="pairwise.complete.obs")  # get correlation of the observations (this is the pck,k+l term in  Yue and Wang 2002 Equation 17a.). This is the first term in manuscript Eqn 10.
   
  } 
  #####################################################################
  
  CorMatrix <- CorMatrix[sites.all, sites.all]
  Cov <- Cov[sites.all, sites.all] # this ensure sites align for the multiplication step
  CovTerm <- sqrt(Cov) * CorMatrix  
  SumCovTerm <- sum(CovTerm[lower.tri(CovTerm)], na.rm=T) # this is the whole added adjustment term from Yue and Wang 2002 Equation 17a.
  # note need to use na.rm = T because it can occur that a site has zero variance (all observations same value or censored) [a warning is issued in this instance]
  
  CorrectedVarTAU <- 1/M^2 * (SumVarTAU + (2 * SumCovTerm)) # this is the whole of Yue and Wang 2002 Equation 17a. This is the whole of manuscript equation 9. 
  
  UncorrectedConfidence<-getCT(TAU=TAU,VarTAU = VarTAU)
  Confidence <- getCT(TAU=TAU, VarTAU=CorrectedVarTAU) # get direction confidence based on the adjusted variance
  
  ConfCat=AssignConfCat(Confidence,CatType = "Direction")
  UncorrectedConfCat=AssignConfCat(UncorrectedConfidence,CatType = "Direction")
  
  SummaryResults <- data.frame(M=M, TAU=TAU, VarTAU=VarTAU, CorrectedVarTAU=CorrectedVarTAU, DT = DT,
                              CT = Confidence,ConfCat=ConfCat, UncorrectedCT = UncorrectedConfidence,UncorrectedConfCat=UncorrectedConfCat) 
  ObservationCorrelations <- CorMatrix[lower.tri(CorMatrix)]
  
  if(DirType=="Improve"){
    if(x$analyte[1] %in% excl){
      SummaryResults[,]<-NA
    }
    SummaryResults$DT=ifelse(ModalDirection==1, "Degrading", "Improving") 
    if(x$analyte[1] %in% Reverse){
      SummaryResults$DT=ifelse(ModalDirection==0, "Degrading", "Improving") 
    }

  }
  
  return(list(SummaryResults=SummaryResults, ObservationCorrelations=ObservationCorrelations ))
} # end  function 

#Make histograms of the between site correlations for regional analysis
PlotRegionalCor<-function(x){
  ########################
  # cross correlations
  ########################
  ObservationCorrelations <- ldply(names(x), function(y) { # x=names(AdjustedPdValues1)[1]
    theData <- x[[y]]$ObservationCorrelations
    theData <- data.frame(analyte = y, Correlations = theData)
    return(theData)
  })
  
  ObservationCorrelations$analyte <- factor(ObservationCorrelations$analyte)
  
  
  p <- ggplot(ObservationCorrelations, aes(x=Correlations)) +
    facet_wrap( ~ analyte,drop=T) +
    geom_histogram(aes(y=100*after_stat(count)/tapply(after_stat(count),after_stat(PANEL),sum)[after_stat(PANEL)]))+geom_vline(xintercept = 0, col="red", lty=2)+
    xlab("Pearson correlation coefficient") + ylab("Relative frequency (%)") + theme_bw(base_size = 12); p 
 return(p) 
}

#Make summary plot of the regional trend confidence and direction
PlotRegionalConfDir<-function(x,mymain=NA,DirType="Decrease"){
  x<-x[!is.na(x$TAU),]
  TauPlot <- x[, c("analyte", "TAU","M", "DT", "CT", "ConfCat")]
  TauPlot$Status <- "Corrected"
  TauPlot2 <- x[, c("analyte", "TAU", "M", "DT", "UncorrectedCT", "UncorrectedConfCat")]
  TauPlot2$Status <- "Uncorrected"
  names(TauPlot2) <- names(TauPlot)
  
  TauPlot <- rbind.data.frame(TauPlot, TauPlot2)
  shape_leg<-c(24,25)
  if(DirType=="Decrease"){
  TauPlot$Direction <- factor(TauPlot$DT, levels = c("Increasing", "Decreasing"))
  names(shape_leg)<-c("Increasing", "Decreasing")
  }else{
    TauPlot$Direction <- factor(TauPlot$DT, levels = c("Degrading", "Improving"))
    names(shape_leg)<-c("Degrading", "Improving")
  }
  TauPlot$Status <- factor(TauPlot$Status, levels = c("Corrected", "Uncorrected"))
  TauPlot$X_lab <- paste0(TauPlot$analyte, " (", TauPlot$M, ")")
  # TauPlot$X_lab <- factor(TauPlot$X_lab, levels = levels(CatConfSitesStak$Strip)) # get order correct
  
  TauPlot$ConfCat <-factor(TauPlot$ConfCat ,  levels = c("As likely as not",
                                           "Likely",
                                           "Very likely",
                                           "Highly likely"))
  
  
 p1 <- ggplot(data=TauPlot[TauPlot$Status == "Corrected", ], 
                                                         aes(y=TAU,x=X_lab,shape=Direction,fill=ConfCat))+
    geom_point(size=7,show.legend = T)+scale_shape_manual(values=shape_leg)+ 
    scale_fill_viridis(discrete = T, direction = 1, drop=F) +
    theme(panel.background=element_rect(fill="white"),panel.border=element_rect(fill=NA,colour="grey70",linewidth=1.2), legend.position="bottom",
          legend.key= element_rect(fill="white"),axis.text.x=element_text(angle=45,vjust=0.5),
          plot.title=element_text(size=13,hjust=0.5), panel.grid.major.y = element_line(color = "grey80", linewidth = 0.5,linetype = 2),
          panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5,linetype = 2),legend.box="vertical") +  
    ylab(expression(Aggregate~trend~strength~"("~hat(Tau)~")")) + xlab("Variable")+
    guides(fill = guide_legend(title.position = "top", override.aes=list(shape = 21),order=0),shape=guide_legend(order=1))+ylim(c(0.5,1))+
    labs(fill=expression(Confidence~aggregate~trend~at~regional~level~"(C"^Tau~")"), shape=expression(Direction~"(D"^Tau~")"))
 if(!is.na(mymain)){p1<-p1+ggtitle(mymain)}
  
  return(p1)
}
#################################################################################
# function to get aggregate trend strength, direction and confidence in direction 
#################################################################################
getCT <- function(TAU, VarTAU) {   
  Z_05 <- (TAU-0.5)/sqrt(VarTAU)
  CT<-pnorm(Z_05)
  return(CT)
}


#######################################################################
#  Assign Categorical confidence trend is decreasing (IPCC defintions)
#######################################################################


AssignConfCat<-function(x,CatType="Improve",Reverse=c("CLAR","MCI"),NoImproveDir=c("pH"),nCat="Simple"){
  #CatType options are: Direction, Improve or Decrease
  #Reverse is a list of variables for which increasing trends indicate improvement (only relevant when using "Improve")
  #NoImproveDir is a list of variables for which increasing trends are not clearly assocated with improvement or degradation  (only relevant when using "Improve")
  #nCat options are Full (full IPCC) or Simple (simplified categories)
  
  if(CatType=="Direction"){  #Returns of confidence in trend direction
    ifelse(is.vector(x),P<-x,P<-x$C)
    
    if(nCat=="Full"){
      mybreaks=c(0.49,0.67, 0.9, 0.95, 0.99, 1.01)
      mycatlabels=c("As likely as not", "Likely", "Very likely", "Highly likely", "Extremely likely", "Virtually certain")
    }else{
      mybreaks=c(0.49,0.67, 0.9, 0.95, 1.01)
      mycatlabels=c("As likely as not","Likely", "Very likely", "Highly likely")
    }
    
  }else{
    P<-x$Cd
    if(nCat=="Full"){
      mybreaks =c(-0.01, 0.01, 0.05, 0.1, 0.33,0.67, 0.9, 0.95, 0.99, 1.01)
      #These breaks and labels are based on guidance in IPCC report
      mycatlabels = c("Exceptionally unlikely",   "Extremely unlikely","Very unlikely", "Unlikely","As likely as not",
                      "Likely", "Very likely","Extremely likely","Virtually certain")
    }else{
      mybreaks =c(-0.01,  0.1, 0.33,0.67, 0.9, 1.01)
    }
    if(CatType=="Improve"){ #Returns confidence that trend is 
      if(nCat=="Simple"){
        mycatlabels=c("Very likely degrading", "Likely degrading","Low confidence",
                      "Likely improving", "Very likely improving")
      }
      
      if(!is.na(Reverse[1])){P[x$analyte %in% Reverse]<-1-P[x$analyte %in% Reverse]}
      if(!is.na(NoImproveDir[1])){P[x$analyte %in% NoImproveDir]<-NA}
      
      
      
    }else if (CatType=="Decrease"){
      if(nCat=="Simple"){
        mycatlabels=c("Very likely increasing", "Likely increasing","Low confidence",
                      "Likely decreasing", "Very likely decreasing")
      }
    }else{
      #Catch in case a CatType that doesn't exist is specified
      stop("CatType must be specified as Direction, Improve, or Decrease")
    }
  }
  
  
  ConfCats <-cut(P, breaks =mybreaks,labels=mycatlabels)
  ConfCats <-as.character(ConfCats )
  ConfCats [is.na(ConfCats )]<-"Not Analysed"
  ConfCats <-factor(ConfCats ,  levels = c(mycatlabels,"Not Analysed"))
  return (ConfCats)
}  




