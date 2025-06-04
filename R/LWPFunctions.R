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
# Most Recent update: September 2018; changes to calcualation of Sen slopes
# and probability trend was decreasing
#######################################################################
#######################################################################
# Please read accompanying documentation: LWPTrendsHelp_2018_09_18.pdf and
# note the disclaimer therein for the use of these scripts
#
# See also companion example files:
#       RunLWPTrendsExample_Sept18.R
#       LWPTrends_SampleData.rdata
#
#######################################################################
#######################################################################
# Functions for analysis of water quality trends
# These function replicate, as closely as possible, the Mann-Kendall and Seasonal Kendall tests in the TimeTrends software.
# Differences between this code and TT occur due to differences in the handling of censored values when calculating S and its variance.
# This is due to the key function is "GetKendal", which computes the Kendall test for censored data.
# The function is adapted from the "cenken" function in the NADA package
# The function returns the Kendall statistic and variance as described by Helsel (2012) "Statistics for censored environmental data using MINITAB and R", John Wiley & Sons.

# Another  difference is these functions do not allow  analysis of more than one value in each time period. If there are more than one value in a
# time period, the median value is used.

# the main analysis functions are as follows:
# Impute.upper - provides imputed values for left Censored data
# Impute.lower - provides imputed values for right censored data
# RemoveAlphaDetect - removes non detect ("<" and ">") and produces a dataframe with a numeric value and a boolean defining if value is censored.
# GetFlowAdjustment - flow adjusts values
# InspectData - options ot view data as timeseries or a heat plot showing WQ data or censored occasions
# NonSeasonalTrendAnalysis - wrapper function that calls MannKendall and SenSlope
# SeasonalTrendAnalysis - wrapper function that calls SeasonalKendall and SeasonalSenSlope
# MannKendall - Mann-Kendall test
# SenSlope - Calculates the annual sen slope
# SeasonalKendall - Seasonal Kendall test.
# SeasonalSenSlope - Calculates the seasonal sen slope
#Additonal Trend Aggregation functions  (applied to trend outputs across one variable at many sites)
# AnalysePIT - provides proportion of improving trends (PIT Statistic)
# FaceValueCounter - counts numbers of increasing and decreasing sites
# ImprovementConfCat - Assigns a cateorical probability of improvment to each trend (based on IPCC categories)
# ImprovementConfCatLAWA - Assigns a cateorical probability of improvment to each trend (based on LAWA categories)


# Suggested procedure for trend analysis.
# 1. Evaluate the data for the time-period of interest. Assess if there is sufficient data using MakeTimeSeries, PlotData
# 2. determine the appropriate time-period increment (seasons) e.g. months or quarters and define Year and Season in the dataframe.
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
if(is.na(match("plyr", installed.packages()[,1]))) stop("these functions require that the R package plyr is installed. Install plyr and try again.")
require(plyr)
#require(gridExtra)
if(is.na(match("NADA", installed.packages()[,1]))) stop("these functions require that the R package NADA is installed. Install NADA and try again.")
require(NADA)

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
Impute.lower <- function(x, ValuesToUse = "RawValue", forwardT="log", reverseT="exp", doPlot=F) {
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
      myData$i1Values <- myData[, ValuesToUse]
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

    if (any(class(myros) == c("ros", "lm"))) {  # if a model is fitted then contiune

      SortValues <- sort(usedValues, na.last = T)  # sort order of values # na.last to preserve the NA values on sorting

      NewROS <- cbind.data.frame(SortValues,as.data.frame(myros))

      # the plot helps to show why imputed replacements for censored values are often greater than
      # non-censored values.
      if(doPlot)  {
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
RemoveAlphaDetect <-  function(x, AdjustDL = 1, AdjustUpper = 1) {   #  removes non detect ("<" and ">") and produces a dataframe with a numeric value and a boolean defining if value is censored.
  if(is.numeric(x)) xNumeric <- x   # check values may be already numeric
  if(is.factor(x)) xNumeric <- unfactor(x) # converts to numeric (but note the values < and > are returned as NA)
  if(is.character(x)) xNumeric <- as.numeric(x) # converts to numeric (but note the values < and > are returned as NA)
  isNonDetLT <- grep("<", x)  # which values are below detection
  isNonDetGT <- grep(">", x)  # which values are above detection
  DetectionLimLT <- as.numeric(sapply(as.character(x[isNonDetLT]), function(x) DougSplit(x, myPart=2, split= "<", fixed = T ))) # returns the numeric values
  DetectionLimGT <- as.numeric(sapply(as.character(x[isNonDetGT]), function(x) DougSplit(x, myPart=2, split= ">", fixed = T ))) # returns the numeric values
  xNumeric[isNonDetLT] <- DetectionLimLT * AdjustDL # /2 # replaces original below detection values with the detection limit.(or halve it??)
  xNumeric[isNonDetGT] <- DetectionLimGT * AdjustUpper # /2 # replaces original below detection values with the detection limit.(or halve it??)
  Censored <- rep(FALSE, times = length(x))
  Censored[isNonDetLT] <- TRUE
  Censored[isNonDetGT] <- TRUE
  CenType <-  rep("not", times = length(x))  # classification of the type of censoring
  CenType[isNonDetLT] <- "lt" # less than censored
  CenType[isNonDetGT] <- "gt" # greater than censored
  CenType <- factor(CenType, levels = c("gt", "lt", "not"))
  cat("*** It is safe to ignore warning message ***")

  Censored <- ifelse(Censored == "TRUE", T, F) # ensure this a binary
  return(data.frame(RawValue = xNumeric, Censored=Censored, CenType=CenType)) # raw values means NOT Flow-Adjusted and no < or > values
}

#######################################################################
# Add extra date information onto original dataset
#######################################################################
#This function adds on some extra date formats onto the dataframe for subsequent analysis
GetMoreDateInfo<-function(Data,firstMonth=1){
  if(firstMonth>1){
    yearshift<-as.numeric(strftime(as.Date(paste("2000",firstMonth,"1",sep="-")), format = "%j"))
    MyNumMonthString<-c(NumMonthString[firstMonth:12],NumMonthString[1:firstMonth-1])
    if(firstMonth>3){
      MyQtrString<-c(NumQtrString[(ceiling(firstMonth/3)):4],NumQtrString[1:(ceiling(firstMonth/3)-1)])
    }else{MyQtrString<-NumQtrString}
  }else{
    MyNumMonthString<-NumMonthString
    MyQtrString<-NumQtrString
    yearshift<-0
  }

  Data$Year <- as.numeric(format(Data$myDate, "%Y"))
  if(firstMonth!=1){Data$CustomYear <- as.numeric(format(Data$myDate + (yearshift), "%Y"))} # Customised start year
  Data$Month <- format(Data$myDate, "%b")    # abbreviated  months

  Data$Month <- factor(format(Data$myDate, "%b"), levels = MyNumMonthString)
  #Data$yrmon <- paste(Data$CustomYear,as.character(Data$Month), sep="-")

  Data$Qtr <- factor(quarters(Data$myDate),levels=MyQtrString)
  #Data$YrQtr   <- paste(Data$CustomYear, as.character(Data$Qtr), sep = "*")

  if (is.null(Data$npID[1])==F){
    Data$npID <- as.character(Data$npID)
    Data$sIDnpID <- paste(Data$sID, Data$npID, sep = "*")}

  return(Data)}

#######################################################################
# Define season strings
#######################################################################
SeasonString <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#
NumMonthString <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
NumQtrString <- c("Q1","Q2","Q3","Q4")
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
                         Covariate = "Flow", Span = c(0.5, 0.75), doPlot = F, plotpval=F,plotloglog=F,...) {  # myFrame <- StackFrame

  # make sure the data are in date order!
  myFrame <- myFrame[order(myFrame$myDate), ]
  OriginalValues <- myFrame[,ValuesToAdjust]  # NOTE this function ignores the censored values an treats all values at face value

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
  #OUT$npID<-myFrame[,"npID"]
  #OUT$sID<-myFrame[,"sID"]
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
    myFrame2 <- myFrame[myFrame$HasFlowAndData, ]                  # try Covariate adjustment
    myFrame2$XXValues <- myFrame2[,ValuesToAdjust]
    myFrame2$XXFlow <- myFrame2[,Covariate]


    hasZeroValues <- ifelse(sum(myFrame2$XXValues <= 0) >0 | sum(myFrame2$XXFlow <= 0)>0, T, F)  # are there any zero or negative Covariate or variable values
    ZerosHere <- myFrame2$XXValues == 0                                                         # keep an index of zero values
    myFrame2$XXValues[myFrame2$XXValues == 0] <-  1e-10                                         # make zeros equal to small values
    myFrame2$XXFlow[myFrame2$XXFlow <= 0] <- 1e-10                                         # make zeros equal to small values
    myx<-seq(from=min(myFrame2$XXFlow),to=max(myFrame2$XXFlow),length.out=100) #Xvalues for plotting

    # fit GAM
    if(any(method == "Gam")) {
      myGam <- gam(XXValues ~ s(XXFlow), data = myFrame2)
      concHat <- predict(myGam, newdata = myFrame2)
      Residuals <- myFrame2$XXValues - concHat
      FlowAdjvalues <- median(OriginalValues) + Residuals
      R2 <- round(100* (1 - (sum(residuals(myGam)^2)/sum((mean((myFrame2$XXValues)) - (myFrame2$XXValues))^2))))
      R<-cor.test(concHat,myFrame2$XXValues,method="pearson")$estimate
      p<-cor.test(concHat,myFrame2$XXValues,method="pearson")$p.value
      myy<-predict(myGam, newdata = data.frame(XXFlow=myx))
      AdjustList[["Gam"]] <- list(FlowAdjvalues = Residuals, Residuals = Residuals, concHat = concHat, Fitted = fitted(myGam), R2 = R2,R=R,p=p,method="Gam",myy=myy,myx=myx)
    }

    # fit log-log
    if(any(method == "LogLog")) {
      myLogLog <- lm(log10(XXValues) ~ log10(XXFlow), data = myFrame2)
      concHat <- 10^predict(myLogLog, newdata = myFrame2)
      Residuals <- myFrame2$XXValues - concHat
      FlowAdjvalues <- median(OriginalValues) + Residuals
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
          FlowAdjvalues <- median(OriginalValues) + Residuals
          R2 <- round(100 * (1 - (sum(residuals(Lowess)^2)/sum((mean((myFrame2$XXValues)) - (myFrame2$XXValues))^2))))
          R<-cor.test(concHat,myFrame2$XXValues,method="pearson")$estimate
          p<-cor.test(concHat,myFrame2$XXValues,method="pearson")$p.value
          myy<-predict(Lowess, newdata = data.frame(XXFlow=myx))
          AdjustList[[thisName]] <- list(FlowAdjvalues = FlowAdjvalues, Residuals = Residuals, concHat = concHat, Fitted = fitted(Lowess), R2 = R2,R=R,p=p,method=thisName,myy=myy,myx=myx)
        }
      }
    }


    if(doPlot) {
      if (plotloglog){par(mfrow = c(1,2), mar=c(4,5,1,1))}
      myCols <- c("red", "darkgreen", "cyan", "orange", "blue")
      myPch <-  c(4, 1, 16, 2, 3)

      myFits<-data.frame(do.call("rbind", lapply(AdjustList, function(df) cbind.data.frame(df$myx,df$myy,df$method))))
      names(myFits)<-c("myx","myy","Method")
      myFits$myx<-as.numeric(as.character(myFits$myx))
      myFits$myy<-as.numeric(as.character(myFits$myy))
      myMethods <- levels(myFits$Method)

      plot(x = myFrame2[!ZerosHere,"XXFlow"], y = myFrame2[!ZerosHere,"XXValues"], pch=16, cex=0.5, xlab="Flow", ylab = "Value", ...)
      for (i in 1:length(myMethods))  lines(x=myFits$myx[myFits$Method == myMethods[i]], y = myFits$myy[myFits$Method == myMethods[i]], col=myCols[i])

      if(plotpval) {

        ModR2_p <- ldply(AdjustList, function(x) {
          PVal <- "NOT SIGNIFICANT p>0.05"
          if(x$p < 0.05) PVal <- "SIGNIFICANT p<0.05"
          return(paste0("R2 =", round((x$R)^2*100), "% ", PVal))  #Note - R2 here is calcualted in teh original unit space
        })

        ModR2_p$lab <- paste(ModR2_p[,1], ModR2_p[,2])
        legend("top",  legend=ModR2_p$lab, col = myCols, text.col = "black", lty = c(1), pch = -1, bg = 'white', inset = .05, cex=0.8)

      }  else {# end plotpval
        legend("top",  legend=myMethods, col = myCols, text.col = "black",  lty = c(1), pch = -1, bg = 'white', inset = .05, cex=0.8)
      }

      if(plotloglog){
        plot(x = myFrame2[!ZerosHere,"XXFlow"], y = myFrame2[!ZerosHere,"XXValues"], pch=16, cex=0.5, xlab="Flow", ylab = "Value",log="xy", ...)
        for (i in 1:length(myMethods))  lines(x=myFits$myx[myFits$Method == myMethods[i]], y = myFits$myy[myFits$Method == myMethods[i]], col=myCols[i])

      }

      # myResiduals<-data.frame(do.call("rbind", lapply(AdjustList, function(df) cbind.data.frame(df$Fitted,df$Residuals,df$method))))
      # names(myResiduals)<-c("myFitted","myResiduals","Method")
      # myResiduals$myFitted<-as.numeric(as.character(myResiduals$myFitted))
      # myResiduals$myResiduals<-as.numeric(as.character(myResiduals$myResiduals))
      #
      # plot(x = myResiduals$myFitted, y = myResiduals$myResiduals, pch=myPch, cex=0.5, xlab="Fitted Values", ylab = "Residuals", col= "white", ...)
      # for (i in 1:length(myMethods))  points(x=myResiduals$myFitted[myResiduals$Method == myMethods[i]], y = myResiduals$myResiduals[myResiduals$Method == myMethods[i]], col=myCols[i], pch=myPch[i], cex=0.5)
      # abline(0,0, lty=3)
      # legend("top",  legend=myMethods, col = myCols, text.col = "black", pch = myPch, bg = 'white', inset = .05, cex=0.8)

    }

    for(i in 1:NoAdjusts) OUT[,i][myFrame$HasFlowAndData] <- AdjustList[[i]]$Residuals[order(myFrame$myDate)[myFrame$HasFlowAndData]] # put in original date order only for oginal data with value and flow.
    for(i in 1:NoAdjusts) OUT[,NoAdjusts+i][myFrame$HasFlowAndData] <- AdjustList[[i]]$R
    for(i in 1:NoAdjusts) OUT[,NoAdjusts*2+i][myFrame$HasFlowAndData] <- AdjustList[[i]]$p


    op <- options(warn=0) # set back to default warnings.
  }

  return(OUT)  # pass back the frame with all data including the Covariate adjustment
}
#######################################################################
# Inspect the data.
#######################################################################

InspectData <- function(x, plotType = c("TimeSeries", "Matrix"), PlotMat = c("RawData","Censored"),
                        Year = "Year",
                        StartYear = 1990, EndYear = 2015,
                        doPlot = TRUE, FlowtoUse=NULL,
                        ...) { # ... arguments passed to the plots

  # this converts the data to an object that represents a time series for the period specified by StartYear and EndYear and for the time increment set by seasons
  # the function adds NAs for gaps and takes the median value where there are multiple observations in a season
  # input data must have columns specifying the season (name =  "Season") and the year (name = "Year")
  # the SeasonString, must be specified in the environment
  # the function returns a summary dataframe and a plot if doPlot = TRUE
  if(is.null(x$Flow)&!is.null(FlowtoUse)) {x$Flow<-x[,FlowtoUse]}else{x$Flow<-NA}
  range(x$Year)

  myYears <-  StartYear:EndYear #
  YearsInPeriod <- length(myYears)
  SeasonString<-levels(x$Season)
  mySeasonYearFrame <- expand.grid(list(Season = SeasonString, Year = myYears))
  mySeasonYearFrame$SeasonYear <- paste(mySeasonYearFrame$Season, mySeasonYearFrame$Year, sep = "-")
  mySeasonYearFrame$SeasonYear <- factor(mySeasonYearFrame$SeasonYear, levels = mySeasonYearFrame$SeasonYear)

  x$SeasonYear <- paste(x[["Season"]], x[[Year]], sep = "-")
  x$SeasonYear <- factor(x$SeasonYear, levels = mySeasonYearFrame$SeasonYear)

  x <- x[x[[Year]] >= StartYear & x[[Year]] <= EndYear, ] # trims down the data to the year range

  # this takes the median of observations in Season (so a time increment only has one value)
  # if the median of the raw values is less than OR EQUAL to the max censored value, call the value censored, otherwise NOT censored.
  Data <-  ddply(x, "SeasonYear", function(y) {   # y=x[x$SeasonYear == "Nov-2009", ]
    Median <- median(y$RawValue)
    if(sum(y$Censored) == 0) {
      Censored = FALSE
    } else {
      MaxCensored <- max(y$RawValue[y$Censored],na.rm = T)
      Censored <- ifelse(Median <= MaxCensored, TRUE, FALSE)
    }
    return(data.frame(V1 = Median, Censored = Censored,Flow<-median(y$Flow)))
  })
  # this determines censor type for Data  based on CenType that has highest prevalence in the season
  Data$CenType <- ddply(x, "SeasonYear", function(y) names(which.max(table(y$CenType))))[,2]

  PlotData <- merge(mySeasonYearFrame, Data, by="SeasonYear", all=TRUE) # this over-writes to include all years and seasons specified by StartYear & EndYear
  PlotData$Row <- 1:nrow(PlotData)

  if(plotType == "TimeSeries") {    # plots the time series object as a time series
    dateseq<-seq(from=as.Date(paste0(StartYear,"-01-01")),to=as.Date(paste0(EndYear,"-12-31")),by="month")
    plot(x=x$myDate,y=x$RawValue, xaxt="n", xlab = "Date", ylab = "Value", pch=4,...)
    axis.Date(1, at=dateseq, labels=format(dateseq, "%b-%Y"))
    points(x= x$myDate[!is.na(x$Censored) & x$Censored],y= x$RawValue[!is.na(x$Censored) & x$Censored], col="red", pch=16)
    legend("topleft",  legend=c("Observation","Censored"), col = c("black","red"), text.col = "black",  pch = c(4, 16), bg = 'white', inset = .05)

  } else  {      # plots the time series object as a heat map
    myFreq <- length(SeasonString)
    nYears <- nrow(PlotData)/myFreq
    if(PlotMat == "RawData") {
      myMat <- ((matrix(PlotData$V1, nrow=myFreq, ncol=nYears, dimnames = list(SeasonString, myYears))))
    } else {
      myMat <- ((matrix(PlotData$Censored, nrow=myFreq, ncol=nYears, dimnames = list(SeasonString, myYears))))
    }
    NiceImage(myMat,myYlab=Year, ...)
  }
  # assess the data
  hasData <- !is.na(PlotData$V1) # dates with data values

  nOccasions <-  nrow(PlotData)
  nData <- sum(hasData)             # NB, this is the unique years and season combinations.
  nMissing <- nOccasions - nData
  firstYear <- min(PlotData$Year[hasData], na.rm=T) # first year with data
  lastYear  <- max(PlotData$Year[hasData], na.rm=T) # last year with data
  nyear <- countUnique(PlotData$Year[hasData]) # number of actual years with data
  propCen <- sum(PlotData$Censored[hasData])/nData
  if(!is.null(PlotData$Flow)) {
    nFlow <- sum(!is.na(PlotData$Flow[hasData]))
  } else {  nFlow = 0 }

  return(data.frame(nOccasions=nOccasions, YearsInPeriod=YearsInPeriod, nData = nData, nMissing = nMissing, firstYear = firstYear, lastYear = lastYear, nYears = nyear, propCen = propCen, nFlow = nFlow))
}

#Generates heat plots for Inspect data
NiceImage <- function(myMat, myPrint = F, yTextCex = 1, myYlab="Years",...) {
  par(xpd = F)
  image(myMat, x = 0:nrow(myMat), y = 0:ncol(myMat), axes = F, ylab=myYlab, xlab="Seasons", ...)
  abline(v = seq(0,nrow(myMat), length = nrow(myMat) + 1) )
  abline(h = seq(0,ncol(myMat), length = ncol(myMat) + 1) )
  par(xpd = T)
  myX <- seq(0.5, nrow(myMat)-0.5, by = 1)
  myY <- seq(0.5, ncol(myMat)-0.5, by = 1)
  text(x = myX, y = -0.2, dimnames(myMat)[[1]], pos = 2, srt = 90)
  text(x = 0, y = myY, dimnames(myMat)[[2]], pos = 2, cex = yTextCex)#
  if(myPrint) {
    for(n in 1:nrow(myMat)) {
      for(m in 1:ncol(myMat)) {
        text(n-0.5, m-0.5, round(myMat[n,m], digits = 2), cex = 0.75)
      }
    }
  }
  return(NULL)
}



#######################################################################
# Trend analysis utility functions
#######################################################################
#These pieces of code were previously imbedded within the Kendall and Sen
# functions, but were repeated.  They have been extracted in order to ensure
# that there is consistency in their application if/when updates are required

ApplyHiCensor <-function(x, ValuesToUse = "RawValue",HiCensor=TRUE){
  if (is.logical(HiCensor)){
    CenMax <- max(x[, ValuesToUse][x$CenType == "lt"], na.rm=T)
  }else{CenMax=HiCensor}

  x$Censored[x[, ValuesToUse] < CenMax] <- TRUE
  x$CenType[x[, ValuesToUse] < CenMax] <- "lt"  # and treat these as all censored.
  x[, ValuesToUse][x[, ValuesToUse] < CenMax] <- CenMax # set all values below the max censored value, to the max censored
  return(x)
}

GetSeasonYear<-function(x,ValuesToUse = "RawValue",Year="Year"){

  # make sure the data are in date order!
  x <- x[order(x$myDate), ]

  if(x$Season[1]==x[1,Year]){
    x$SeasonYear<-x$Season
    x$myFreq<-1
  }else{
    x$SeasonYear <- paste(x[["Season"]], x[[Year]], sep = "-")
    x$myFreq<-length(unique(x$Season))
  }

  # remove values and censored values that are NA.
  take <- !is.na(x[, ValuesToUse]) & !is.na(x$Censored)
  x <- x[take, ] # drop the NA values

  return(x)
}

MedianForSeason <-function(x, ValuesToUse = "RawValue",Year="Year"){

  # this takes the median of observations and DATES in Season (so a time increment only has one value)
  # if the median of the raw values is less than OR EQUAL to the max censored value, call the value censored, otherwise NOT censored.

  Data <-  ddply(x, "SeasonYear", function(y) {
    Median <- median(y[, ValuesToUse]) # get the median of the values
    NewDate <- median(y[, "myDate"]) # get the median of the dates
    if(sum(y$Censored) == 0) {
      Censored = FALSE
    } else {
      MaxCensored <- max(y[, ValuesToUse][y$Censored],na.rm = T)
      Censored <- ifelse(Median <= MaxCensored, TRUE, FALSE)
    }
    return(data.frame(V1 = Median, NewDate = NewDate, Censored = Censored, myYear=y[1,Year], Season=y$Season[1]))
  })
  # this determines censor type for Data  based on CenType that has highest prevalence in the season
  Data$CenType <- ddply(x, "SeasonYear", function(y) names(which.max(table(y$CenType))))[,2]# y = x[x$myYearSeason == "2015-Nov", ]

  return(Data)
}

GetAnalysisNote<-function(Data, ValuesToUse = "RawValue",IsSeasonal=FALSE,SecondTierTest=FALSE,...){
  #This function is used to check that there are sufficient data at several points in teh analysis
  #to contintue with the calcualtions, and if there are not to provide an error message describing the issue
  AnalysisNote="ok"
  #Set generic levels here about nubmers of unique values and non-cenosred values required
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
            longRun <- max(unlist(rle(diff(theValues))[1]))/length(theValues) > 0.5  # this catches seasons with one value to avoid warning.
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
  TrendClass <- "Indeterminant" # default assumption is  that there is insufficient data to define trend direction
  if(A3$Probability >= 0.95)  TrendClass <- "Decreasing"
  if(A3$Probability <= 0.05)  TrendClass <- "Increasing"
  if(is.na(A3$Probability))  TrendClass <- "Not Analysed"

  TrendDirection<- "Indeterminant" # rare - but dooes occur!
  if(A3$S < 0 )  TrendDirection <- "Decreasing"
  if(A3$S > 0)  TrendDirection <- "Increasing"
  if(is.na(A3$S))  TrendDirection <- "Not Analysed"
  return(data.frame(TrendCategory=TrendClass,TrendDirection=TrendDirection))
}

GetInterObservationSlopes<-function(Data,myPrecision=NULL){
  # number of Time Increments between observations, this excludes comparing observations in the same year+season
  AllTimeIncs <- outer(as.numeric(Data$NewDate), as.numeric(Data$NewDate), `-`)/365.25 # the time increment in years
  SameYearSeason <- outer(as.numeric(Data$NewDate), as.numeric(Data$NewDate), "==")
  AllTimeIncs[SameYearSeason] <- NA # remove any increments withing the same year + season

  # take each observation and compute differences for all other observations
  AllDifferences <- outer(Data$V1, Data$V1, `-`)
  AllDifferences[SameYearSeason] <- NA # remove any distances withing the same year + season

  CenLab <- outer(Data$CenType, Data$CenType, 'paste')
  Slopes <- AllDifferences/AllTimeIncs

  OUTPUTS <- data.frame(Slopes=as.vector(Slopes [lower.tri(Slopes, diag = FALSE)]),
                        CensorLabel=as.vector(CenLab [lower.tri(CenLab, diag = FALSE)]) )

  OUTPUTS <- OUTPUTS[!is.na(OUTPUTS$Slopes),]
  return(OUTPUTS)
}

PlotTrend<-function(x,x1,mymain=NULL,Intercept=NA,ValuesToUse = "RawValue",AnnualSenSlope=NA,
                    Lci=NA,Uci=NA,IsSeasonal=FALSE,Ymed=NA,Tmed=NA,Percent.annual.change=NA,Probability=NA,legend.pos="top"){
  if(is.null(x$myDate)) x$myDate<-x$NewDate

  Ymed<-median(x[,ValuesToUse])
  Interceptu <- Ymed - (Uci*Tmed)
  Interceptl <-Ymed - (Lci*Tmed)

  ylimlow <- min(0,min(min(x[, ValuesToUse], na.rm=T),min(x1[, ValuesToUse], na.rm=T)))
  ylimhi <- max(max(x[, ValuesToUse], na.rm=T),max(x1[, ValuesToUse], na.rm=T))*1.1
  if(ylimlow<0){ylimlow<-ylimlow-0.05*(ylimhi-ylimlow)}

  plot(x1$myDate, x1[,ValuesToUse], ylab = "Value",  xlab = "Time", ylim=c(ylimlow,ylimhi),
       pch=21,col="gray68",cex=1.2) #
  points(x1$myDate[x1$Censored], x1[, ValuesToUse][x1$Censored], col="tomato1", pch=21,cex=1.2)
  grid(nx = NULL, ny = NULL, col = "lightgray")
  # value at end of time period
  T1 <- Intercept + as.numeric(diff(range(x$myDate)))/365.25 * AnnualSenSlope
  T1l <- Interceptl + as.numeric(diff(range(x$myDate)))/365.25 * Lci
  T1u <- Interceptu + as.numeric(diff(range(x$myDate)))/365.25 * Uci

  segments(x0=min(x$myDate), y0 = Intercept, x1 = max(x$myDate), y1=T1, col="blue", lwd=2)
  segments(x0=min(x$myDate), y0 = Interceptu, x1 = max(x$myDate), y1=T1u, col="dodgerblue", lwd=2,lty=2)
  segments(x0=min(x$myDate), y0 = Interceptl, x1 = max(x$myDate), y1=T1l, col="dodgerblue", lwd=2,lty=2)

  points(x$myDate, x[, ValuesToUse], col="black", pch=16)
  points(x$myDate[x$CensoredOrig], x[, ValuesToUse][x$CensoredOrig], col="red", pch=16)

  SenLabel<-ifelse(IsSeasonal,"Annual Seasonal Sen Slope","Annual Sen Slope")
  if(is.null(mymain)==FALSE){title(main=mymain)}


  legend(paste0(legend.pos,"left"),  legend=c(paste("% ",SenLabel," = ", round(Percent.annual.change,1),"%"),
                                              paste(SenLabel," = ", signif(AnnualSenSlope,3)),
                                              paste("Probability decreasing = ",round(Probability,3))),
         text.col = "black",  pch = c(NA, NA), bg = 'white', inset = .05,cex=0.8)

  legend(paste0(legend.pos,"right"),  legend=c("Trend","90% C.I.", "Non-censored","Censored","Raw Observations"),
         col = c("blue","dodgerblue","black","red","gray68"),  text.col = "black",
         lty = c(1,2, -1,-1,-1), pch = c(-1, -1,16, 16,21), bg = 'white', inset = .05,cex=0.8)

}
#######################################################################
#Wrapper functions
#######################################################################
NonSeasonalTrendAnalysis<-function(...){

  A1<-MannKendall(...)
  if(A1$AnalysisNote =="ok"){#Then carry on and do the Sen Slope test
    A1$AnalysisNote <- NULL
    A2<-SenSlope(Probability=A1$Probability,...)
  }else{
    A2<-data.frame(Median=NA, Sen_VarS=NA, AnnualSenSlope=NA, Intercept=NA,
                   Lci=NA, Uci=NA,Sen_Probability=NA, Probabilitymax=NA, Probabilitymin=NA,
                   Percent.annual.change=NA)
  }
  A<-cbind.data.frame(A1,A2)

  #Add on trend category and trend direction
  if(is.na(A$Probability)){
    A$TrendCategory<-"Not Analysed"
    A$TrendDirection<-"Not Analysed"
  }else{
    A<-cbind.data.frame(A,GetTrendDirectionandClass(A))}

  return(A)

}

SeasonalTrendAnalysis<-function(...){
  #Do Seasonal Kendall Test
  A1<-SeasonalKendall(...)

  if(A1$AnalysisNote =="ok"){
    #Contiue on to do Seasonal Sen Slope if Seasonal Kendall Test was sucessful
    A1$AnalysisNote <- NULL
    A2<-SeasonalSenSlope(Probability=A1$Probability,...)

  }else{
    A2<-data.frame(Median=NA, Sen_VarS=NA, AnnualSenSlope=NA, Intercept=NA,
                   Lci=NA, Uci=NA,Sen_Probability=NA, Probabilitymax=NA, Probabilitymin=NA,
                   Percent.annual.change=NA)
  }
  A<-cbind.data.frame(A1,A2)
  if(is.na(A$Probability)){
    A$TrendCategory<-"Not Analysed"
    A$TrendDirection<-"Not Analysed"
  }else{
    A<-cbind.data.frame(A,GetTrendDirectionandClass(A))}

  return(A)
}


#######################################################################
#Seasonality Test
#######################################################################
SeasonalityTest <-  function(x, ValuesToUse = "RawValue", HiCensor=FALSE,do.plot=TRUE,...) { # x = DataForTest

  #Following the treatment suggested by Helsel (2012), all censored values and values less than
  #the highest non-detect (<) are assigned the same low value and treated as ties.
  #If the censored values are > then all non-detects and values higher than the
  #lowest non-detect (>) are assigned the same high value and are treated as ties.

  Observations <- sum(!is.na(x[, ValuesToUse]))

  if(!is.factor(x[, "Season"])) x[, "Season"] <- factor(as.character(x[, "Season"]), levels = SeasonString) # Season needs to be a factor - put in correct order

  if(Observations > 5 & length(unique(x[, ValuesToUse]))>1) { # check there is minimum data and it has some variation

    #If required, apply the hi-censor
    if(HiCensor) x<-ApplyHiCensor(x,ValuesToUse,HiCensor)

    if(do.plot==TRUE){boxplot(x[, ValuesToUse] ~ x[, "Season"], col= "lightblue", range=0, ...)}

    if(length(table(x[, "Season"][!is.na(x[, ValuesToUse])]))>1)
    { # check there is more than one season represented
      #if(!is.factor(x$Season))  x$Season<-factor(as.character(x[, "Season"]),levels=unique(x[,"Season"]))

      KW <- kruskal.test(x[, ValuesToUse] ~ factor((x[, "Season"])))
      OUT <- data.frame(Observations=Observations, KWstat = KW$statistic, pvalue = KW$p.value)

      if(do.plot==TRUE){legend("topleft",  legend=c(paste("Observations =", Observations),   paste("Kruskal Wallis statistic = ", round(OUT$KWstat,3)), paste("P value =", round(OUT$pvalue,3))),
                               text.col = "black",  pch = c(NA, NA), bg = 'white', inset = .05)}
    } else {
      OUT <- data.frame(Observations=Observations, KWstat = NA, pvalue = NA)
    }

  } else {
    OUT <- data.frame(Observations=Observations, KWstat = NA, pvalue = NA)
  }
  return(OUT)
}

#######################################################################
# Mann-Kendall test.
#######################################################################
MannKendall <- function(x, ValuesToUse = "RawValue", Year = "Year", HiCensor=FALSE,...) { #
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have Observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the season (name =  "Season") and the year (name = "Year" )

  # CHECK TO SEE IF DATA IS OK TO CONTINUE
  AnalysisNote <- GetAnalysisNote(x,ValuesToUse)

  if(AnalysisNote !="ok") { # if any TRUE dont do the test
    KendalTest <- data.frame(nObs = nrow(x), VarS=NA, S = NA, D = NA, tau = NA, Z=NA, p=NA,Probability=NA)
  } else {

    #..............................................................................#
    #Tidy data and add on SeasonYear column
    x<-GetSeasonYear(x,ValuesToUse,Year)

    #if the HiCensor option is true and there are censored values
    if(HiCensor!=FALSE & sum(x$CenType == "lt") > 0)   x<-ApplyHiCensor(x,ValuesToUse,HiCensor)

    x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])

    #Take medians over each season
    Data<-MedianForSeason(x,ValuesToUse,Year)
    #..............................................................................#

    # CHECK TO SEE IF DATA OK TO CONTINUE
    AnalysisNote <- GetAnalysisNote(Data,ValuesToUse="V1",SecondTierTest = TRUE)

    if(AnalysisNote != "ok" ) { # if any TRUE dont do the test
      KendalTest <- data.frame(nObs = nrow(x), VarS=NA, S = NA, D = NA, tau = NA, Z=NA, p=NA,Probability=NA)
    } else {

      # organise the ordering of data
      Data <- Data[order(Data$NewDate), ] # restores the time ordering

      KendalTest <- GetKendal(x = Data)
      KendalTest$Probability<-1-0.5*KendalTest$p
      KendalTest$Probability[which(KendalTest$S>0)]<-0.5*KendalTest$p[which(KendalTest$S>0)]
      names(KendalTest)[names(KendalTest)=="vars"]<-"VarS"
      KendalTest<-as.data.frame(KendalTest)
      if(is.na(KendalTest$Probability)==T) AnalysisNote = "Not Analysed"
    }
  }  # end of first if statement

  KendalTest$AnalysisNote<-AnalysisNote
  KendalTest$prop.censored <- sum(x$Censored)/nrow(x) # proportion of values that were censored
  KendalTest$prop.unique <- length(unique(x[, ValuesToUse]))/nrow(x) # proportion of values that were censored
  KendalTest$no.censorlevels<-length(unique(x[x$Censored=="TRUE",ValuesToUse]))

  return(KendalTest)
}

#######################################################################
# Sen slope test
#######################################################################
SenSlope <- function(x, ValuesToUse = "RawValue", ValuesToUseforMedian="RawValue",HiCensor=F, Year = "Year",
                     doPlot=FALSE,mymain=NULL, ...) {  # x=DataForTest
  # Calculates the annual Sen slope and its 90% CIs - for 95% confidence in trend direction.
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the season (name =  "Season") and the year (name = "Year" )

  # CHECK TO SEE IF DATA OK TO CONTINUE
  AnalysisNote <- GetAnalysisNote(x, ValuesToUse = ValuesToUse,...)

  if(AnalysisNote != "ok") { # if any TRUE dont do the test and return NA values
    Median<-NA; VarS<-NA; AnnualSenSlope<-NA; Intercept<-NA; Lci<-NA; Uci<-NA;
    Sen_Probability<-NA; Probabilitymin<-NA;Probabilitymax<-NA; Percent.annual.change<-NA
  } else {
    #..............................................................................#
    # For the purposes of of computing the sen slope, replace "lt" censored values with half
    # the censored values and "gt" values with 1.1 the censored value (to aoiv ties with non-cenosred obs)
    x[x$CenType=="lt",ValuesToUse]<-x[x$CenType=="lt",ValuesToUse]*0.5
    x[x$CenType=="gt",ValuesToUse]<-x[x$CenType=="gt",ValuesToUse]*1.1

    #Tidy data and add on SeasonYear column
    x<-GetSeasonYear(x,ValuesToUse,Year)
    myFreq<-x$myFreq[1]; x$myFreq<-NULL
    x1<-x;x1$V1<-x1[,ValuesToUse]# Save raw data here for plotting later

    #if the HiCensor option is NOTfalse and there are censored values
    if(HiCensor!=FALSE & sum(x$CenType == "lt") > 0)   x<-ApplyHiCensor(x,ValuesToUse,HiCensor)

    x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])

    Median <- median(x[,  ValuesToUseforMedian], na.rm=T) #The median of the data

    #Take medians over each season
    Data<-MedianForSeason(x,ValuesToUse,Year)
    #..............................................................................#
    # do some tests to make data robust; must have at least 5 values, there must be >3 unique values.
    # CHECK TO SEE IF DATA OK TO CONTINUE
    AnalysisNote <- GetAnalysisNote(Data,ValuesToUse="V1",SecondTierTest = TRUE)

    if(AnalysisNote != "ok") { # if any TRUE dont do the test and return NA values
      VarS<-NA; AnnualSenSlope<-NA; Intercept<-NA; Lci<-NA; Uci<-NA;  AnalysisNote<-GetAnalysisNote(Data,ValuesToUse="V1",RunLength=RunLength)
      Sen_Probability<-NA; Probabilitymin<-NA;Probabilitymax<-NA; Percent.annual.change<-NA
    }
    else {

      Data <- Data[order(Data$NewDate), ] # restores the time ordering

      TheSlopes <- GetInterObservationSlopes(Data) # calculated slopes
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
      Data$Censored <- FALSE #
      Data$CenType <- "not"

      # estimate confidence intervals for 1-2*Alpha (where alpha = 0.05)
      Nprime <- length(TheSlopes$Slopes)             # number of calculated slopes
      N <- sum(!is.na(Data$V1))            # number of data points (years by freq less missing values)
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
      # R0 <- floor(approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0)$y)
      R0 <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0)$y

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
      Probabilitymax <- pnorm(Z1minusAlphamax)
      Probabilitymin <- pnorm(Z1minusAlphamin)

      # get intercept.
      Time <- (1:length(Data$V1))-1
      Ymed <-  median(Data$V1, na.rm=T) # the median of the measurements that are used to compute slope.
      Tmed <- median(Time)/myFreq  # the median of the time
      Intercept <- Ymed - (AnnualSenSlope*Tmed)

      Percent.annual.change=AnnualSenSlope/abs(Median) * 100 #have to use abs, as sometimes median -ve after flow adjustment

      if(doPlot) {
        PlotTrend(Data,x1,Intercept=Intercept,AnnualSenSlope=AnnualSenSlope,Lci=Lci,Uci=Uci,mymain=mymain,ValuesToUse = "V1",
                  IsSeasonal=FALSE,Ymed=Ymed,Tmed=Tmed,Percent.annual.change=Percent.annual.change,...)
      }

    } # end else
  } # end of first ifelse statement

  return(data.frame(Median=Median, Sen_VarS=VarS, AnnualSenSlope=AnnualSenSlope, Intercept=Intercept,
                    Lci=Lci, Uci=Uci, AnalysisNote=AnalysisNote, Sen_Probability=Sen_Probability, Probabilitymax, Probabilitymin,
                    Percent.annual.change=Percent.annual.change))
} # endSenSlope function


#######################################################################
# Seasonal Kendall test.
#######################################################################
SeasonalKendall <- function(x, ValuesToUse = "RawValue", Year = "Year", HiCensor=F, ...) {
  # Calculates the seasonal Kendall test
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the season (name =  "Season") and the year (name = "Year" )

  AnalysisNote <- GetAnalysisNote(x,ValuesToUse)   # CHECK TO SEE IF DATA OK TO CONTINUE

  if(AnalysisNote != "ok") { # if any TRUE dont do the test and return NA values
    VarS<-NA; S<-NA; D<-NA; tau<-NA;  Z<-NA; p<-NA; n<-length(unique(x[, ValuesToUse])) ;Probability=NA
  } else {

    #Tidy data and add on SeasonYear column
    x<-GetSeasonYear(x,ValuesToUse,Year)

    #if the HiCensor option is not FALSE and there are censored values
    if(HiCensor!=FALSE & sum(x$CenType == "lt") > 0)   x<-ApplyHiCensor(x,ValuesToUse,HiCensor)
    x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])

    #Take medians over each season
    Data<-MedianForSeason(x,ValuesToUse,Year)
    # CHECK TO SEE IF DATA OK TO CONTINUE
    AnalysisNote <- GetAnalysisNote(Data,ValuesToUse = "V1",IsSeasonal = TRUE,SecondTierTest = TRUE)

    if(AnalysisNote != "ok") { # dont do the test if not sufficient variation, non censored values or long runs
      VarS<-NA; S<-NA;  D<-NA;  tau<-NA;  Z<-NA;  p<-NA;  n<-length(unique(x[, ValuesToUse]));Probability=NA
    } else {

      Data <- Data[order(Data$NewDate), ] # restores the time ordering
      # get the years and seasons
      Data$Season <- unpaste(Data$SeasonYear, sep="-")[[1]]
      Data$Year <- as.numeric(unpaste(Data$SeasonYear, sep="-")[[2]])

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

        Probability<-1-0.5*p
        Probability[which(S>0)]<-0.5*p[which(S>0)]
      } else {
        AnalysisNote="Insufficent data to complete Seasonal Mann Kendall"
        Z <- NA
        p <- NA
        Probability<-NA
      }
    } # end of second if-else
  } # end of first if-else
  prop.censored <- sum(x$Censored)/nrow(x) # proportion of values that were censored
  prop.unique <- length(unique(x[, ValuesToUse]))/nrow(x) # proportion of values that were censored
  no.censorlevels<-length(unique(x[x$Censored=="TRUE",ValuesToUse]))
  return(data.frame(nObs=n, S=S, VarS=VarS, D=D, tau=tau, Z=Z, p=p,Probability=Probability,
                    AnalysisNote=AnalysisNote,prop.censored=prop.censored,prop.unique=prop.unique,no.censorlevels=no.censorlevels))
}

#######################################################################
#Seasonal Sen slope.
#######################################################################
SeasonalSenSlope <- function(x, ValuesToUse = "RawValue",ValuesToUseforMedian="RawValue", Year = "Year",
                             HiCensor=F, doPlot=FALSE, mymain=NULL,...) {
  # Calculates the seasonal sen slope and its 90% CIs - for 95% confidence in trend direction.
  # input data of dates and observations must have RawValue and Censored columns produced by RemoveAlphaDetect
  # input data must have observation dates in myDate being a date vector of class "Date"
  # input data must have columns specifying the season (name "Season") and the year  (names "Season"Year")

  AnalysisNote <- GetAnalysisNote(x, ValuesToUse = ValuesToUse,...)   # CHECK TO SEE IF DATA OK TO CONTINUE
  if(AnalysisNote != "ok") { # if not "ok" dont do the test
    Median<-NA; VarS<-NA;   AnnualSenSlope<-NA;  Intercept<-NA;  Lci<-NA;  Uci<-NA;
    Sen_Probability<-NA;  Probabilitymax <- NA;   Probabilitymin <- NA;  Percent.annual.change<-NA
  } else {

    # For calculating the sen slope, replace "lt" censored values with half the censored values
    x[x$CenType=="lt",ValuesToUse]<-x[x$CenType=="lt",ValuesToUse]*0.5
    x[x$CenType=="gt",ValuesToUse]<-x[x$CenType=="gt",ValuesToUse]*1.1

    #Tidy data and add on SeasonYear column
    x<-GetSeasonYear(x,ValuesToUse,Year)
    myFreq<-x$myFreq[1]; x$myFreq<-NULL
    x1<-x ; x1$V1<-x1[,ValuesToUse]#Saved here to use for plotting later

    #if the HiCensor option is true and there are censored values
    if(HiCensor!=FALSE & sum(x$CenType == "lt") > 0)   x<-ApplyHiCensor(x,ValuesToUse,HiCensor)
    x[,ValuesToUse]<-as.numeric(x[,ValuesToUse])

    Median <- median(x[,  ValuesToUseforMedian], na.rm=T) #The median of the data

    #Take medians over each season
    Data<-MedianForSeason(x,ValuesToUse,Year)
    Data$Season<-factor(as.character(Data[, "Season"]),levels=unique(Data[,"Season"])) #Adjust the factors for the seasons for the analysis to allow completely missing seasons
    # CHECK TO SEE IF DATA OK TO CONTINUE
    AnalysisNote <- GetAnalysisNote(Data,ValuesToUse = "V1",IsSeasonal = TRUE,SecondTierTest = TRUE)

    # dont do the test if not sufficient variation, non censored values or long runs
    if(AnalysisNote != "ok") {
      VarS<-NA;  AnnualSenSlope<-NA;  Intercept<-NA;  Lci<-NA;  Uci<-NA;
      Sen_Probability<-NA;  Probabilitymax <- NA;   Probabilitymin <- NA;  Percent.annual.change<-NA
    } else {

      Data <- Data[order(Data$NewDate), ] # restores the time ordering
      # get the years and seasons
      Data$Season <- unpaste(Data$SeasonYear, sep="-")[[1]]
      Data$Year <- as.numeric(unpaste(Data$SeasonYear, sep="-")[[2]])
      Data$myDate<-Data$NewDate

      TheSlopes <- ddply(Data, c("Season"), function(y) GetInterObservationSlopes(y))
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
      Data$Censored <- FALSE # Note the calculation of VarS and Senslope does not include censored values, but these fields are required by GetKendall
      Data$CenType <- "not"

      # estimate confidence intervals for 1-2*Alpha (where alpha = 0.05)
      Nprime <- length(TheSlopes$Slopes)             # number of calculated slopes
      N <- sum(!is.na(Data$V1))            # number of data points (years by freq less missing values)
      VarS <- SeasonalKendall(x=Data, ValuesToUse="V1")$VarS # get the variance using the Seasonal Kendall test function (ignore censored values)
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
      R0 <- approx(y=RankOfSlopes, x=sort(TheSlopes$Slopes), xout = 0)$y
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
      Probabilitymax <- pnorm(Z1minusAlphamax)
      Probabilitymin <- pnorm(Z1minusAlphamin)

      # get intercept.
      Time<-(1:length(Data$V1))-1 # need all dates despite some being censored.
      Ymed <-  median(Data$V1, na.rm=T) # the median of the measurements that are used to compute slope.
      Tmed <- median(Time)/myFreq  # the median of the time
      Intercept <- Ymed - (AnnualSenSlope*Tmed)

      Percent.annual.change = AnnualSenSlope/abs(Median)*100 #using abs of median as sometimes median is -ve after flow adjustment

      if(doPlot) { # graphics.off(); x11()
        PlotTrend(Data,x1,Intercept=Intercept,AnnualSenSlope=AnnualSenSlope,Lci=Lci,Uci=Uci,mymain=mymain,ValuesToUse = "V1",
                  IsSeasonal=TRUE,Ymed=Ymed,Tmed=Tmed,Percent.annual.change=Percent.annual.change,...)
      }
    } # end of second if-else statement
  }  # end of first if-else statement


  return(data.frame(Median=Median, Sen_VarS=VarS, AnnualSenSlope=AnnualSenSlope, Intercept=Intercept,
                    Lci=Lci, Uci=Uci, AnalysisNote=AnalysisNote,
                    Sen_Probability=Sen_Probability, Probabilitymax = Probabilitymax,  Probabilitymin = Probabilitymin,
                    Percent.annual.change=Percent.annual.change))
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

  xx <- x$V1
  cx <- x$Censored
  yy <- seq(1:length(x$V1)) # y is time and only needs to be  a sequentially increasing value
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

  return(list(nObs = n, S = kenS, vars=varS,  D = J, tau = kenS/J, Z=Z, p=p))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#######################################################################
#                   TREND AGGREGATION FUNCTIONS
#######################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#######################################################################
#  Calculate aggregate proportion of improving trends (PIT Statistic)
#     Analytical Solution
#######################################################################
# this function assesses the PIT statistic for a single water quaility variable, analytically
AnalysePIT <- function(x,               # dataframe of trend assessement results that  contains columns named "npID" and  "Probability", which is the probability the true trend was decreasing
                       myRound = 1,
                       Reverse = c("CLAR","MCI")) { # variables for which decreasing = degrading

  x <- x[!is.na(x$Probability), ] # remove NA values that result from NotAnalysed
  N <- nrow(x)

  PIT <- round((sum(x$Probability >0.5) + sum(x$Probability == 0.5)/2)/N *100, myRound)
  if(x$npID[1] %in% Reverse){PIT<-100-PIT}  #If an exception, take the complement of the PIT to make it probabilitie of improving trend

  VarP <- 1/N^2 * sum(x$Probability * (1-x$Probability))
  sdPIT <- round(sqrt(VarP)*100, myRound)
  return(data.frame(n.Sites=N, PIT = PIT, sdPIT= sdPIT))
}
#######################################################################
#  Calculate aggregate proportion of improving trends (PIT Statistic)
#     Solved by Monte Carlo simulation
#######################################################################
# this function assesses the PIT statistic for a single water quality variable using a Monte Carlo simulation
MonteCarloPIT <- function(x,   # dataframe of trend assessement results that contains columns named "npID" and  "Probability"
                          myRound = 1,  # for rounding the results
                          Realise = 1000,#  number of realisations to run
                          Reverse = c("CLAR","MCI")) { # variables for which decreasing = degrading

  x<-x[!is.na(x$Probability),]   #Remove NA values
  N <- nrow(x)
  # the MC simulation
  SIM <- lapply(x[, "Probability"],function(y) rbinom(Realise,1,as.numeric(y))) # Bernoulli (success/fail) outcome at site
  SIM <- ldply(SIM, function(x) return(x))
  mySimulationResults<-data.frame(PIT = round(100*mean(colSums(SIM)/N), myRound), sdPIT = round(100*sd(colSums(SIM)/N),myRound)) # count successes and divide by total sites
  mySimulationResults$n.Sites <- nrow(x)
  mySimulationResults <- mySimulationResults[, c( "n.Sites", "PIT", "sdPIT")]
  if(x$npID[1] %in% Reverse){mySimulationResults$PIT<-100-mySimulationResults$PIT}
  return(mySimulationResults)
}

#######################################################################
#  Calculate counts of increasin and decreasing trends
#######################################################################
#This provides a face value count of in INCREASING or DERCREASING trends,.  Note this is not improvement
FaceValueCounter <- function(x,  myRound = 1,Reverse=Reverse) {  # for rounding the results) { # x = dataframe of trend assessement results that  contains column named "AnnualSenSlope"

  x <- x[!is.na(x$AnnualSenSlope), ] # remove NA values that result from NotAnalysed
  n <- nrow(x)
  if(n>0) {
    n.Decreasing <- sum(x$AnnualSenSlope < 0) + round(sum(x$AnnualSenSlope == 0)/2)  # assign zero trends 50:50, myRound)
    n.Increasing <- sum(x$AnnualSenSlope > 0) + round(sum(x$AnnualSenSlope == 0)/2)  # assign zero trends 50:50, myRound)

    PropDecreasing <- round(n.Decreasing/n * 100, myRound)
    PropIncreasing <- round(n.Increasing/n * 100, myRound)

  } else {
    n.Decreasing <- NA
    n.Increasing <-NA
    PropDecreasing <- NA
    PropIncreasing <- NA
  }
  if(x$npID[1] %in% Reverse){
    PropImproving<-PropIncreasing
    PropDegrading<-PropDecreasing
    n.Improving<-n.Increasing
    n.Degrading<-n.Decreasing
  }else{
    PropDegrading<-PropIncreasing
    PropImproving<-PropDecreasing
    n.Degrading<-n.Increasing
    n.Improving<-n.Decreasing
  }
  return(data.frame(n=n, n.Degrading = n.Degrading, n.Improving= n.Improving, PropDegrading=PropDegrading, PropImproving=PropImproving)) #
}

#######################################################################
#  Assign Categorical improvement confidence (IPCC defintions)
#######################################################################
#Function to provide caterogrical confidence categories for the probability of improvement for trend assessment results
ImprovementConfCat<-function(x,Reverse=c("CLAR","MCI")){ # x = dataframe of trend assessement results that  contains column named "Probability"
  # Reverse = vector of variable names where increasing trends indicate improvement

  P<-x$Probability

  if(!is.na(Reverse[1])){P[x$npID %in% Reverse]<-1-P[x$npID %in% Reverse]}

  ConfCats <-cut(P, breaks =c(-0.01, 0.01, 0.05, 0.1, 0.33,0.67, 0.9, 0.95, 0.99, 1.01),
                 labels = c("Exceptionally unlikely",    #These breaks and labels are based on guidance in IPCC report
                            "Extremely unlikely",
                            "Very unlikely",
                            "Unlikely",
                            "As likely as not",
                            "Likely",
                            "Very likely",
                            "Extremely likely",
                            "Virtually certain"))
  ConfCats <-as.character(ConfCats )
  ConfCats [is.na(ConfCats )]<-"Not Analysed"

  ConfCats <-factor(ConfCats ,  levels = c("Exceptionally unlikely",
                                           "Extremely unlikely",
                                           "Very unlikely",
                                           "Unlikely",
                                           "As likely as not",
                                           "Likely",
                                           "Very likely",
                                           "Extremely likely",
                                           "Virtually certain",
                                           "Not Analysed"))
  return (ConfCats)
}


#######################################################################
#  Assign Categorical improvement confidence (LAWA defintions)
#######################################################################
#Function to provide caterogrical confidence categories for the probability of improvement for trend assessment results
ImprovementConfCatLAWA<-function(x,Reverse=c("CLAR","MCI")){ # x = dataframe of trend assessement results that  contains column named "Probability"
  # Reverse = vector of variable names where increasing trends indicate improvement

  P<-x$Probability

  if(!is.na(Reverse[1])){P[x$npID %in% Reverse]<-1-P[x$npID %in% Reverse]}

  ConfCats <-cut(P, breaks =c(0, 0.1, 0.33,0.67, 0.9, 1.0),
                 labels = c("-2",    #These breaks and labels are based on outputs used by LAWA
                            "-1",
                            "0",
                            "1",
                            "2"))

  ConfCats <-as.numeric(as.character(ConfCats ))
  ConfCats [is.na(ConfCats )]<-"Not Analysed"

  return (ConfCats)
}
