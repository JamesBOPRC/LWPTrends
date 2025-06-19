####################################################################
#### method for searching for a start point for trend analysis. ####
####################################################################

ContinuousTS <- function (data, gap)
{

  require(dplyr)

  names(data) <- c("SiteID", "Time", "Value")
  Min_date <- as.Date(Sys.time(),tz="etc/GMT+12")
  Max_date <- as.Date(max(data$Time),tz="etc/GMT+12")


  dat_sub <- as.Date(data$Time,tz="etc/GMT+12")
  dat_sub <- dat_sub[order(dat_sub, decreasing = T)]
  dat_diff <- as.numeric(diff(dat_sub))
  if (Max_date - max(dat_sub) < (5 * 365)) {
    if (min(dat_diff) < (-gap)) {
      first_exceed <- which(dat_diff < (-gap))
      first_exceed_date <- dat_sub[min(first_exceed)]
      Min_date <- if_else(first_exceed_date < Min_date,
                          first_exceed_date, Min_date)
    }
    else {
      oldest_date <- min(dat_sub)
      Min_date <- if_else(oldest_date < Min_date,
                          oldest_date, Min_date)
    }
  }

  return(Min_date)
}

#####################################################################
##### Basic Trend Analysis.  Allows you to fix flow adjustment, #####
##### and seasonal adjustment.  There is no 'automation'.       #####
##### adjustment is either on or off for each covariate         #####
##### This was previously known as 'Trend_Calc'                 #####
#####################################################################

#Data should be in the format - "Site","LocationName","Time","FLOW", PARAMETERS 1:X
# you need to rename the parameters in capitals....
# myDate must be present in the timeframe

Fixed_Trend_Analysis <- function(Data,
                                 param_colnames,
                                 siteID_col,
                                 time_col = "Time",
                                 Timeframe=NA,
                                 Seasonal = F,
                                 Co_variate_Adjust=F,
                                 Co_variate_col=NA,
                                 Data_gap = 3,
                                 log_parameters = c("ECOLI","TSS","FC","ENT")){


  library(plyr)
  library(lubridate)
  library(NADA)
  library(gam)
  library(rlang)

  #call data
  Data <- Data

  #define parameter column names
  param_colnames <- param_colnames

  #set the siteID as a factor
  #Data <-   Data %>% mutate(!!sym(siteID_col) := as.factor(!!sym(siteID_col)))
  Sites <-  Data %>% pull(siteID_col) %>% unique()

  #create an empty final DF
  Final_Df <- NULL

  #for each site
  for(i in Sites){

    SiteName=i

    #for each defined parameter (as column names)
    for(p in param_colnames){

      #if flow adjustment is required - create dataframe with siteID, time and flow
      if(Co_variate_Adjust==T){
        #make a dataset that

        if(is.na(Co_variate_col)){stop("Covariate adjustment is true but you haven't defined a co-variate column.")}



        current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p,Co_variate_col)]
      }else{


        current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p)]
      } #end loop

      #ensure no missing data for trend analysis
      current_dat <- current_dat[complete.cases(current_dat),]

      #this is to ensure that there is actually data
      if(nrow(current_dat)==0){
        cat(paste0("There is no ",p," data at ",i,"/n"))

        #this is the other case if no data.
        Trend_Analysis <- data.frame( "nObs" = "NOT ENOUGH DATA","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                      "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                      "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                      "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA, "Season_ADJ"= NA,"Season_ADJ_p"= NA ,"Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"= NA,
                                      "Site"= i,"analyte"= p,"Start_Year"= Start_Year,"End_Year"= End_Year,"Time_Period"=Timeframe)


        #this loop assumes a suitable dataframe
      } else if(nrow(current_dat)>0){


        #no definition so give all data
        if(is.na(Timeframe)){
          current_dat <-current_dat %>% filter(as.Date(current_dat[[time_col]],tz="etc/GMT+12")>as.Date(ContinuousTS(current_dat, gap = Data_gap * 365), tz = "etc/GMT+12"))
          Start_Year <-  format(min(current_dat[[time_col]],na.rm=T),"%Y")
          End_Year <-  format(max(current_dat[[time_col]],na.rm=T),"%Y")

        }else{
          #otherwise subset to the current date minus timeframe as a number
          current_dat <-current_dat %>% filter(as.Date(current_dat[[time_col]],tz="etc/GMT+12")>
                                                 (as.Date(max(current_dat[[time_col]]), tz = "etc/GMT+12")-years(Timeframe)))

          #and then check to see if the dataset is continuous
          current_dat <-current_dat %>% filter(as.Date(current_dat[[time_col]],tz="etc/GMT+12")>
                                                 as.Date(ContinuousTS(current_dat, gap = Data_gap * 365), tz = "etc/GMT+12"))

          End_Year <-  format(max(current_dat[[time_col]],na.rm=T),"%Y")
          Start_Year <-  format(min(current_dat[[time_col]],na.rm=T),"%Y")
        }

        #does the parameter need to be log scaled?
        if(p %in% log_parameters){
          current_dat[,3] <- log10(current_dat[,3]+1)
          cat(paste0("Parameter ",p," has been log10 transformed at ", i,"\n"))

          }

        #only run trend analysis if you have more than 10 values
        if(nrow(current_dat[,1:3])>=10){

          #create a new column called myDate
          current_dat <- current_dat %>% mutate(myDate = as.Date(current_dat[[time_col]],tz="etc/GMT+12"))
          current_dat <- GetMoreDateInfo(current_dat)

          #remove alpha detects
          current_dat <-  RemoveAlphaDetect(current_dat,ColToUse = p)

          #inspect data
          current_dat <- InspectTrendData(current_dat)



#if the input forces covariate adjustment
          if(Co_variate_Adjust==T){

            #run covariate adjustment models
            adj <- AdjustValues(current_dat[complete.cases(current_dat),], method = c("Gam","LOESS"),main=paste(i, ": Co_variate_Adjustment_",Timeframe,sep=""), ValuesToAdjust = 'RawValue', Covariate = Co_variate_col, Span = c(0.7,0.8,0.9),doPlot = F, plotpval=T)

            #work out the smallest p value
            minflowp <- min(adj$LOESS0.7_p[1],adj$Gam_p[1])

            #name of the CV adjustment model used
            CVmod <- colnames(adj[which(adj[1,]==minflowp)])

            #removing _p from string
            CVmod_data <- substr(CVmod,1,nchar(CVmod)-2)

            #merge current dat with adjusted values
            current_dat<-merge(current_dat[complete.cases(current_dat),],adj[,c('myDate',CVmod_data)])

            cat(paste0("Parameter ",p," will be adjusted by ",Co_variate_col,  " at ", i," with a pvalue = ",as.numeric(round(minflowp)),"\n"))


            seasonp <- GetSeason(current_dat)[[2]][1,3]
    #if the input forces seasonal adjustment
            if(Seasonal==T){
              #carry out a seasonal trend analysis
              cat(paste0("Seasonal adjustment for ",p," will be carried out at ", i," with a pvalue of ",as.numeric(round(seasonp)), "\n"))
              current_dat <- GetSeason(current_dat)[[1]]
              Trend_Analysis <- SeasonalTrendAnalysis(current_dat,ValuesToUse=CVmod_data,Timeframe,doPlot=F)
              Trend_Analysis$Season_ADJ <- "YES"
              Trend_Analysis$Season_ADJ_p <- GetSeason(current_dat)[[2]][1,3]

            #covariate adjustment but no seasonal adjustment
            }else{
              cat(paste0("No seasonal adjustment for ",p," will be carried out at ", i,"\n"))
              Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,ValuesToUse=CVmod_data,doPlot=F)
              Trend_Analysis$Season_ADJ <- "NO"
              Trend_Analysis$Season_ADJ_p <- "N/A"
            }
    #metadata to state that it was covariate adjusted
            Trend_Analysis$CV_ADJ <- "YES"
            Trend_Analysis$CV_MOD <- CVmod_data
            Trend_Analysis$CV_MOD_p <- minflowp

#if the input does not force covariate adjustment
          }else{
            cat(paste0("No covariate adjustment for ",p," will be carried out at ",i,"\n"))

            seasonp <- GetSeason(current_dat)[[2]][1,3]
            if(Seasonal==T){
              #carry out a seasonal trend analysis

              current_dat <- GetSeason(current_dat)[[1]]
              cat(paste0("Seasonal adjustment for ",p," will be carried out at ", i," with a pvalue of ",as.numeric(round(seasonp)), "\n"))
              Trend_Analysis <- SeasonalTrendAnalysis(current_dat,doPlot=F)#this doesn't need a 'valuesToUse' as it will default to RawValue
              Trend_Analysis$Season_ADJ <- "YES"
              Trend_Analysis$Season_ADJ_p <- GetSeason(current_dat)[[2]][1,3]

            }else{
              cat(paste0("No seasonal adjustment for ",p," will be carried out at ", i,"\n"))
              Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,doPlot=T)
              Trend_Analysis$Season_ADJ <- "NO"
              Trend_Analysis$Season_ADJ_p <- "N/A"
            }

#metadata to state that it was not covariate adjusted
            Trend_Analysis$CV_ADJ <- "NO"
            Trend_Analysis$CV_MOD <- "N/A"
            Trend_Analysis$CV_MOD_p <- "N/A"


          }

          #additional metadata
          Trend_Analysis$Site <- i
          Trend_Analysis$analyte <- p
          Trend_Analysis$Start_Year <-Start_Year
          Trend_Analysis$End_Year <- End_Year
          Trend_Analysis$Time_Period <- Timeframe

        }else{
          #this is the other case if not more than 10 values - dataframe of NA values.
          Trend_Analysis <- data.frame( "nObs" = "NOT ENOUGH DATA","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                        "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                        "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                        "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA, "Season_ADJ"= NA,"Season_ADJ_p"= NA ,"Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"= NA,
                                        "Site"= i,"analyte"= p,"Start_Year"= Start_Year,"End_Year"= End_Year,"Time_Period"=Timeframe)

        }
        #merge together the holder and the trend analysis output
        Final_Df <- rbind(Final_Df,Trend_Analysis)

      }
    }#end of parameter loop
  }#end of site loop
  return(Final_Df)
  detach(plyr)
}#end of function



############################################################################################################

#############################################################################
##### Flexible Trend Analysis.                                          #####
##### The rules for this are flow adjust if FlowMod p-value is < 0.05   #####
#####       and seasonally adjust if ADJ_p <0.05.                       #####
#############################################################################

#season will be automatically detected.
Automated_Trend_Analysis <- function(Data,
                                     param_colnames,
                                     siteID_col,
                                     time_col = "Time",
                                     Timeframe="All", #or a number that will be converted to years
                                     Co_variate_col=NA,
                                     Data_gap = 3,
                                     min_values = 20,
                                     log_parameters = c("ECOLI","TSS","FC","ENT")){


  library(plyr)
  library(lubridate)
  library(NADA)
  library(gam)
  library(rlang)


  #call data
  Data <- Data

  #define parameter column names
  param_colnames <- param_colnames

  #set the siteID as a factor
  #Data <-   Data %>% mutate(!!sym(siteID_col) := as.factor(!!sym(siteID_col)))
  Sites <-  Data %>% pull(siteID_col) %>% unique()

  #create an empty final DF
  Final_Df <- NULL

  #for each site
  for(i in Sites){

    SiteName=i

    #for each defined parameter (as column names)
    for(p in param_colnames){

      #if there is a coariate column, put this in the dataset
      if(is.na(Co_variate_col)){
        current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p)] %>% filter(complete.cases(.))
        cat(paste0(p," data at ",i, " will not be adjusted as there is no defined co-variate\n"))


      }else{#if not then don't, and don't bother with flow adjustment
          current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p,Co_variate_col)]%>% filter(complete.cases(.))
          cat(paste0("Attempt to adjust ",p," by ",Co_variate_col,  " at ", i,"\n"))

      }


      #this is to ensure that there is actually data this will fail if there is an empty dataframe.
      if(nrow(current_dat)==0){
        cat(paste0("There is no ",p," data at ",i))

        #this is the other case if no data.
        Trend_Analysis <- data.frame( "nObs" = "NOT ENOUGH DATA","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                      "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                      "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                      "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA, "Season_ADJ"= NA,"Season_ADJ_p"= NA ,"Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"= NA,
                                      "Site"= i,"analyte"= p,"Start_Year"= Start_Year,"End_Year"= End_Year,"Time_Period"=Timeframe)


      } else if(nrow(current_dat)>0){

        #if timeframe was set to 'all' - SET current data frame to the correct timescale
        if(is.na(Timeframe)){
          current_dat <-current_dat %>% filter(as.Date(current_dat[[time_col]],tz="etc/GMT+12")>
                                                 as.Date(ContinuousTS(current_dat, gap = Data_gap * 365), tz = "etc/GMT+12"))
          Start_Year <-  format(min(current_dat[[time_col]],na.rm=T),"%Y")
          End_Year <-  format(max(current_dat[[time_col]],na.rm=T),"%Y")

        }else{
          #otherwise subset to the current date minus timeframe as a number
          current_dat <-current_dat %>% filter(as.Date(current_dat[[time_col]],tz="etc/GMT+12")>
                                                 (as.Date(max(current_dat[[time_col]]), tz = "etc/GMT+12")-years(Timeframe)))

          #and then check to see if the dataset is continuous
          current_dat <-current_dat %>% filter(as.Date(current_dat[[time_col]],tz="etc/GMT+12")>
                                                 as.Date(ContinuousTS(current_dat, gap = Data_gap * 365), tz = "etc/GMT+12"))

          End_Year <-  format(max(current_dat[[time_col]],na.rm=T),"%Y")
          Start_Year <-  format(min(current_dat[[time_col]],na.rm=T),"%Y")
        }
        #does the parameter need to be log scaled?
        if(p %in% log_parameters){
          current_dat[,3] <- log10(current_dat[,3]+1)
          cat(paste0("Parameter ",p," has been log10 transformed at ", i,"\n"))

          }

        #only run trend analysis if you have more than min_values
        if(nrow(current_dat[,1:3])>=min_values){

          #create a new column called myDate
          current_dat <- current_dat %>% mutate(myDate = as.Date(current_dat[[time_col]],tz="etc/GMT+12"))
          current_dat <- GetMoreDateInfo(current_dat)

          #remove alpha detects
          current_dat <-  RemoveAlphaDetect(current_dat,ColToUse = p)

          #inspect data
          current_dat <- InspectTrendData(current_dat)

          if(!is.na(Co_variate_col)){

            #### work out if the data should be Co_variate adjusted ####
            #run covariate adjustment models
            adj <- AdjustValues(current_dat[complete.cases(current_dat),], method = c("Gam","LOESS"),main=paste(i, ": Co_variate_Adjustment_",Timeframe,sep=""), ValuesToAdjust = 'RawValue', Covariate = Co_variate_col, Span = c(0.7,0.8,0.9),doPlot = F, plotpval=T)

            #work out the smallest p value
            minflowp <- min(adj$LOESS0.7_p[1],adj$Gam_p[1]) #based on the lowess 0.7 model or gam...

            ######## NEED TO INSERT SOME OUTPUT FOR CREATING PLOTS THAT SHOW WHAT THE ADJUSTMENT LOOKS LIKE ####

            if(minflowp <= 0.05){

              cat(paste0("A relationship between ",p," and ",Co_variate_col,  " detected at ", i," with a pvalue = ",as.numeric(round(minflowp)),"\n"))

              #name of the CV adjustment model used
              CVmod <- colnames(adj[which(adj[1,]==minflowp)])

              #removing _p from string
              CVmod_data <- substr(CVmod,1,nchar(CVmod)-2)

              #merge current dat with adjusted values
              current_dat<-merge(current_dat[complete.cases(current_dat),],adj[,c('myDate',CVmod_data)])

              Co_variate_Adjust = T
            }else{
              cat(paste0("No relationship between ",p," and ",Co_variate_col,  " detected at ", i,". Pvalue = ",as.numeric(round(minflowp)),"\n"))
              Co_variate_Adjust = F
            }
          }else{
            Co_variate_Adjust = F
          }

        #### work out if the data should be seasonally adjusted ####
            #add appropriate season to data
          current_dat <- GetSeason(current_dat)[[1]]
            #run seasonality test
            seasonp <- GetSeason(current_dat)[[2]][1,3]

            if(seasonp <= 0.05){
              cat(paste0("Seasonality detected in ",p," at ", i," with a pvalue of ",as.numeric(round(seasonp)), "\n"))

              ######## NEED TO INSERT SOME OUTPUT FOR CREATING PLOTS THAT SHOW WHAT THE SEASONAL DIFFERENCES LOOK LIKE ####
              if(Co_variate_Adjust == T){

                Trend_Analysis <- SeasonalTrendAnalysis(current_dat,ValuesToUse=CVmod_data,Timeframe,doPlot=F)
                Trend_Analysis$CV_ADJ <- "YES"
                Trend_Analysis$CV_MOD <- CVmod_data
                Trend_Analysis$CV_MOD_p <- minflowp

              }else{# if not co-variate adjusted

                Trend_Analysis <- SeasonalTrendAnalysis(current_dat,Timeframe,doPlot=F)
                #metadata to state that it was not covariate adjusted
                Trend_Analysis$CV_ADJ <- "NO"
                Trend_Analysis$CV_MOD <- "N/A"
                Trend_Analysis$CV_MOD_p <- "N/A"

              }

              Trend_Analysis$Season_ADJ <- "YES"
              Trend_Analysis$Season_ADJ_p <- seasonp

            }else{ #if there is no seasonality in data
              cat(paste0("No seasonality detected in ",p," at ", i,". Pvalue = ",as.numeric(round(seasonp,4)),"\n"))

              if(Co_variate_Adjust == T){

                Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,ValuesToUse=CVmod_data,doPlot=F)
                Trend_Analysis$CV_ADJ <- "YES"
                Trend_Analysis$CV_MOD <- CVmod_data
                Trend_Analysis$CV_MOD_p <- minflowp

              }else{# if not co-variate adjusted

                Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,doPlot=F)
                #metadata to state that it was not covariate adjusted
                Trend_Analysis$CV_ADJ <- "NO"
                Trend_Analysis$CV_MOD <- "N/A"
                Trend_Analysis$CV_MOD_p <- "N/A"

              }

              Trend_Analysis$Season_ADJ <- "NO"
              Trend_Analysis$Season_ADJ_p <- "N/A"

            }

            #additional metadata
            Trend_Analysis$Site <- i
            Trend_Analysis$analyte <- p
            Trend_Analysis$Start_Year <-Start_Year
            Trend_Analysis$End_Year <- End_Year
            Trend_Analysis$Time_Period <- Timeframe


            }else{
              #this is the other case if not more than 10 values - dataframe of NA values.
              Trend_Analysis <- data.frame( "nObs" = "NOT ENOUGH DATA","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                            "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                            "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                            "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA, "Season_ADJ"= NA,"Season_ADJ_p"= NA ,"Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"= NA,
                                            "Site"= i,"analyte"= p,"Start_Year"= Start_Year,"End_Year"= End_Year,"Time_Period"=Timeframe)

            }



      }

      #merge together the holder and the trend analysis output
      Final_Df <- rbind(Final_Df,Trend_Analysis)

      }#end of parameter loop
    }#end of site loop
  return(Final_Df)
  detach(plyr)
  }#end of function

########

Outliers_Z_Score <- function(dataset,probability = "0.99",output = "values",
                             type = "z"){

  #c("z", "t", "chisq", "iqr", "mad")

  names(dataset) <- c("LocationName", "Site", "Time",  "Parameter", "Value")

  library(tidyr)
  Final_Dataset <- NA
  Final_Outliers <- NA


  #assume LocationName is consistent
  for (i in unique(dataset$LocationName)){

    #these are locked for this project.  As columns.
    for (x in unique(dataset$Parameter)){

      Data_to_check <- dataset[dataset$LocationName == i & dataset$Parameter == x, ]

      Data_to_check <- Data_to_check[complete.cases(Data_to_check),]
      Outliers <- Data_to_check[scores(Data_to_check[,"Value"], type=type, prob=probability),]

      Output_Data <- anti_join(Data_to_check,Outliers)

      if(nrow(Output_Data >0)){
        Output_Data$LocationName = i
        Output_Data$Parameter = x
        names(Output_Data) <- c("LocationName","Site","Time","Parameter","Value")
        Final_Dataset <- rbind(Final_Dataset, Output_Data)
      }


      if(nrow(Outliers >0)){
        Outliers$LocationName = i
        Outliers$Parameter = x
        names(Outliers) <- c("LocationName","Site","Time","Parameter","Value")
        Final_Outliers <- rbind(Final_Outliers, Outliers)
      }

    }


  }

  Final_Dataset <- Final_Dataset[-1,]
  Final_Outliers <- Final_Outliers[-1,]

  if(output == "values"){
    return(Final_Dataset)
  } else if (output == "outliers"){
    return(Final_Outliers)
  }


}



