
#method for searching for a start point for trend analysis.
ContinuousTS <- function (data, gap)
{

  require(dplyr)

  names(data) <- c("SiteID", "Time", "Value")
  Min_date <- as.Date(Sys.time())
  Max_date <- as.Date(max(data$Time))


  dat_sub <- as.Date(data[, 2],tz="+12")
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


################################################################################################################
##### Basic Trend Analysis.  Allows you to fix flow adjustment, and seasonal adjustment.  No  'automation'.
#Trend_Calc

Fixed_Trend_Analysis <- function(Data,param_colnames,siteID_col,
                       time_col = "Time",
                       Timeframe="All",
                       Seasonal = F,
                       Co_variate_Adjust=F,
                       Co_variate_col=NA,
                       log_parameters = c("ECOLI","TSS","FC","ENT")){


  require(plyr)
  require(lubridate)
  require(NADA)
  require(gam)

  Data <- Data
  param_colnames <- param_colnames
  Data[,siteID_col] <- as.factor(Data[,siteID_col])
  Final_Df <- NULL

  #for each site
  for(i in levels(Data[,siteID_col])){

    SiteName=i

    #for each defined parameter (as column names)
    for(p in param_colnames){

      #if flow adjustment is required
      if(Co_variate_Adjust==T){
        #make a dataset that
        current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p,flow_col)]
      }else{
        current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p)]
      }

      current_dat <- current_dat[complete.cases(current_dat),]

      if(nrow(current_dat)>0){


        if(Timeframe == "All"|Timeframe=="all"){
          current_dat <-current_dat[ as.Date(current_dat[,time_col]) > as.Date(ContinuousTS_2(current_dat,gap = 3*365)),]
          Start_Year <-  format(ContinuousTS_2(current_dat,gap = 3*365),"%Y")
          End_Year <-  format(max(current_dat[,time_col]),"%Y")
        }else{
          current_dat <-current_dat[as.Date(current_dat[,time_col]) > as.Date(max(current_dat[,time_col]))-(365*as.numeric(Timeframe)),]

          End_Year <-  format(max(current_dat[,time_col]),"%Y")
          Start_Year <-  format(max(current_dat[,time_col])- years(as.numeric(Timeframe)),"%Y")

        }


        if(p %in% log_parameters){current_dat[,3] <- log10(current_dat[,3]+1)}

        #only run trend analysis if more than 10 values
        if(nrow(current_dat[,1:3])>=10){

          current_dat$myDate <- as.Date(as.character(current_dat[,time_col]),"%Y-%m-%d")
          current_dat <- GetMoreDateInfo(current_dat)



          current_dat$Season<-current_dat$Month
          SeasonString<-levels(current_dat$Season)

          NewValues <- RemoveAlphaDetect(current_dat[,p])
          current_dat <- cbind.data.frame(current_dat, NewValues)



          #png(filename="test.png")
          #p1 <- InspectData(current_dat, StartYear = Start_Year, EndYear = End_Year, FlowtoUse = "FLOW",plotType = "TimeSeries", main= "Ex 1a Time Series of Monthly Data")
          #dev.off()

          #InspectData(current_dat, StartYear = Start_Year, EndYear = End_Year, plotType = "Matrix",PlotMat="RawData",main="Matrix of Values: monthly data")
          #InspectData(current_dat, StartYear = Start_Year, EndYear = End_Year, plotType = "Matrix",PlotMat="Censored",main="Matrix of censoring: monthly data")

          pst <- SeasonalityTest(current_dat,main=paste(i, ": Raw Monthly Data",sep=""),do.plot = F)

          if(Co_variate_Adjust==T){
            adj <- AdjustValues(current_dat[complete.cases(current_dat),], method = c("Gam","LOESS"),main=paste(i, ": Co_variate_Adjustment_",Timeframe,sep=""), ValuesToAdjust = 'RawValue', Covariate = Co_variate_col, Span = c(0.7,0.8,0.9),doPlot = F, plotpval=T)
            minflowp <- min(adj$LOESS0.7_p[1],adj$Gam_p[1])


            #          if(minflowp<0.05){

            flowmod <- colnames(adj[which(adj[1,]==minflowp)])
            flowmod_data <- substr(flowmod,1,nchar(flowmod)-2)
            current_dat<-merge(current_dat[complete.cases(current_dat),],adj[,c('myDate',flowmod_data)])

            png(filename=paste(i,"_",p,"_Co_variate_adjust_",Timeframe,".png",sep = ""))
            AdjustValues(current_dat[complete.cases(current_dat),], method = c("Gam","LOESS"),main=paste(i,"_",p,": Co_variate_Adjustment_",Timeframe,sep=""), ValuesToAdjust = 'RawValue', Covariate = Co_variate_col, Span = c(0.7,0.8,0.9),doPlot = T, plotpval=T)
            dev.off()

            if(Seasonal==T){
              png(filename=paste(i,"_",p,"_Seasonal_Trend_Plot_",Timeframe,".png",sep = ""))
              Trend_Analysis <- SeasonalTrendAnalysis(current_dat,ValuesToUse=flowmod_data,mymain=paste(i,"_",p,"_Co_variate_Adjusted: Seasonal Trend Analysis_",Timeframe,sep=""),doPlot=T)
              dev.off()

              Trend_Analysis$Season_ADJ <- "YES"
              Trend_Analysis$Season_ADJ_p <- pst$pvalue

              # #non FA TA as well
              # Trend_Analysis_2 <- SeasonalTrendAnalysis(current_dat, mymain=paste(i,"_",p,"_Raw Data: Seasonal Trend Analysis",sep=""),doPlot=T)
              # Trend_Analysis_2$Season_ADJ <- "YES"

            }else{
              png(filename=paste(i,"_",p,"_Non_Seasonal_Trend_Plot_",Timeframe,".png",sep = ""))
              Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,ValuesToUse=flowmod_data,mymain=paste(i,"_",p,"_Co_variate_Adjusted: Non-Seasonal Trend Analysis_",Timeframe,sep=""),doPlot=T)
              dev.off()
              Trend_Analysis$Season_ADJ <- "NO"
              Trend_Analysis$Season_ADJ_p <- "N/A"

            }

            Trend_Analysis$Flow_ADJ <- "YES"
            Trend_Analysis$Flow_MOD <- flowmod_data
            Trend_Analysis$Flow_MOD_p <- minflowp


          }else{#if not flow adjust = T

            if(Seasonal==T){
              png(filename=paste(i,"_",p,"_Trend_Plot_",Timeframe,".png",sep = ""))
              Trend_Analysis <- SeasonalTrendAnalysis(current_dat, mymain=paste(i,"_",p,"_Raw Data: Seasonal Trend Analysis_",Timeframe,sep=""),doPlot=T)
              dev.off()
              Trend_Analysis$Season_ADJ <- "YES"
              Trend_Analysis$Season_ADJ_p <- pst$pvalue

            }else{
              png(filename=paste(i,"_",p,"_Trend_Plot_",Timeframe,".png",sep = ""))
              Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,mymain=paste(i,"_",p,"_Raw Data: Non-Seasonal Trend Analysis_",Timeframe,sep=""),doPlot=T)
              dev.off()

              Trend_Analysis$Season_ADJ <- "NO"
              Trend_Analysis$Season_ADJ_p <- "N/A"
            }


            Trend_Analysis$Flow_ADJ <- "NO"
            Trend_Analysis$Flow_MOD <- "N/A"
            Trend_Analysis$Flow_MOD_p <- "N/A"



          }


          Trend_Analysis$Site <- i
          Trend_Analysis$Param <- p
          Trend_Analysis$Start_Year <-Start_Year
          Trend_Analysis$End_Year <- End_Year
          Trend_Analysis$Time_Period <- Timeframe

        }else{
          Trend_Analysis <- data.frame( "nObs" = "NOT ENOUGH DATA","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                        "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                        "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                        "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA, "Season_ADJ"= NA,"Season_ADJ_p"= NA ,"Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"= NA,
                                        "Site"= NA,"Param"= NA,"Start_Year"= NA,"End_Year"= NA,"Time_Period"=NA)


          Trend_Analysis$Site <- i
          Trend_Analysis$Param <- p
          Trend_Analysis$Start_Year <-Start_Year
          Trend_Analysis$End_Year <- End_Year
          Trend_Analysis$Time_Period <- Timeframe



        }



        Final_Df <- rbind(Final_Df,Trend_Analysis)


      }
    }
  }
  return(Final_Df)
}



############################################################################################################


#### The rules for this are flow adjust if FlowMod p-value is < 0.05
#### and seasonally adjust if ADJ_p <0.05.

Flexible_Trend_Analysis <- function(Data,param_colnames,siteID_col,
                                     time_col = "Time",
                                     Timeframe="All",
                                      Co_variate_col=NA,
                                     season = "Quarter",
                                    log_parameters = c("ECOLI","TSS","FC","ENT")){


  require(plyr)
  require(lubridate)
  require(NADA)
  require(gam)

  Data <- Data
  param_colnames <- param_colnames
  Data[,siteID_col] <- as.factor(Data[,siteID_col])
  Final_Df <- NULL

  for(i in levels(Data[,siteID_col])){

    SiteName=i

    for(p in param_colnames){

      #two different data files

      current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p,flow_col)]


      if(sum(!is.na(current_dat[,p]))>20){


        if(Timeframe == "All"|Timeframe=="all"){
          current_dat <-current_dat[ as.Date(current_dat[,time_col]) > as.Date(ContinuousTS_2(current_dat,gap = 3*365)),]


          Start_Year <-  format(ContinuousTS_2(current_dat,gap = 3*365),"%Y")
          End_Year <-  format(max(current_dat[,time_col]),"%Y")

          Data_Start <- min(current_dat[,time_col])
          Data_End <- max(current_dat[,time_col])



        }else{
          current_dat <-current_dat[as.Date(current_dat[,time_col]) > as.Date(max(current_dat[,time_col]))-(365*as.numeric(Timeframe)),]

          End_Year <-  format(max(current_dat[,time_col]),"%Y")
          Start_Year <- format(max(current_dat[,time_col])- years(as.numeric(Timeframe)),"%Y")
          Data_Start <- min(current_dat[,time_col])
          Data_End <- max(current_dat[,time_col])

        }

        if(p %in% log_parameters){current_dat[,3] <- log10(current_dat[,3]+1)}

        current_dat$myDate <- as.Date(as.character(current_dat[,time_col]),"%Y-%m-%d")
        current_dat <- GetMoreDateInfo(current_dat)


        if(season == "month"|season=="MONTH"|season=="Month"){

          current_dat$Season<-current_dat$Month
          SeasonString<-levels(current_dat$Season)
        }else if(season == "quarter"|season=="QUARTER"|season=="Quarter"){
          current_dat$Season<-current_dat$Qtr
          SeasonString<-levels(current_dat$Season)
        }else{
          current_dat$Season<-current_dat$CustomYear
          SeasonString<-levels(current_dat$Season)
        }


        NewValues <- RemoveAlphaDetect(current_dat[,p])
        current_dat <- cbind.data.frame(current_dat, NewValues)


        current_dat_FA <- current_dat[complete.cases(current_dat),]
        adj <- AdjustValues(current_dat_FA[complete.cases(current_dat_FA),], method = c("Gam","LOESS"),main=paste(i, ": Flow Adjustment",sep=""), ValuesToAdjust = 'RawValue', Covariate = Co_variate_col, Span = c(0.7,0.8,0.9),doPlot = F, plotpval=T)
        minflowp <- min(adj$LOESS0.7_p[1],adj$Gam_p[1])

        #which dataset to use?
        if(is.na(minflowp)==F & minflowp<0.05){

          flowmod <- colnames(adj[which(adj[1,]==minflowp)])
          flowmod_data <- substr(flowmod,1,nchar(flowmod)-2)
          current_dat_FA<-merge(current_dat_FA[complete.cases(current_dat_FA),],adj[,c('myDate',flowmod_data)])

          #create plot
          #          png(filename=paste(i,"_",p,"_Flow_adjust",".png",sep = ""))
          #          AdjustValues(current_dat[complete.cases(current_dat),], method = c("Gam","LOESS"),main=paste(i,"_",p,": Flow Adjustment",sep=""), ValuesToAdjust = 'RawValue', Covariate = Co_variate_col, Span = c(0.7,0.8,0.9),doPlot = T, plotpval=T)
          #          dev.off()

          Co_variate_Adjust = T

        }else{

          Co_variate_Adjust = F
        }

        #make sure have both datasets
        current_dat <- current_dat[complete.cases(current_dat[,-4]),-4]

        #Check for enough data to run analysis (20 values)
        if(Co_variate_Adjust ==T){
          Run_Analysis <- ifelse(nrow(current_dat_FA[,1:3])>=20,TRUE,FALSE)
        }else{
          Run_Analysis <- ifelse(nrow(current_dat[,1:3])>=20,TRUE,FALSE)
        }

        if(Run_Analysis == T){


          if(Co_variate_Adjust ==T){
            pst <- SeasonalityTest(current_dat_FA,main=paste(i, ": Co-variate Adjusted Monthly Data",sep=""),do.plot = F)
          }else{
            pst <- SeasonalityTest(current_dat,main=paste(i, ": Raw Monthly Data",sep=""),do.plot = F)
          }


          if(pst$pvalue<0.05){

            if(Co_variate_Adjust==T){
              Trend_Analysis <- SeasonalTrendAnalysis(current_dat_FA,ValuesToUse=flowmod_data,mymain=paste(i,"_",p,"_Co_variate_Adjusted: Seasonal Trend Analysis",sep=""),doPlot=F)

              if(Trend_Analysis$TrendCategory == "Not Analysed"){
                Trend_Analysis <- SeasonalTrendAnalysis(current_dat,ValuesToUse="RawValue",mymain=paste(i,"_",p,"_Co_variate_Adjusted: Seasonal Trend Analysis",sep=""),doPlot=F)

                Trend_Analysis$Flow_ADJ <- "Attempted"
                Trend_Analysis$Flow_MOD <- flowmod_data
                Trend_Analysis$Flow_MOD_p <- minflowp


              }else{

                Trend_Analysis$Flow_ADJ <- "YES"
                Trend_Analysis$Flow_MOD <- flowmod_data
                Trend_Analysis$Flow_MOD_p <- minflowp

              }



            }else{
              Trend_Analysis <- SeasonalTrendAnalysis(current_dat,ValuesToUse="RawValue",mymain=paste(i,"_",p,"_Raw Data: Seasonal Trend Analysis",sep=""),doPlot=F)
              Trend_Analysis$Flow_ADJ <- "NO"
              Trend_Analysis$Flow_MOD <- "N/A"
              Trend_Analysis$Flow_MOD_p <- "N/A"


            }

            Trend_Analysis$Season_ADJ <- "YES"
            Trend_Analysis$Season_ADJ_p <- pst$pvalue
            #Trend_Analysis <- Trend_Analysis[,c(1:8,10:18,9,19:29)]


          }else{

            if(Co_variate_Adjust==T){
              Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat_FA,ValuesToUse=flowmod_data,mymain=paste(i,"_",p,"_Co_variate_Adjusted: Non-Seasonal Trend Analysis",sep=""),doPlot=F)

              Trend_Analysis$Flow_ADJ <- "YES"
              Trend_Analysis$Flow_MOD <- flowmod_data
              Trend_Analysis$Flow_MOD_p <- minflowp

            }else{
              Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,ValuesToUse="RawValue",mymain=paste(i,"_",p,"_Raw Data: Non-Seasonal Trend Analysis",sep=""),doPlot=F)
              Trend_Analysis$Flow_ADJ <- "NO"
              Trend_Analysis$Flow_MOD <- "N/A"
              Trend_Analysis$Flow_MOD_p <- "N/A"
            }

            Trend_Analysis$Season_ADJ <- "NO"
            Trend_Analysis$Season_ADJ_p <- "N/A"

          }

        }

        #think this is incorrect

        Trend_Analysis$Site <- i
        Trend_Analysis$Param <- p
        Trend_Analysis$Start_Year <-Start_Year
        Trend_Analysis$End_Year <- End_Year
        Trend_Analysis$Time_Period <- Timeframe
        Trend_Analysis$Data_Start <- Data_Start
        Trend_Analysis$Data_End <- Data_End
        Trend_Analysis$Season <- season


      }else{
        Trend_Analysis <- data.frame( "nObs" = "INVALID DATASET","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                      "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                      "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                      "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA,
                                      "Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"=NA,"Season_ADJ"= NA, "Season_ADJ_p"= NA,"Site"= NA,"Param"= NA,"Start_Year"= NA,"End_Year"= NA,"Data_Start"=NA,"Data_End"=NA,"Season"=NA)

        Trend_Analysis$Site <- i
        Trend_Analysis$Param <- p
        Trend_Analysis$Start_Year <-Start_Year
        Trend_Analysis$End_Year <- End_Year
        Trend_Analysis$Time_Period <- Timeframe
        Trend_Analysis$Season <- season



      }






      Final_Df <- rbind(Final_Df,Trend_Analysis)


    }
  }

  return(Final_Df)

}


########################################################################################################

#Trend analysis with no adjustments for seasonality or co-variates.  This is just based on the raw data.
Basic_Trend_Analysis <- function(Data,param_colnames,siteID_col,
                                       time_col = "Time",
                                       Timeframe="All",
                                       season = "Quarter",
                                 log_parameters = c("ECOLI","TSS","FC","ENT")){
  require(plyr)
  require(lubridate)
  require(NADA)
  require(gam)

  Data <- Data
  param_colnames <- param_colnames
  Data[,siteID_col] <- as.factor(Data[,siteID_col])
  Final_Df <- NULL

  for(i in levels(Data[,siteID_col])){

    SiteName=i

    for(p in param_colnames){

      #two different data files

      current_dat <- Data[Data[[siteID_col]]==i,c(siteID_col,time_col,p)]


      if(sum(!is.na(current_dat[,p]))>20){


        if(Timeframe == "All"|Timeframe=="all"){
          current_dat <-current_dat[ as.Date(current_dat[,time_col]) >= as.Date(ContinuousTS_2(current_dat,gap = 3*365)),]


          Start_Year <-  format(ContinuousTS_2(current_dat,gap = 3*365),"%Y")
          End_Year <-  format(max(current_dat[,time_col]),"%Y")

          Data_Start <- min(current_dat[,time_col])
          Data_End <- max(current_dat[,time_col])



        }else{
          current_dat <-current_dat[as.Date(current_dat[,time_col]) >= as.Date(max(current_dat[,time_col]))-(365*as.numeric(Timeframe)),]

          End_Year <-  format(max(current_dat[,time_col]),"%Y")
          Start_Year <- format(max(current_dat[,time_col])- years(as.numeric(Timeframe)),"%Y")
          Data_Start <- min(current_dat[,time_col])
          Data_End <- max(current_dat[,time_col])

        }


        if(p %in% log_parameters){current_dat[,3] <- log10(current_dat[,3]+1)}


        current_dat$myDate <- as.Date(as.character(current_dat[,time_col]),"%Y-%m-%d")
        current_dat <- GetMoreDateInfo(current_dat)


        if(season == "month"|season=="MONTH"|season=="Month"){

          current_dat$Season<-current_dat$Month
          SeasonString<-levels(current_dat$Season)
        }else if(season == "quarter"|season=="QUARTER"|season=="Quarter"){
          current_dat$Season<-current_dat$Qtr
          SeasonString<-levels(current_dat$Season)
        }else{
          current_dat$Season<-current_dat$CustomYear
          SeasonString<-levels(current_dat$Season)
        }


        NewValues <- RemoveAlphaDetect(current_dat[,p])
        current_dat <- cbind.data.frame(current_dat, NewValues)



        #make sure have both datasets
        current_dat <- current_dat[complete.cases(current_dat),]

        #Check for enough data to run analysis (20 values)

        Run_Analysis <- ifelse(nrow(current_dat[,1:3])>=20,TRUE,FALSE)


        if(Run_Analysis == T){

          pst <- SeasonalityTest(current_dat,main=paste(i, ": Raw Monthly Data",sep=""),do.plot = F)

          if(pst$pvalue<0.05){
            Trend_Analysis <- SeasonalTrendAnalysis(current_dat,ValuesToUse="RawValue",mymain=paste(i,"_",p,"_Raw Data: Seasonal Trend Analysis",sep=""),doPlot=F)
            Trend_Analysis$Flow_ADJ <- "NO"
            Trend_Analysis$Flow_MOD <- "N/A"
            Trend_Analysis$Flow_MOD_p <- "N/A"
            Trend_Analysis$Season_ADJ <- "YES"
            Trend_Analysis$Season_ADJ_p <- pst$pvalue
            #Trend_Analysis <- Trend_Analysis[,c(1:8,10:18,9,19:29)]
          }else{
            Trend_Analysis <- NonSeasonalTrendAnalysis(current_dat,ValuesToUse="RawValue",mymain=paste(i,"_",p,"_Raw Data: Non-Seasonal Trend Analysis",sep=""),doPlot=F)
            Trend_Analysis$Flow_ADJ <- "NO"
            Trend_Analysis$Flow_MOD <- "N/A"
            Trend_Analysis$Flow_MOD_p <- "N/A"
            Trend_Analysis$Season_ADJ <- "NO"
            Trend_Analysis$Season_ADJ_p <- "N/A"
          }

          Trend_Analysis$Site <- i
          Trend_Analysis$Param <- p
          Trend_Analysis$Start_Year <-Start_Year
          Trend_Analysis$End_Year <- End_Year
          Trend_Analysis$Time_Period <- Timeframe
          Trend_Analysis$Data_Start <- Data_Start
          Trend_Analysis$Data_End <- Data_End
          Trend_Analysis$Season <- season

        }else{
          Trend_Analysis <- data.frame( "nObs" = "INVALID DATASET","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                        "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                        "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                        "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA,
                                        "Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"=NA,"Season_ADJ"= NA, "Season_ADJ_p"= NA,"Site"= NA,"Param"= NA,"Start_Year"= NA,"End_Year"= NA,"Data_Start"=NA,"Data_End"=NA,"Season"=NA)

          Trend_Analysis$Site <- i
          Trend_Analysis$Param <- p
          Trend_Analysis$Start_Year <-Start_Year
          Trend_Analysis$End_Year <- End_Year
          Trend_Analysis$Time_Period <- Timeframe
          Trend_Analysis$Season <- season

        }

        Final_Df <- rbind(Final_Df,Trend_Analysis)

      }else{

        Trend_Analysis <- data.frame( "nObs" = "INVALID DATASET","S"= NA,"VarS"= NA,"D"= NA,"tau"= NA,"Z"=NA,"p"= NA,"Probability"= NA,
                                      "prop.censored"= NA,"prop.unique"= NA,"no.censorlevels"= NA,"Median"= NA,"Sen_VarS"= NA,"AnnualSenSlope"= NA,
                                      "Intercept"= NA,"Lci"= NA,"Uci"= NA,"AnalysisNote"= NA,"Sen_Probability"= NA,"Probabilitymax"= NA,"Probabilitymin"= NA,
                                      "Percent.annual.change"= NA,"TrendCategory"= NA,"TrendDirection"= NA,
                                      "Flow_ADJ"= NA,"Flow_MOD"= NA,"Flow_MOD_p"=NA,"Season_ADJ"= NA, "Season_ADJ_p"= NA,"Site"= NA,"Param"= NA,"Start_Year"= NA,"End_Year"= NA,"Data_Start"=NA,"Data_End"=NA,"Season"=NA)

        Trend_Analysis$Site <- i
        Trend_Analysis$Param <- p
        Trend_Analysis$Start_Year <-Start_Year
        Trend_Analysis$End_Year <- End_Year
        Trend_Analysis$Time_Period <- Timeframe
        Trend_Analysis$Season <- season

      }

    }
  }
  return(Final_Df)
}

###############################################################

ImprovementConfCat_BOPRC <- function (x, Reverse = c("CLAR","Final_Clarity" ,"MCI","QMCI","ASPM","Bottom_DO","MidHyp_DO",
                                                  "Native"))
{
  P <- x$Probability
  if (!is.na(Reverse[1])) {
    P[x$npID %in% Reverse] <- 1 - P[x$npID %in% Reverse]
  }
  ConfCats <- cut(P, breaks = c(-0.01, 0.01, 0.05, 0.1, 0.33,
                                0.67, 0.9, 0.95, 0.99, 1.01), labels = c("Exceptionally unlikely",
                                                                         "Extremely unlikely", "Very unlikely", "Unlikely",
                                                                         "As likely as not", "Likely", "Very likely",
                                                                         "Extremely likely", "Virtually certain"))
  ConfCats <- as.character(ConfCats)
  ConfCats[is.na(ConfCats)] <- "Not Analysed"
  ConfCats <- factor(ConfCats, levels = c("Exceptionally unlikely",
                                          "Extremely unlikely", "Very unlikely", "Unlikely",
                                          "As likely as not", "Likely", "Very likely",
                                          "Extremely likely", "Virtually certain",
                                          "Not Analysed"))
  return(ConfCats)
}
