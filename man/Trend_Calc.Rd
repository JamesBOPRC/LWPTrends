\name{Trend_Calc}
\alias{Trend_Calc}
\title{Basic Trend Calculation}
\usage{
Trend_Calc(Data,param_colnames,siteID_col,
                       time_col = "Time",
                       Timeframe="All",
                       Seasonal = F,
                       Flow_Adjust=F,
                       flow_col=NA)
}
\description{
This trend analysis function performs trend analysis in the same way as all the other available functions, however the user can only choose to apply co-variate and seasonal adjustment to the entire dataset, i.e., all data is either adjusted or non-adjusted.
}
\examples{

}
