\name{Trends_for_Aggregation_RAW}
\alias{Trends_for_Aggregation_RAW}
\title{Automated Trend Analysis without Flow Data}
\usage{
Trends_for_Aggregation_RAW(Data,param_colnames,siteID_col,
                                       time_col = "Time",
                                       Timeframe="All",
                                       season = "Quarter")
}
\description{
This function does not require flow data.  Data is seasonally adjusted based on ADJ_p <0.05.  This is useful for lake and estuarine data etc.
}
\examples{

}
