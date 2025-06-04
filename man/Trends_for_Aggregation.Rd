\name{Trends_for_Aggregation}
\alias{Trends_for_Aggregation}
\title{Automated Trend Analysis}
\usage{
Trends_for_Aggregation(Data,param_colnames,siteID_col,
                                     time_col = "Time",
                                     Timeframe="All",
                                     flow_col=NA,
                                     season = "Quarter")
}
\description{
Automated trend analysis that flow adjustes if the FlowMod p-value is < 0.05, and seasonally adjusts if ADJ_p <0.05.  This function was used for the Rivers Report trend aggregation.
}
\examples{

}
