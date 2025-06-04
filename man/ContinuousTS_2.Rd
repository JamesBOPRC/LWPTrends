\name{ContinuousTS_2}
\alias{ContinuousTS_2}
\title{Calculate Date of Timeseries Gap}
\usage{
ContinuousTS_2(data,gap)
}
\description{
Timeseries datasets often have gaps in them that can affect trend analysis.  This function will assess the timeseries and return the date of the first sample collected after a break of a specified duration.

For example, you may decide that three years is the maximum allowable gap.  This function will return the latest date where any gaps are less than three years in duration.
}
\examples{

}
