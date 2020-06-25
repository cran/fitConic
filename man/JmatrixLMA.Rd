\name{JmatrixLMA}
\alias{JmatrixLMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate a Jacobian Matrix
%%  ~~function to do ... ~~
}
\description{
Calculate the Jacobian matrix with the original dataset and the
 current version of fitted data.  This is not intended for 
 external use. It is called from \code{fitConic}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
JmatrixLMA(XY, parA, XYproj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{
  The original input dataset
%%     ~~Describe \code{XY} here~~
}
  \item{parA}{
  The current set of ABCDEF quadratic equation coefficients.
%%     ~~Describe \code{parA} here~~
}
  \item{XYproj}{
  The current calculated dataset based on the latest iteration
   of the coefficient set. 
%%     ~~Describe \code{XYproj} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{Res }{ residuals based on the norm of XY - XYproj}
 \item{J }{matrix of values for each input data point corresponding
 to the terms in the general quadratic Ax^2 + Bxy + Cy^2 +Dx + Ey +F }
%% ...
}
\references{
This is a copy of \code{\link[conicfit]{JmatrixLMA} } with some
  validation steps added. 

%% ~put references to the literature/web site here ~
}
\author{
Carl Witthoft <carl@witthoft.com>
%%  ~~who you are~~
}
