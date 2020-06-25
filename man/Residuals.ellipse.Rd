\name{Residuals.ellipse}
\alias{Residuals.ellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Residual Error For Current Coefficients
%%  ~~function to do ... ~~
}
\description{
This function is not intended for external use. It is called from \code{fitConic} when iterating to find the best-fit ellipse. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Residuals.ellipse(XY, parG)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY}{
  The x,y dataset
%%     ~~Describe \code{XY} here~~
}
  \item{parG}{
  The "G-parameter" set for the current iteration.
%%     ~~Describe \code{parG} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{RSS }{Figure of merit, the 'norm' of the difference between the input XY data and the output "XYproj" data generated. }
 \item{XYproj }{ Calculated dataset to be used in generating the  Jacobian matrix for the next iteration of \code{fitConic}  }
%% ...
}
\references{
This is a slightly modified (and debugged) version of \code{\link[conicfit]{Residuals.ellipse} }
%% ~put references to the literature/web site here ~
}
\author{
Carl Witthoft <carl@witthoft.com>
%%  ~~who you are~~
}

