\name{fhyp}
\alias{fhyp}
\alias{fhypopt}
\alias{costparab}
\alias{costparabxy}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal Functions to Perform Bootstrap Fitting Operations
%%  ~~function to do ... ~~
}
\description{
These functions are not intended for external use. 
\code{fhyp} and \code{fhypopt} support the parent function \code{bootHyperbola} by providing functions for \code{optimize} to use.
The functions \code{costparab} \code{costparabxy} similarly provide functions for \code{optim} to use inside the function \code{fitParabola} .

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fhyp(xy, b3, Ang)
costparabxy(theta, xy)
costparab(theta, xy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{
A Nx2 array of data
%%     ~~Describe \code{xy} here~~
}
  \item{b3}{
  Three of the parameters describing a hyperbola. These three are the "other parameters" fed to \code{optim}
%%     ~~Describe \code{b3} here~~
}
  \item{Ang}{
The initial angle of rotation, also optimized during the process.
%%     ~~Describe \code{Ang} here~~
}
  \item{theta}{
  The angle of rotation of the parabola for this run of \code{\link[stats]{optimize}}
%%     ~~Describe \code{theta} here~~
}
}

\value{
various combinations of "cost" values, i.e. Figure of Merit, to determine optimal set of coefficients, along with datasets where necessary. 
..
}

\author{
Carl Witthoft <carl@witthoft.com>
%%  ~~who you are~~
}

