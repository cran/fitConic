\name{bootHyperbola}
\alias{bootHyperbola}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Attempt a Crude Fit of Data to a Hyperbola
%%  ~~function to do ... ~~
}
\description{
This function is not intended for direct use. It attempts to generate an approximate fit of a data set to a hyperbola, returning a parameter set for use in intializing the main function \code{conicFit} . 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootHyperbola(x, y = NULL, maxiter = 10000, ...)
}

\arguments{
  \item{x}{ vector of x-values, or a Nx2 array of x and y values. In the latter case, the input \code{y} is ignored. 
}
  \item{y}{ vector of y-values. 
}
  \item{maxiter}{ A 'safety' limiter on the number of iterations to try before giving up.

}
  \item{\dots}{ possible other arguments to be passed to future upgrades
}
}

\value{
\item{parA}{the new 6-parameter set defining the non-rotated conic.}
\item{parAr}{the new 6-parameter set defining the rotated conic.}
\item{theta}{the angle of rotation between ParA and ParAr}
\item{fitdat}{the information returned from \code{\link[stats]{optim}}}
}
\author{
%%  ~~who you are~~
 Carl Witthoft,  \email{carl@witthoft.com}
}

\seealso{
\code{\link[fitConic]{fitConic}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

