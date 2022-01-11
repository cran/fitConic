\name{doWeights}
\alias{doWeights}

\title{
 Function to Apply Weights to Data
}
\description{
This function applies an integer weight set to an array of (x,y) data points.  It normally is only called from \code{fitConic} but can be applied directly to a dataset if desired. 
 
}
\usage{
doWeights(XY, weights)
}

\arguments{
  \item{XY}{
A Nx2 array of data representing (x,y) pairs
}

\item{weights}{
 A vector of weights the same length as the number of rows in XY. At this time, only nonnegative integer values are allowed.  Doubles are rounded and negative values are set to zero. A zero weight will remove the matching data value from the dataset. 
}
}
\value{
A new Nx2 array. Basically, each row in the input XY is repeated \code{weights[j]} times. 
..
}

\author{
Carl Witthoft <carl@witthoft.com>
%%  ~~who you are~~
}

