\name{fitConic}
\alias{fitConic}
\title{
Fit Data to A Conic Section Curve 
}
\description{
This function fits data to an ellipse, hyperbola, or parabola.
 It can do so without any initial conditions, or can accept
  initial parameter values when known.  
}
\usage{
fitConic(X, Y = NULL, parInit = NULL, conicType = c("e", "h", "p"),
 LambdaIni = 1, epsilonP = 1e-06, epsilonF = 1e-06, IterMAX = 20000)
}
\arguments{
  \item{X}{ vector of x-values, or a Nx2 array of x and y values.
   In the latter case, the input \code{y} is ignored. 
}
  \item{Y}{ vector of y-values. 
}
  \item{parInit}{
Optional. A vector either of six values representing an initial guess
 at the "ABCDEF" coefficients of the quadratic, or five values
  representing an initial guess at the "hvab,theta" coefficients.
  In the latter case, a value of either "e" or "h" is required for
 \code{conicType}.  See the Details section for more information. 
}
  \item{conicType}{
If \code{parInit} is either NULL or the "hvab,theta" option,
 \code{conicType} is required. Enter either "e", "h", or "p" for
  fitting to ellipse, hyperbola, or parabola. 
}
  \item{LambdaIni}{
A control parameter used in the fitting algorithm. Typically
 there is no reason to change from the default value. 
%%     ~~Describe \code{LambdaIni} here~~
}
  \item{epsilonP}{
A tolerance value to determine whether convergence has occurred. 
  %%     ~~Describe \code{epsilonP} here~~
}
  \item{epsilonF}{
    A tolerance parameter for determining when to adjust lambda
     away from the input value \code{LambdaIni}.

}
  \item{IterMAX}{
A "safety" value to avoid loop thrashing when convergence isn't taking place. 

}
}
\details{
  \code{ParInit}, when supplied is either a 6-value set representing the standard quadratic form  Ax^2 + Bxy + Cy^2 +Dx + Ey +F = 0  or a 5-value set representing the "hvab,theta" form ((x-h)cosA +(y-v)sinA)^2/a^2 + ((x-h)sinA-(y-v)cosA)^2/b^2 = 1 . In the latter case the value \code{conicType} is required, because ellipses and hyperbolas have a different sign for the y-term. 
In most cases, the bootstrapper tools work well enough to allow the main algorithm to fit to an ellipse or hyperbola. However, "knowledge is power." If you have a good idea approximately what the \code{ParIni} values are,
 entering them will help avoid convergence to the wrong local minimum. 
The algorithm branch which fits data to parabolas does not use or need
 initialization, as it uses a RANSAC-type search to find 
 the best rotation angle, and then does a simple quadratic polynomial fit. 

}
\value{

\code{parA }{vector of the six "ABCDEF" coefficients
\code{RSS }{'root sum square' figure of merit describing the relative fit quality}
\code{iters }{number of iterations at convergence}
\code{exitCode }{1 means ellipse, 2 means hyperbola, 3 means parabola.
 If other values show up (possibly -1, 0, 4), most likely the
  dataset led to a degenerate case such as a line fit. }

}
}
\references{
\url{https://people.cas.uab.edu/~mosya/cl/ } for information
 on the original "LMA" fitting algorithm. 
 \url{https://math.stackexchange.com/questions/426150} and 
 \url{https://math.stackexchange.com/questions/2800817 }
 for various related equations concerning conic sections.
 \url{https://en.wikipedia.org/wiki/Ellipse} for several parameter
 conversion formulas

}
\author{
Carl Witthoft <carl@witthoft.com>

}
\seealso{  
\code{\link[fitConic]{createConic}} , \code{\link[fitConic]{fitParabola} }
}
\examples{

##-create a hyperbola, add noise
Ang =  0.42 #radians
xh <- seq(-20,20,by=0.1)
parAxyh <- c(0, 1, 0, -2, 4, -15 )
parAxyhr <- rotateA(parAxyh, Ang)$parA
newxyr <-createConic(xh,parAxyhr)
newxyrn <- createConic(xh,parAxyhr,ranFun=rnorm, noise= 0.05)
plot(newxyr, t = 'l',asp=TRUE)
points(newxyrn, pch = '.', cex = 3)
# Now find the hyperbola for that dataset
hypfitr <-fitConic(newxyrn, conicType = 'h')
hypdatr <- createConic(xh, hypfitr$parA)
lines(hypdatr, col='red')
}
