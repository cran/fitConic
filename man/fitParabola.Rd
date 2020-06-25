\name{fitParabola}
\alias{fitParabola}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Data to Parabola 
%%  ~~function to do ... ~~
}
\description{
This function fits a data set to a parabola, including any rotation angle. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitParabola(x, y = NULL, searchAngle = c(-pi/2, pi/2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{ vector of x-values, or a Nx2 array of x and y values. In the latter case, the input \code{y} is ignored. 
%%     ~~Describe \code{x} here~~
}
  \item{y}{ vector of y-values. 
%%     ~~Describe \code{y} here~~
}
  \item{searchAngle}{
  Optional pair of angles, in radians, defining the limits of the search range to find the rotation angle of the parabola. Usually the default range \code{-pi/2:+pi/2} works acceptable. 
%%     ~~Describe \code{searchAngle} here~~
}
  \item{\dots}{
  For possible future expansion to pass to additional features. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{
\code{fitParabola} starts by doing a RANSAC-style search to find the optimum rotation angle. Once that is chosen, the data are rotated by that angle and a simple polynomial fit to the (rotated) vertical parabola is done.
 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{vertex }{ calculated vertex of the parabola}
 \item{theta }{angle of rotation relative to a vertical parabola }
 \item{parA }{ the "ABCDEF" coefficients of the fitted parabola}
 \item{parQ }{ the coefficients of the derotated parabola's simple quadratic polynomial, highest power first}
  \item{cost }{ final value of the "cost" parameter used for optimization}
}
\references{
Some of the code is based on \url{http://www.mathworks.com/matlabcentral/answers/80541}
%% ~put references to the literature/web site here ~
}
\author{
Carl Witthoft <carl@witthoft.com>

%%  ~~who you are~~
}
\note{
When the function \code{fitConic} is called with instructions to fit to a parabola, it passes the inputs to \code{fitParabola} and does nothing else. For parabolic data, then, either function will give the same result. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[fitConic]{createConic}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Create vertical parabola with some noise
parP <-c(.5,0,0,2,-1,4)
xp <- seq(-5,5,by=0.05)
 partest <-createConic(xp,param = parP,ranFun = rnorm, noise = 1)
 plot(partest, pch= '.',asp=TRUE, cex=3)
 # rotate the data
 partestr <-xyrot(partest,theta = -.35)
 points(partestr,col='green',pch='.',cex=3)
 # do the fit
 parfit <-fitParabola(partestr)
 points(parfit$vertex,pch='X',col='blue') 
 parout <- createConic(xp,parfit$parA)
 lines(parout,col='red')

}
