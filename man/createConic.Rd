\name{createConic}
\alias{createConic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create A Conic Section Dataset Based on Parameter Set
%%  ~~function to do ... ~~
}
\description{
Given a vector of x-values and a parameter set defining a conic section, produce an array of x- and y- values, optionally with noise added, for the specified conic section. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createConic(x, param, conicType, ranFun = NULL, noise = 1, seedit = NULL, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Vector of (real) values 
%%     ~~Describe \code{x} here~~
}
  \item{param}{
  Either a 6-value set representing the standard quadratic form  Ax^2 + Bxy + Cy^2 +Dx + Ey +F = 0  or a 5-value set representing the "hvab,theta" form ((x-h)cosA +(y-v)sinA)^2/a^2 + ((x-h)sinA-(y-v)cosA)^2/b^2 = 1 . In the latter case the value \code{conicType} is required. 
%%     ~~Describe \code{param} here~~
}
  \item{conicType}{
  Either the character "e" for ellipse or "h" for hyperbola. Only required if the "hvab,theta" form is used in \code{param} . 
%%     ~~Describe \code{conicType} here~~
}
  \item{ranFun}{
  If random noise is to be added to the calculated y-values, provide a vectorized function which takes a single input (x). See Details. 
%%     ~~Describe \code{ranFun} here~~
}
  \item{noise}{
  Optional argument to multiply the output of \code{ranFun} .
%%     ~~Describe \code{noise} here~~
}
  \item{seedit}{
  Optional argument to set a starting seed for \code{ranFun} to use. 
%%     ~~Describe \code{seedit} here~~
}
  \item{tol}{
  A (small) value used to decide whether various parameter terms are so small that they should be zero. This is used to facilitate distinguishing, e.g., parabolas from hyperbolas. 
%%     ~~Describe \code{tol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

When supplied \code{ranFun} is used as follows.
y <- y + ranFun(y)*noise  . Make sure any function supplied fits that form (no other input argument required; only a vector returned).


}
\value{
An N x 2 array of the x,y pairs.  Warning: since there are often two possible y-values for a given x-value (these being quadratic equations), the array does contain duplicate x-values. This may "annoy" some other packages' functions which don't allow that sort of repeated value. If this presents a problem, I'd recommend applying a very small amount of noise to the x-values in this output. 

}

\author{
Carl Witthoft <carl@witthoft.com>
%%  ~~who you are~~
}


\examples{

# create noisy ellipse
parGr <- c(-2.3,4.2,5,3,pi/4)
xe <-seq(-8,9,by=.05)
elipGrn <- createConic(xe, parGr, 'e',ranFun=rnorm, noise=0.25)
elipGr <- createConic(xe,  parGr, 'e')
plot(elipGrn, pch='.',cex = 4, asp = TRUE)  #, xlim = c(-5,8), ylim = c(0,7))
lines(elipGr,col='green')
}

