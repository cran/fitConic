\name{rotateA}
\alias{rotateA}
\alias{derotateA}
\alias{xyrot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rotate Conic Section Equation Parameters Or A Dataset, With Respect To X-Y Axes.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{rotateA} Takes as input "parA," the 6 values of the general quadratic Ax^2 + Bxy + Cy^2 +Dx + Ey +F = 0 , and applies a rotation angle to the coefficient set. 
\code{derotateA} calculates the rotation angle required to change the conic section defined by 'parA' into one that is orthogonal to the cartesian axes. 
\code{xyrot} is a simple function to rotate the coordinate system by theta. 
}
\usage{
rotateA(parA, theta)
derotateA(parA, ACmin = 1e-05)
xyrot(x, y = NULL, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parA}{
%%     ~~Describe \code{parA} here~~
the 6 values of the general quadratic Ax^2 + Bxy + Cy^2 +Dx + Ey +F = 0 
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
the angle, in radians, to rotate the conic section. 
}
  \item{ACmin}{
A tolerance parameter for deciding that the product of 
 parameters A and C is actually zero (in which case the 
 type of conic section is more likely a parabola or a degenerate case)
}
  \item{x}{
 Either a vector of x-coordinates or a Nx2 array of x and 
 y coordinates, in which case the y-input is ignored
}
  \item{y}{
 A vector of y-coordinates.
}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
\code{derotateA} uses the following standard formula to calculate the angle.
Derotate means to remove the xy term, i.e. force B = 0 . Some algebra shows that  cot(2theta) = (A-C)/B  and thus tan(2theta) = B/(A-C)

For \code{xyrot}, the internal \code{xy.coords} is used. If you enter only a vector for x and nothing for y, this will feed the new vectors 1:N for x and x-input for y to the rotator, which is probably not useful. 
}
\value{
For \code{derotateA}, 
\item{parA}{the new 6-parameter set defining the derotated conic.}
\item{theta}{the derived angle by which the parameter set was rotated}
For \code{rotateA}
\item{parA}{the new 6-parameter set defining the rotated conic.}
\item{theta}{the angle by which the parameter set was rotated}
For \code{xyrot} {a Nx2 array of the x,y coordinates of the rotated data set. }

}

\author{
 Carl Witthoft,  \email{carl@witthoft.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createConic}}
%% ~~objects to See Also as \code{\link{createConic}}, ~~~
}
\examples{
# make an ellipse and derotate it
parGr <- c(-2.3,4.2,5,3,pi/4)
xe <-seq(-8,9,by=.05)
elipGr <- createConic(xe,  parGr, 'e')
plot(elipGr, t= 'l',  asp = TRUE)
# convert to ABCDEF form
parAr <- GtoA(parGr,'e')
elipAr <- createConic(xe,parAr$parA) 
points(elipAr,pch='.',col='red')
# remove rotation angle
 parAd <- derotateA(parAr$parA)
 # returns theta = pi/4, how much the ellipse had been rotated by
 elipAd <-createConic(xe,parAd$parA)
lines(elipAd)
# rotate back
parAdr <- rotateA(parAd$parA, parAd$theta)
elipAdr <-createConic(xe,parAdr$parA)
lines(elipAdr,lty=3, lwd = 3, col='green')


}

