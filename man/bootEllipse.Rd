\name{bootEllipse}
\alias{bootEllipse}

\title{
Simple, Medium-Quality Ellipse Fitting Function
}
\description{
This function generates a half-decent fit to the source data. It is intended only for internal use, to bootstrap the higher-quality \code{fitConic} function. 
}
\usage{
bootEllipse(x, y = NULL, ...)
}
\arguments{
  \item{x}{ vector of x-values, or a Nx2 array of x and y values. In the latter   case, the input \code{y} is ignored. 
}
  \item{y}{ vector of y-values. }
  \item{\dots}{ possible other arguments to be passed to future upgrades
}
}
\details{
This can be used as a Q&D ellipse fitting algorithm, but is intended only for internal use by \code{fitConic}, providing that function with an initial estimate for the ellipse's defining parameter set. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{parA }{6-element set with estimate of the "ABCDEF" coefficients for the general quadratic equation}
 \item{centroid }{estimate of the ellipse's centroid}
%% ...
}
\references{
This is a revision of the function EllipseDirectFit in package conicfit by Jose Gama, with minor upgrades. 
 Original MATLAB code by: Nikolai Chernov \url{http://www.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method}
A. W. Fitzgibbon, M. Pilu, R. B. Fisher, "Direct Least Squares Fitting of Ellipses", IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
 Halir R, Flusser J (1998) Proceedings of the 6th International Conference in Central Europe on Computer Graphics and Visualization, 
 Numerically stable direct least squares fitting of ellipses (WSCG, Plzen, Czech Republic), pp 125132.

}
\author{
 Carl Witthoft,  \email{carl@witthoft.com}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[fitConic]{fitConic}} , \code{\link[fitConic]{createConic} }
}


