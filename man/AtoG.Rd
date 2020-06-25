\name{AtoG}
\alias{AtoG}
\alias{GtoA}
\alias{parab3toA}
\alias{FEDtoA}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{  A Set Of Functions To Convert Among Various Conic-Section-Defining Parameter Sets. 
%%  ~~function to do ... ~~
}
\description{
\code{AtoG} Convert from full quadratic "ABCDEF"  to focus, axis, angle "hvab theta" parameters.
\code{GtoA} Convert from "hvab theta" to "ABCDEF" parameters.
\code{parab3toA} Simple conversion from  a + bx + cx^2 to "ABCDEF" parameters.
\code{FEDtoA} Convert focus, eccentricity, and directrix to "ABCDEF" parameters. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AtoG(parA, tol = 1e-06)
GtoA(parG, conicType = c("e", "h"))
parab3toA(ADF, theta = 0)
FEDtoA(focus = c(0, 0), directrix = c(1, 0, 1), eccentricity = 0.5)
}

\arguments{
  \item{parA}{ The six coefficients in the quadratic Ax^2 + Bxy + Cy^2 +Dx + Ey +F = 0 
}
  \item{tol}{ A small value, used to check whether small coefficient values might be actually zero. See "Details."  
}
  \item{parG}{ a five-element vector "h,v,a,b,theta" . See "Details" for the standard equation form for this. 
%%     ~~Describe \code{parG} here~~
}
  \item{conicType}{
%%     ~~Describe \code{conicType} here~~
Because the 'hvab' equation has a sign difference for ellipses vs. hyperbolas, it is necessary to indicate which kind of input is intended.  See "Details." 
}
  \item{focus}{ location of the conic sections focus.
%%     ~~Describe \code{focus} here~~
}
  \item{directrix}{ the 3-element directrix. 
%%     ~~Describe \code{directrix} here~~
}
  \item{eccentricity}{ the eccentricity of the conic section. 
%%     ~~Describe \code{eccentricity} here~~
}
  \item{ADF}{ The A,D,F coeffients in the standard quadratic. Thus, the x^2 term, the x term, and the constant term.
%%     ~~Describe \code{ADF} here~~
}
  \item{theta}{ An angle by which the entire parabola is to be rotated. 
%%     ~~Describe \code{theta} here~~
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~

The \code{tol} input for  \code{AtoG} checks two conditions. First, is B practically zero, in which case B is set to exactly zero, implying no rotation of the conic section. Second, is B^2 - 4*A*C almost zero, implying that the conic is probably a parabola, and conversion to 'hvab' form is not useful.

The "hvab" form for describing an ellipse or a hyperbola looks like [Center(1:2), Axes(1:2)/2] angle A, to fill the equation 
  
 ((x-h)cosA +(y-v)sinA)^2/a^2 + ((x-h)sinA-(y-v)cosA)^2/b^2 = 1 
The length of the axes are 2*a, 2*b .

A discussion of the focus/directrix/eccentricity form of a conic section is rather lengthy and not presented here.  One short introduction can be found at 
\url{https://en.wikipedia.org/wiki/Conic_section#Eccentricity,_focus_and_directrix}

}
\value{
for \code{AtoG}, 
\item{parG}{ c(h,v,a,b,theta)}
\item{exitCode}{a value used in fitConic. 1,2, or 3 for ellipse, hyperbola, parabola}
\item{conicType}{matching exitCode with a char "e", "h", or "p" }
for \code{GtoA}
\item{parA}{the ABCDEF coefficients of the general quadratic}
\item{exitCode}{ a value used in fitConic. 1,2, or 3 for ellipse, hyperbola, parabola}
\item{conicType}{matching exitCode with a char "e", "h", or "p" }
for \code{FEDtoA}, the ABCDEF coefficients of the general quadratic
for \code{parab3toA},
\item{parA}{the ABCDEF coefficients of the general quadratic}
\item{exitCode}{always numeric 3,  a value used in fitConic}
	\item{conicType}{always char  "p"}
.
}


