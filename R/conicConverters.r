# Conversion of geometric parameters of an ellipse
# parG <- [Center(1:2), Axes(1:2), Angle]'
# to algebraic parameters
# have to tell it whether ellipse or hyperbola , i.e. x^2 +/- y^2 = 1
# Not applicable to parabolas. 
GtoA<-function(parG, conicType = c('e','h'))  {
if (length(parG) != 5) {
	stop('parG must have five values')
}
conicType = tolower(conicType)
if (conicType == 'e') {
	thesign = 1
	exitCode = 1
}   else {
		if (conicType == 'h') {
			thesign = -1
			exitCode  = 2
			} else stop( 'conicType must be "e" or "h"')
	}
if(length(parG) != 5) stop('parG must be 5 elements, has ',length(parG)) 

k <- cos(parG[5])
s <- sin(parG[5])
a <- parG[3]
b <- parG[4]
h <- parG[1]
v <- parG[2]
# this is same as Wikipedia ellipse page if you divide thru by (a^2*b^2)
P <- (k/a)^2 + thesign*(s/b)^2
Qmat <- (s/a)^2 + thesign*(k/b)^2
R <- 2*k*s*(1/a^2 - thesign*1/b^2)
#  the factor of 2 is hidden inside "R" terms
# don't cbind() because that makes a matrix
parA <- c(A=P, B=R, C=Qmat, D=-2*P*h-R*v, E=-2*Qmat*v-R*h, F=P*h^2+Qmat*v^2+R*h*v-1)
return(invisible(list(parA=parA, conicType = conicType, exitCode = exitCode))) 
}


# convert  standard parabola quadratic to parA with optional rotation added
#output is ABCDEF, six values
parab3toA <-function(ADF,theta = 0){
# check input
if (length(ADF) != 3) {
	stop(c('ADF must be 3 coeffs Ax^2 + Dx + F but is ',length(ADF),' long'))
	}
sintheta <- sin(-theta)
costheta <- cos(-theta) 
# remember we push Y to left side, so E = -1 before rotating
# rotate as user expects, i.e. not xyrot direction
# so x-prime = x*cost y*sintheta, y-prime = x*sintheta + y*cost
A = ADF[1]*costheta^2
A[2] = -2*ADF[1]*costheta*sintheta
A[3] = ADF[1]*sintheta^2
A[4] = ADF[2]*(costheta -sintheta)
A[5] = -1 *( ADF[2]*sintheta + costheta)
A[6] = ADF[3]
return(invisible(list(parA = A, conicType = 'p', exitCode = 3) ))
}



#   formulas based on equations in Wikipedia ellipse page. 
# check for parabola, i.e. B^2 == 4AC to within "tol"
# 'tol' also used to check B == 0, i.e. not rotated
AtoG <-function(parA, tol = 1e-6) {
#check for 6 values 
if (length(parA) !=6) stop('parA must have length 6')
A = parA[1]
B = parA[2]
C = parA[3]
D = parA[4]
E = parA[5]
F = parA[6]
# for safety, esp. in finding theta
if (abs(B) < tol) B = 0

# useful intermediate values
radical =  B^2 - 4*A*C
# sign flip for hyperbola vs. ellipse
ACB <- sqrt( (A-C)^2 + B^2) * -sign(radical)
# for safety,
if (abs(radical) < tol) radical <- 0 # practically a parabola
# for compatibility with fitConic()  generate exitCode
if (radical < 0){
	conicType = 'e'
	exitCode = 1
	}  else {
		if (radical > 0) {
		exitCode = 2
		conicType = 'h'
		} else {
		exitCode = 3 #parabola; don't try to convert
		conicType = 'p'
		}
	}
# only do work when radical !=0
if (radical != 0) {
#	run thru  equations for a,b,h,v 
	abone <- 2*(A*E^2 + C*D^2 - B*D*E + F * radical)
	abtwo <- A + C
	a <- sign(radical) * sqrt(abs(abone * ( abtwo + ACB )) )/radical
	b <-  sign(radical) * sqrt(abs(abone * ( abtwo - ACB )) )/radical
	h <- (2*C*D - B*E)/radical
	v <- (2*A*E - B*D)/radical
} else {
	#  parabola, so return nothing
		h <- v <- a <-b <- NaN
	}
# now do theta 
if(B ==0) theta = pi/2*(1- ((A-C) < 0 ) ) else{
	theta = atan( (C - A - ACB)/B ) 
	}
return(invisible(list(parG = c(h,v,a,b,theta), exitCode = exitCode, conicType =  conicType)))
}



# Derotate means to remove the xy term, i.e. force B = 0, 
# find theta:  cot(2theta) = (A-C)/B  
# tan(2theta) = B/(A-C)
# cos(theta) = cos(atan(2theta)/2)
derotateA <-function(parA, ACmin = 1e-5){
if (length(parA) != 6){
	stop( c('parA must have exactly 6 elements, has ',length(parA)) )
	}
#TODO: put in checks for  other things that make angle calc'ns blow up due to input being unrotated or degenerate conic 
if ( abs(parA[1]-parA[3]) < ACmin) {
	warning('A -C near zero, may be degenerate case')
	}
# wolfram says (at least for ellipse), if A>C add pi/2 to theta. I don't see
# that working here.  
# flip sign of theta depending on parity of A, C 
elhyp <- -sign(parA[2]^2 - 4*parA[1]*parA[3]) 
if (elhyp == 0 ) elhyp = 1  # sneak parabola thru
theta = elhyp* atan((parA[2])/(parA[1]-parA[3]))/2
tcos = cos(theta) #  + pi/2*(parA[1]>parA[3]) 
tsin = sin(theta) # + pi/2*(parA[1]>parA[3]) 
Qout <- parA[1]*tcos^2 + parA[2] * tsin * tcos + parA[3] * tsin^2
Qout[2] <- -2*parA[1] * tsin * tcos - parA[2] * tsin^2 + parA[2] * tcos^2 + 2*parA[3] * tsin * tcos

Qout[3] <- parA[1]*tsin^2 - parA[2]*tsin*tcos + parA[3] * tcos^2
Qout[4] <- parA[4]*tcos + parA[5] * tsin
Qout[5] <- -parA[4]*tsin + parA[5]*tcos 
Qout[6] =   parA[6]
#cleanup
names(theta) <- NULL
names(Qout) <- c("A","B","C","D","E","F")
return(invisible(list(parA = Qout, theta = theta ) ))
}

# apply arbitrary rotation to any conic function
rotateA <- function(parA, theta) {
if (length(parA) != 6){
	stop( c('parA must have exactly 6 elements, has ',length(parA)) )
	}
theta = -theta[1] # just in case, and "user expectation" of direction of rotation
tcos = cos(theta) # + pi/2*(parA[1]>parA[3]) 
tsin = sin(theta) #+ pi/2*(parA[1]>parA[3]) 
Qout <- parA[1]*tcos^2 + parA[2] * tsin * tcos + parA[3] * tsin^2
Qout[2] <- -2*parA[1] * tsin * tcos - parA[2] * tsin^2 + parA[2] * tcos^2 + 2*parA[3] * tsin * tcos

Qout[3] <- parA[1]*tsin^2 - parA[2]*tsin*tcos + parA[3] * tcos^2
Qout[4] <- parA[4]*tcos + parA[5] * tsin
Qout[5] <- -parA[4]*tsin + parA[5]*tcos 
Qout[6] =   parA[6]
return(invisible(list(parA = Qout, theta = theta ) ))
}


#Convert from focus-directrix-eccentricity to parA form 
#  ported from info at 
#  https://math.stackexchange.com/questions/2800817
# where generalized directrix D =  ax + by +c = 0 (different a,b,c) and
# focus (xf,yf) , leads to  (e is eccentricity here)
# (x-xf)^2 + (y-yf)^2 = e^2*(D^2)/(a^2+b^2) which can "easily" be converted into parA
# the parA solution interms of focus & directrix & e is
 #  t = a * a + b * b; 
    # A = t - e^2 * (a * a); 
    #  C = t - e^2 * (b * b); 
    #  D = (-2 * t * h) - (2 * e^2 * c * a); 
    #  E = (-2 * t * v) - (2 * e^2 * c * b); 
    # B = -2 * e^2 * a * b; 
    # F = (-e^2 * c * c) + (t * h * h) + (t * v * v); 
FEDtoA <-function(focus = c(0,0), directrix = c(1,0,1), eccentricity = 0.5 ) {
h = focus[1]
v = focus[2]
da = directrix[1]
db = directrix[2]
dc = directrix[3]
ec = eccentricity^2
# sign flip from GFG page
k = (da^2 + db^2)
parA = k - ec*da^2 # A term
parA[2] = -2*ec*da*db  # B term,  and so on
parA[3] = k -ec*db^2
parA[4] = -2*h*k - 2*ec*da*dc
parA[5] = -2*v*k - 2*ec*db*dc
# if dc is zero get degenerate case because F is zero? yes -- not a bug. 
parA[6] = -ec*dc^2 + k*(h^2 + v^2)

return(invisible(parA)) 
}

# have 6 unknowns and 6 eqns so it is invertible...
# May be added in future releases
# AtoFED <-function(parA) {
# }

#need this for a couple funcs
xyrot<-function(x, y=NULL, theta){
	# pairs must be Nx2 matrix w/ x in first column and y in second
	xy <-xy.coords(x,y)
	xrot <- xy$x*cos(theta) + xy$y*sin(theta)
	yrot <- -xy$x*sin(theta) + xy$y*cos(theta)
	return(invisible(cbind(xrot,yrot)))
}