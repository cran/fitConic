 
# See the solved equations at math.SE 
#
#  https://math.stackexchange.com/questions/426150
#  https://math.stackexchange.com/questions/2800817

# IMPORTANT ATTRIBUTION:
#  https://people.cas.uab.edu/~mosya/cl/  and the folks referred to there,
# for fitConicLMA .
#  https://en.wikipedia.org/wiki/Ellipse for several parameter conversion formulas
#
# G-parameter for hyperbola is only for transverse axis - otherwise change to  y-term minus x-term = 1.  
# ellipse is same G-param but sum the two terms. 
#
# G-form: center or vertex is always h,k . Axes are 2*a, 2*b . 
#  for angle A,   ((x-h)cosA +(y-v)sinA)^2/a^2 + ((x-h)sinA-(y-v)cosA)^2/b^2 = 1 covers rotated ellipse

# here's the rules for dif't conic sections:
# Ax^2+Bxy+Cy^2+Dx+Ey+F=0  where A,B,C,D,E and F are constants. 
# If B^2 - 4AC is less than zero, if a conic exists, it will be either a circle or an ellipse. 
# If B^2 - 4AC equals zero, if a conic exists, it will be a parabola. 
# If B^2 - 4AC is greater than zero, if a conic exists, it will be a hyperbola.
#if A = C and B = 0, the equation represents a circle
#
#   B = 0  means axes parallel to X,Y frame.
#

# Fitting a conic to a given set of points (Implicit method) using algebraic parameter
#epsilonP <- tolerance (small threshold)
#epsilonF <- tolerance (small threshold)
#IterMAX <- maximal number of (main) iterations usually 10-20 suffice

#allow parInit = NULL to force a call to a bootstrap. conicType is required
# in this option.  
fitConic <- function(X, Y = NULL, parInit=NULL, conicType=c('e','h', 'p'), LambdaIni=1, epsilonP= 1e-6, epsilonF=1e-6, IterMAX=2e4) {
#This is a highly modified version of fitConicLMA function in package conicfit by Jose Gama
newxy <- xy.coords(X,Y)
XY <-cbind(newxy$x,newxy$y)
# removing NaNs and Infs
goods <- is.finite(XY[ ,1]) & is.finite(XY[ ,2])
XY <- XY[goods,]

if (is.null(parInit) ) {
#		message(c('Bootstrapping type ',conicType[1]))
		switch(tolower(conicType[1] ),
		'e' = {
			fitdir <-bootEllipse(XY)
			parInit = fitdir$parA
			tmp <- AtoG(parInit)
			parGini <- tmp$parG
			parAini <- parInit
			exitCode = 1
		},
		'h' = {
			fitdir <-bootHyperbola(XY)
			parInit = fitdir$parAr
			tmp <- AtoG(parInit)
			parGini <- tmp$parG
			parAini <- parInit
			exitCode = 2
		},
		'p' = {
	#no separate bootstrap required
			exitCode <- 3			
		},
		stop('Unknown conicType: select e[llipse], h[yperbola], or p[arabola]')
		) #end switch
} else {
	# Input validation
	parInit = as.vector(parInit) 
	#test for type of initial parameter set
	if (length(parInit) == 6) {
	# is A-parameters
		tmp <- AtoG(parInit)
		parGini <- tmp$parG
		exitCode <- tmp$exitCode
		parAini <- parInit
	} else if (length(parInit) == 5) {
	# it's parG
		parAini <- as.vector( GtoA(parInit, conicType)$parA )
		parGini <- parInit
		exitCode <- AtoG(parAini)$exitCode
	}else stop('parInit must be 5 (G) or 6 (A) values')
} #end of ifelse parInit is NULL. 

lambda.sqrt <- sqrt(LambdaIni)   #  sqrt(Lambda) is actually used by the code

switch(as.character(exitCode),
	'1' = {
		tmp <- Residuals.ellipse(XY,parGini)
		Fp <- tmp$RSS
		XYproj <- tmp$XYproj  
	},
	'2' = {
		tmp <- Residuals.hyperbola(XY,parGini)
		Fp <- tmp$RSS
		XYproj <- tmp$XYproj 
	},
	'3' = {
# RANSAC-style fiting.  Seems to work better than some matrix methods. 
		pfit <- fitParabola(XY)
		return(invisible(pfit))  
	},
	stop(paste0('invalid conic type ', exitCode)) )
#end of switch 
# calculate the Jacobian matrix
tmp <- JmatrixLMA(XY,parAini,XYproj)  
  Res <- tmp$Res
  J <- tmp$J
parA <- parAini
parG <- parGini
codeTemp <- exitCode
parGTemp <- parG
#  main loop, each run is one (main) iteration
for (iter in 1:IterMAX) {   
#	if( !(iter%%10)) message(c('ten more iters,code is ',exitCode ) )
#	iwhile = 0
# init for diagnostics
#	progress = 10
    while (TRUE) {   #secondary loop - adjusting Lambda (no limit on cycles)

#if( !(iwhile%%10)) message(c('ten more whiles,code is ', exitCode,' progress ' ,progress))
#iwhile <- iwhile + 1
      DelPar <- pracma::mldivide(rbind(J, lambda.sqrt*diag(6)), rbind(-Res, matrix(0,6,1))) # step candidate
        ParTemp <- parA + DelPar
        ParTemp <- ParTemp/norm(ParTemp,'F')
        tmp <- AtoG(ParTemp)
        parGTemp <- tmp$parG
        codeTemp <- tmp$exitCode
        if (codeTemp != exitCode) progress <- 1 else progress <- norm(DelPar,'F')
        if (progress < epsilonP)  break               # stopping rule

		switch(as.character(codeTemp),
			'1' = {
				tmp <- Residuals.ellipse(XY,parGTemp)
				FTemp <- tmp$RSS
				XYprojTemp <- tmp$XYproj  
		},
			'2' = {
				tmp <- Residuals.hyperbola(XY,parGTemp)
				FTemp <- tmp$RSS
				XYprojTemp <- tmp$XYproj   
		},
		{
				lambda.sqrt <- lambda.sqrt*2 #if it's degenerate, increase lambda, recompute the step
				next
				}
		)  #end of switch  

		if (FTemp < Fp*(1-epsilonF/lambda.sqrt)){        #   yes, improvement
			lambda.sqrt <- lambda.sqrt/2   # reduce lambda, move to next iteration
			break
		}else {       #   no improvement
			lambda.sqrt <- lambda.sqrt*2 # increase lambda, recompute the step
			next
		}
    }   #  end of  while (TRUE)
    if (progress < epsilonP)  break # stopping rule
    tmp <- JmatrixLMA(XY,ParTemp,XYprojTemp)
    Res <- tmp$Res
    J <- tmp$J
    parA <- ParTemp
    Fp <- FTemp  # update the iteration
    parG <- parGTemp
    exitCode <- codeTemp
}      #    main loop
RSS <- Fp
iters <- iter
return (invisible(list(parA=as.vector(parA),RSS=RSS,iters=iters,exitCode=exitCode) ))
}#end of fitconic

# Ellipse bootstrapper.
bootEllipse <-function(x, y = NULL, ...){
#This is a copy of the function EllipseDirectFit in package conicfit by Jose Gama, with minor upgrades. 
# Original code by: Nikolai Chernov http://www.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method
# A. W. Fitzgibbon, M. Pilu, R. B. Fisher, "Direct Least Squares Fitting of Ellipses", IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
# Halir R, Flusser J (1998) Proceedings of the 6th International Conference in Central Europe on Computer Graphics and Visualization, 
# Numerically stable direct least squares fitting of ellipses (WSCG, Plzen, Czech Republic), pp 125 - 132.
xy <- xy.coords(x,y)
XY <-cbind(xy$x,xy$y)
centroid <- apply(XY,2,mean)
D1 <- cbind((XY[,1]-centroid[1])^2, (XY[,1]-centroid[1])*(XY[,2]-centroid[2]), (XY[,2]-centroid[2])^2)
D2 <- cbind(XY[,1]-centroid[1], XY[,2]-centroid[2], matrix(1,dim(XY)[1]))
S1 <- t(D1) %*% D1
S2 <- t(D1) %*% D2
S3 <- t(D2) %*% D2
T <- -solve(S3) %*% t(S2)
M <- S1 + S2 %*% T
M <- rbind(M[3,]/2, -M[2,], M[1,]/2)
evec<-eigen(M)$vectors
cond <- 4*evec[1,]*evec[3,]-evec[2,]^2
A1 <- matrix(evec[,which(cond>0)[1]],3)
A <- rbind(A1, T %*% A1)
A4 <- A[4]-2*A[1]*centroid[1]-A[2]*centroid[2]
A5 <- A[5]-2*A[3]*centroid[2]-A[2]*centroid[1]
A6 <- A[6]+A[1]*centroid[1]^2+A[3]*centroid[2]^2+ A[2]*centroid[1]*centroid[2]-A[4]*centroid[1]-A[5]*centroid[2]
A[4] <- A4;  A[5] <- A5;  A[6] <- A6
A <- A / norm(A,'F')
# # general-form conic equation  ax^2 + bxy + cy^2 +dx + ey + f = 0
# a <- A[1];b <- A[2]/2;C <- A[3];d <- A[4]/2;E <- A[5]/2;f <- A[6]
# x0 <- (C*d-b*E)/(b*b-a*C)
# y0 <- (a*E-b*d)/(b*b-a*C)
# semiaxis.a <- sqrt(2*(a*E^2+C*d^2+f*b^2-2*b*d*E-a*C*f)/((b^2-a*C)*(sqrt((a-C)^2+4*b^2)-(a+C))))
# semiaxis.b <- sqrt(2*(a*E^2+C*d^2+f*b^2-2*b*d*E-a*C*f)/((b^2-a*C)*(-sqrt((a-C)^2+4*b^2)-(a+C))))
# #, x0=x0,y0=y0,semiaxis.a=semiaxis.a,semiaxis.b=semiaxis.b
# make output consistent with others
return(invisible(list(parA = as.vector(A), centroid=centroid)  ) )
}


# support functions for fitConicLMA follow
# Not intended for standalone use
Residuals.ellipse <- function(XY,parG)  {
#This is a bug-fixed copy of the function in package conicfit by Jose Gama
Center <- parG[1:2]
Axes <- parG[3:4]
Angle <- parG[5]
n <- dim(XY)[1]
XYproj <- matrix(0,n,2)
tolerance <- 1e-9
#  First handling the circle case  
if (abs((Axes[1]-Axes[2])/Axes[1])<tolerance){
    phiall <- atan2(XY[,1]-Center[1] , (XY[,2]-Center[2]))
    XYproj <- matrix(c(Axes[1] * cos(phiall)+Center[1], Axes[2] * sin(phiall)+Center[2]),n,2)
# this doesn't work because he just reduced XYproj to a 1x2 matrix from nx2 
# All because he put "1" instead of "n" into the matrix() arg for rows. 
   RSS <- norm(XY-XYproj,'F')^2
    return(list(RSS=RSS, XYproj=XYproj))
}
# Now dealing with non-circle ellipses
a <- Axes[1]
b <- Axes[2]
aa <- a^2
bb <- b^2
tol.a <- tolerance * a
tol.b <- tolerance * b
tol.aa <- tolerance * aa
#  Matrix Q for rotating the points and the ellipse to the canonical system
s <- sin(Angle)
cA <- cos(Angle)
Qmat <- matrix(c(cA, -s,s, cA),2,2,byrow=TRUE)
#  data points in canonical coordinates
XY0  <- cbind(XY[,1]-Center[1], XY[,2]-Center[2]) %*% Qmat
XYA <- abs(XY0)
Tini <- apply(cbind(a * (XYA[,1]-a),b * (XYA[,2]-b)),1,max)
#  main loop over the data points
for (i in 1:n){
    u <- XYA[i,1]
 v <- XYA[i,2]
    ua <- u * a
     vb <- v * b
    if (u == 0) z1 <- 1 else z1 <- sign(XY0[i,1])
    if (v == 0) z2 <- 1 else z2 <- sign(XY0[i,2])
    #       does the point lie on the minor axis?
    if (u<tol.a){
        if (XY0[i,2]<0) XYproj[i,] <- cbind(0, -b) else XYproj[i,] <- cbind(0, b)
        next
    }
    #       does the point lie on the major axis?
    if (v<tol.b){
        if (u < (a-bb/a)){
            xproj <- aa * u/(aa-bb)
            XYproj[i,] <- matrix(c(z1 * xproj, z2 * b * sqrt(max(1-(xproj/a)^2,0))),1,2)
        } else {
            XYproj[i,] <- cbind(z1 * a, 0)
        }
        next
    }
    #      generic case: start the iterative procedure
    T <- Tini[i]
    for (iter in 1:100){
        Taa <- T + aa
        Tbb <- T + bb
        PP1 <- (ua/Taa)^2
        PP2 <- (vb/Tbb)^2
        Fp  <- PP1 + PP2 - 1
        if (Fp<0) break
        Fder <- 2 * (PP1/Taa + PP2/Tbb)
        Ratio <- Fp/Fder
        if (Ratio<tol.aa) break
        T <- T + Ratio
    }
    #       compute the projection of the point onto the ellipse
    xproj <- XY0[i,1] * aa/Taa
    yproj <- sign(XY0[i,2]) * b * sqrt(max(1-(xproj/a)^2,0))
    XYproj[i,] <- cbind(xproj, yproj)
} # end the main loop
#    rotate back to the original system
XYproj <- XYproj %*% t(Qmat)
XYproj <- cbind(XYproj[,1]+Center[1], XYproj[,2]+Center[2])
RSS <- norm(XY-XYproj,'F')^2
# message('leaving residuals.ellipse')
list(RSS=RSS,XYproj=XYproj)
}   # Residuals.ellipse

# Not intended for standalone use
Residuals.hyperbola <- function(XY,parG) {
#This is a copy of the function in package conicfit by Jose Gama,
# with some code added to avoid crashes when Inf values are generated
#   Projecting a given set of points onto a hyperbola and computing the distances from the points to the hyperbola
Center <- parG[1:2]   
Axes <- parG[3:4]  
Angle <- parG[5]
n <- dim(XY)[1]
XYproj <- matrix(0,n,2)
tolerance <- 1e-9
a <- Axes[1]
b=Axes[2]
aa <- a^2  
bb <- b^2 
at=sqrt(a)
bt=sqrt(b)
tol.a <- tolerance*a
tol.b <- tolerance*b
tol.aa <- tolerance*aa
#  Matrix Qmat for rotating the points and the hyperbola
s <- sin(Angle)
cA <- cos(Angle)
Qmat <- matrix(c(cA, -s,s, cA),2,2,byrow=TRUE)
#  data points in canonical coordinates
XY0  <- (cbind(XY[,1]-Center[1], XY[,2]-Center[2])) %*% Qmat
XYA <- abs(XY0)
# BUG: get rid of Inf vals
gooda <- is.finite(XYA[,1])&is.finite(XYA[,2])
XYA <- XYA[gooda, ]
# find the inflection points TInf for all given pairs
XYS <- sqrt(XYA)
XYS1 <- matrix(bb*at*XYS[,1]-aa*bt*XYS[,2],ncol=1)
XYS2 <- matrix(at*XYS[,1]+bt*XYS[,2],ncol=1)
TInf <- matrix(XYS1/XYS2,ncol=1) #inflection points
# have to clean this out too
goodf <- is.finite(TInf) 
TInf <- TInf[goodf ]
#  main loop over the data points
for (i in 1:n){
    u <- XYA[i,1]  
v <- XYA[i,2]
    ua <- u*a      
vb <- v*b
    if (u == 0) z1 <- 1 else z1 <- sign(XY0[i,1])
    if (v == 0) z2 <- 1 else z2 <- sign(XY0[i,2])
    #       does the point lie on the major axis?
    if (v<tol.b){
        if (u>a+bb/a){
            xproj <- aa*u/(aa+bb)
            XYproj[i,] <- cbind(z1*xproj, z2*b*sqrt(max((xproj/a)^2-1,0)))
       } else XYproj[i,] <- cbind(z1*a, 0)
        next
    } # end if
    #       does the point lie on the minor axis?
    if (u<tol.a){
        yproj <- bb*v/(aa+bb)
        XYproj[i,] <- cbind(z1*a*sqrt(1+(yproj/b)^2), z2*yproj )
        next
    } # } if
    #     generic case: start the iterative procedure
    T0 <- TInf[i]  # inflection point t
    F0  <- (ua/(T0+aa))^2 - (vb/(-T0+bb))^2 - 1  # the corresponding F(t)
    # if F0>0 then pick the initial point to the right of the root
    # if F0<0 then pick the initial point to the left of the root
if ( is.null(F0) || !length(F0) || is.nan(F0) || !is.finite(F0)) browser()

    if (F0 > 0){  # case1: pick the initial point T until F(T)<0
        for (j in 1:20){
            Tp<-bb-(-T0+bb)/2^j
            Fp<-(ua/(Tp+aa))^2 - (vb/(-Tp+bb))^2 - 1
            if (Fp<0) break
        } # end for j
        for (iter in 1:100){
            Taa <- Tp + aa
            Tbb <- -Tp + bb
            PP1 <- (ua/Taa)^2
            PP2 <- (vb/Tbb)^2
            Fp  <- PP1 - PP2 - 1
            if (Fp>0) break
            Fder <- 2*(PP1/Taa + PP2/Tbb)
            Ratio <- Fp/Fder
            if (abs(Ratio)<tol.aa) break
            Tp <- Tp + Ratio
        } } else { for (j in 1:20) { # case2: pick the initial point T until F(T)>0
            Tp=-aa+(T0+aa)/2^j
            Fp=(ua/(Tp+aa))^2 - (vb/(-Tp+bb))^2 - 1
            if (Fp>0) break
        } # end for j
        for (iter in 1:100){
            Taa <- Tp + aa
            Tbb <- -Tp + bb
            PP1 <- (ua/Taa)^2
            PP2 <- (vb/Tbb)^2
            Fp  <- PP1 - PP2 - 1
            if (Fp<0)  break
            Fder <- 2*(PP1/Taa + PP2/Tbb)
            Ratio <- Fp/Fder
            if (abs(Ratio)<tol.aa) break
            Tp <- Tp + Ratio
        }
    } # end if
    #   compute the projection of the point onto the hyperbola
    if (Taa < 1e-6){
        yproj <- XY0[i,2]*bb/Tbb
        XYproj[i,] <- cbind(sign(XY0[i,1])*a*sqrt(1+(yproj/b)^2), yproj)
    } else { if (Tbb < 1e-6){
        xproj=XY0[i,1]*aa/Taa
        yproj <- sign(XY0[i,2])*b*sqrt(max((xproj/a)^2-1,0))
        XYproj[i,] <- cbind(xproj, yproj)
   } else XYproj[i,] <- cbind(XY0[i,1]*aa/Taa, XY0[i,2]*bb/Tbb)
    }  # end if
} # } the main for-loop
XYproj <- XYproj %*% t(Qmat)
XYproj <- cbind(XYproj[,1]+Center[1], XYproj[,2]+Center[2])
RSS <- norm(XY-XYproj,'F')^2
list(RSS=RSS, XYproj=XYproj)
}   # Residuals.hyperbola

# Not intended for standalone use
JmatrixLMA<-function(XY,parA,XYproj){
#This is a copy of the function in package conicfit by Jose Gama
#Compute the Jacobian matrix(Implicit method)using algebraic parameter)
# safety check
parA <- as.vector(parA)
n <- dim(XY)[1]
Res <- matrix(0,n,1)
X <- matrix(XY[,1],ncol=1)
Y <- matrix(XY[,2],ncol=1) 
Z <- cbind(X^2, X*Y, Y^2, X, Y, matrix(1,n,ncol=1) ) %*% parA
DD <- XY-XYproj
for (i in 1:n) Res[i]  <- sign(Z[i])*norm(matrix(DD[i,]),'F')
D2 <- matrix(c(XYproj,matrix(1,n,1)),ncol=3,byrow=FALSE)
x <- matrix(XYproj[,1],ncol=1)
y <- matrix(XYproj[,2],ncol=1)
xx <- matrix(x^2,ncol=1)
yy <- matrix(y^2,ncol=1)
xy <- matrix(x*y,ncol=1)
# dPar <- matrix(c(xx,xy,yy,x,y,ones(n,1)),1)       #partial derivative of P wrt parA
du <- D2 %*% matrix(c(2*parA[1],parA[2],parA[4]),3)   #partial derivative of P wrt x-coordinate
dv <- D2 %*% matrix(c(parA[2],2*parA[3],parA[5]),3)   #partial derivative of P wrt y-coordinate
eA  <- sqrt(du^2+dv^2)
J <- cbind(xx/eA, xy/eA, yy/eA, x/eA, y/eA,matrix(1,n,1)/eA)
list(Res = Res,J = J)
}
