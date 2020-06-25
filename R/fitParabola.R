#based loosely on:
# http://www.mathworks.com/matlabcentral/answers/80541
#   derotate by finding a good theta, so that the data are now 1:1 with
# respect to x values, then a simple polyfit will do.  
# note that vertex translation in x is implicit in the existence of a linear x-term
#
# Added optional angle-range values to reduce time spent ransac-ing. Use with care,
# as optimizer may 'clamp' at one extreme. 
# in this setup, the minimum is in a very narrow window and there are other local minima
#

fitParabola <-function(x,y=NULL,searchAngle=c(-pi/2, pi/2),... ){
xy <- xy.coords(x,y)
xy<-cbind(xy$x,xy$y)
# optimize theta, starting at zero, with xy data passed into cost functions
# Seems happier with optimize() than optim()
#  bar2 <-optim(0, costparabxy, gr=NULL, xy, method='Brent',lower = searchAngle[1], upper = searchAngle[2])
bar2 <- optimize(costparabxy, searchAngle,xy=xy)
theta <-  bar2$minimum  # bar2$par  for optimize()
finalcost <- costparab(theta,xy) 
coeffs <-finalcost$coeffs  
 
#raw x- vertex, prior to rotation back to source data position
xv <- - coeffs[2]/2/coeffs[3]  
# y = a(x-h)^2 + k, vertex (h,k), which means (-b/2a , plug that into form)
vertex <- xyrot( xv, coeffs[1]+coeffs[2]*xv+coeffs[3]*xv^2, theta)
#return as ABCDEF plus a theta
# parA0 = c(coeffs[3],0,0,coeffs[2],0,coeffs[1])
# then rotate that via theta
costhet = cos(theta)
sinthet = -sin(theta) # watch for neg sign...
parA = coeffs[3]*costhet^2
parA[2] = 2*coeffs[3] * costhet * sinthet
parA[3] = coeffs[3] *sinthet^2
parA[4] = sinthet + coeffs[2]*costhet
parA[5] = coeffs[2]*sinthet - costhet
parA[6] = coeffs[1]
return(list(vertex=vertex, theta=theta, parA=parA, parQ = coeffs,  cost = finalcost$thecost) ) 

}
#end of main function 

# cost function. Not intended for separate use
costparab <- function(theta,xy){
# xyrot is in "coordconverters" in this package.
	rxy <-xyrot(xy, theta = -theta)
	lmout <- lm(rxy[,2] ~ I(rxy[,1]) + I(rxy[,1]^2) )	
# this is a reasonable "cost" value to optimize 
	normres <- norm(as.matrix(lmout$residuals),'F')
	return(list(thecost=normres, coeffs=lmout$coefficients))
}
#wrapper which returns one value to fit optimize() rules
costparabxy <- function(theta,xy) costparab(theta,xy)$thecost





