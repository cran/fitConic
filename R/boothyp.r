# hyperbola bootstrappers, i.e. cheap guesses without initial parameters
# usually locks to theta = pi/2, but output is sufficient for main algorithm
# to converge.
#   maxiter to be passed to optim()
bootHyperbola <-function(x,y= NULL,maxiter = 1e4,...) {
xy <-xy.coords(x,y)
hdat <- cbind(xy$x,xy$y)[order(xy$x),]
# costhyp will compare newy rotated back into xy space
#initialize parameters
b3 = c(1,1,1)
Ang = 0
#wrapper for cost function to feed to optim()
fhypopt <-function(prm=c(b3 = c(b3[1],b3[2],b3[3]),Ang), xy=xy){ fhyp(xy,prm[1:3],prm[4])
}
thefit <- optim(c(b3=b3,Ang= -Ang),fhypopt, xy = hdat, control=list(maxit=maxiter))
# convert to ABCDEF: y = a0 + a1/(x+a2) -->  xy - a0*x + a2*y -(a2*a0 +a1) = 0 
newb <- thefit$par[1:3]   
parA <- c(0, 1, 0, -newb[1], newb[3], -(newb[1] * newb[3] + newb[2]) )
# rotate for final coords to match
# need a neg angle here???
parAr <-rotateA(parA, thefit$par['Ang'])$parA
return(invisible(list(parA=parA, parAr=parAr, theta = thefit$par['Ang'], fitdat = thefit)) )
}
 
# need a 'catcher' for residcost going infinite. Return a big number to 
# make optim() happy
# Not intended for external use
 fhyp <- function (xy, b3, Ang) {
xy <- xy[is.finite(xy[,1]) & is.finite(xy[,2]),]
#play with sign of Ang here too .Orig was positive
xyr <-xyrot(xy, theta = Ang)
newy <- b3[1] + b3[2]/(b3[3] + xyr[,1])
# orig was negative
newyrot <-xyrot(xyr[,1], newy, theta = -Ang)[,2]
residcost <- norm(as.matrix(xy[,2]-newyrot),'F')
if (!is.finite(residcost)) residcost <- 1e5 
return(residcost)
}




