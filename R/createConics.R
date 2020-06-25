#  Some info standards:
#  Conics can be classified (or distinguished) by computing the discriminant,
#   B2 -  4AC , of the general quadratic equation,
# Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
# If  B != 0 , then the conic is rotated such that its major (and minor) axis is no longer parallel to one of the coordinate axes.
# 1. If  B2 - 4AC > 0 , then the conic is a hyperbola. This is equivalent to showing that
# A and C have opposite signs.
# 2. If B2- 4AC = 0 , then the conic is a parabola. This is equivalent to showing that
# either A or C is equal to zero.
#  (because if B,A,C all nonzero we reduce to (x+ky)^2 = j, 
#  which  violates the assumption that A or C is zero)
# 3. If B2 - 4AC < 0 , then the conic is either an ellipse or a circle. This is equivalent to showing that A and C have the same sign. Circles are distinguished from  ellipses when 
# A = C.
# 
#  the "G" parameter set is  CenterX, CenterY, semimajorAxisX, smAxisY, theta .


# Generate any conic section from ABCDEF coefficients OR, as possible, from G-set
createConic <- function(x,param, conicType, ranFun=NULL, noise = 1, seedit = NULL, tol = 1e-6) {
#  quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0 
# for a given x0, use quadratic formula to find all possible y0.
#  cy^2 + y*(b*x0+e) + (a*x0^2 +d*x0 +f) =0
x <- x[is.numeric(x)] # just in case
if (length(param) == 5) {
	if( missing(conicType) || (!length(intersect(tolower(conicType[1]), c('e','h')))) ) stop('G-type params require conicType= "e" or "h" ')
	parA <- GtoA(param, conicType)$parA
	} else {
		parA <- param
		if (length(param) != 6 ) {
			stop( c('param must have 5 or 6 elements but has ',length(param)) )
			}
		}
#  {A,B,C}term are  for  Ay^2 + By + C = 0 form of parA
Aterm <- parA[3]
Bterm <- parA[2]*x + parA[5]
Cterm <- parA[1]*x^2 + parA[4]*x + parA[6]
# But in simple parabola, e.g. where c is zero (and b), it's simply  C we want. 
# split into separate cases for c zero and for a zero; 
if (abs(parA[3]) < tol) {
	# no y^2 term
	if (abs(parA[2]) < tol && abs(parA[5] < tol)) {
	# denom (e + bx) is always zero
		y <-  Cterm  
		outs <- cbind(x,y)
		} else {
		# denom (e+bx) will be zero at a couple singularities at worst
			y <- -Cterm/Bterm  # not +C/B 
			outs <- cbind(x,y)			
		}
} else {
# this else refers back to the first if
		theroot<-suppressWarnings( sqrt(Bterm^2-4*Aterm*Cterm) )
		y1 <- (-Bterm + theroot)/(2*Aterm)
		y2 <- (-Bterm - theroot)/(2*Aterm) 
# rev() puts data in nice order for line plotting
		y <- c(y1,rev(y2) )
		x <- c(x,rev(x)) # must match lengths
# removing NaNs
# hmmm - don't want to do this because it changes the size of the output
# and that might be bad in scripts which want to compare input to output
		# goods <- !is.nan(y) & !is.nan(x)
		# x <- x[goods]
		# y <- y[goods]
		outs <- cbind(x,y)
	}
outs <- unique(outs)  # in case of dupes in y1,y2 calculation
#  noise calculator
if (is.function(ranFun)) {
	if(!is.null(seedit)) set.seed(seedit)
	outs[,2] <- outs[,2] + ranFun(length(outs[,2]) ) * noise
	}
colnames(outs)<-c('x','y')
return(invisible(outs))
}

