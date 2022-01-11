# doWeights helper function
# Normally only used internally but can be used to modify a dataset in general
#Note that XY must be a N X 2 array 
doWeights <- function(XY, weights) {

	if(length(weights) != length(XY[,1])){
		stop('weights vs. data length mismatch')
	}
	weights = round(weights) # just in case
	if(sum(weights < 0 )) {
		message('negative weights not allowed, setting to zero')
		weights[weights < 0 ] <- 0 
	}
	#now extend XY
	fakerle = list(values = XY[, 1], lengths = weights)
	xy1tmp <- inverse.rle(fakerle)
	fakerle = list(values = XY[, 2], lengths = weights)
	xy2tmp <- inverse.rle(list(values = XY[,2], lengths=weights))
	XY <- cbind(xy1tmp,xy2tmp)
	return(invisible(XY))
}