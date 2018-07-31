predictions <- function(classif, x=matrix(0,nrow=1,ncol=1)){
	res <- prediction(classif, x)
	return(res)
}