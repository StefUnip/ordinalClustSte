predictions <- function(classif, x=matrix(0,nrow=1,ncol=1)){
	seed = get_seed()
	res <- prediction(classif, x, seed)
	return(res)
}