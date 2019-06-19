boscoclust <- function(x=matrix(0,nrow=1,ncol=1), idx_list=c(1), 
	kr, kc, init, nbSEM, nbSEMburn, 
	nbRepeat=1, nbindmini, m=0, percentRandomB=0){
	idx_list <- idx_list - 1 # patch for indexes
	checkParamsCoclust(x, init, nbSEM, nbSEMburn)
	seed = get_seed()
	res <- coclust(xMat=x, myList=idx_list, kr, kc, init, nbSEM, nbSEMburn, 
		nbRepeat, nbindmini, m=m, percentRandomB=percentRandomB, seed=seed)
	if(length(res@icl)==0){
		warning('Error: probably empty clusters')
	}
	return(res)
}