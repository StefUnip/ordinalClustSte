bosclust <- function(x=matrix(0,nrow=1,ncol=1), idx_list=c(1), kr, init, 
					nbSEM, nbSEMburn, nbindmini, m=0, percentRandomB=0){
	idx_list <- idx_list - 1 # patch for indexes
	checkParamsClust(x, init, nbSEM, nbSEMburn)
	seed = get_seed()
	res <- clust(xMat=x, myList=idx_list, kr, init, nbSEM, nbSEMburn, nbindmini, 
				m=m, percentRandomB=percentRandomB, seed=seed)
	if(length(res@icl)==0){
		warning('Error: probably empty clusters')
	}
	return(res)
}