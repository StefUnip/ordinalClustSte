bosclust <- function(x=matrix(0,nrow=1,ncol=1), idx_list=c(0), kr, init, nbSEM, nbSEMburn, nbindmini, m=0){
	checkParamsClust(x, init, nbSEM, nbSEMburn)
	res <- clust(xMat=x, myList=idx_list, kr, init, nbSEM, nbSEMburn, nbindmini, m=m)
	if(length(res@icl)==0){
		warning('Error: probably empty clusters')
	}
	return(res)
}