boscoclust <- function(x=matrix(0,nrow=1,ncol=1), idx_list=c(0), kr, kc, init, nbSEM, nbSEMburn, nbRepeat=1, nbindmini, m=0){
	checkParamsCoclust(x, init, nbSEM, nbSEMburn)
	res <- coclust(xMat=x, myList=idx_list, kr, kc, init, nbSEM, nbSEMburn, nbRepeat, nbindmini, m=m)
	if(length(res@icl)==0){
		warning('Error: probably empty clusters')
	}
	return(res)
}