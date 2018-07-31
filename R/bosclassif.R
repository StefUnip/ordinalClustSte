bosclassif <- function(x=matrix(0,nrow=1,ncol=1), y, idx_list=c(0), kr, kc=0, init, nbSEM, nbSEMburn, nbindmini, m=0){
	checkParamsClassif(x, init, nbSEM, nbSEMburn)
	if(kc[1]!=0){
		res <- classif(xMat=x, y, myList=idx_list, kr, kc, init, nbSEM, nbSEMburn, 
			nbindmini, m=m)
	}
	else{
		res <- classifM(xMat=x, y, myList=idx_list, kr, init, nbSEM, nbSEMburn, 
			nbindmini, m=m)
	}
	if(length(res@icl)==0){
		warning('Error: probably empty clusters')
	}
	return(res)
}