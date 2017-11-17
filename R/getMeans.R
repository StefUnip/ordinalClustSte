getMeans <-
function(M){
  if(is.vector(M)){
    return(mean(M))
  }
  else{
    return(colMeans(M))
  }
}
