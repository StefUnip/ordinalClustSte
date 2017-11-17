logsum <-
function (logx){
  # astuce pour le calcul des log de sommme
  #Log(a+b+c) = log(a) + log(1 + exp(log(b)-log(a)) + exp(log(c)-log(a)))
  if(length(logx)==1){
    return(logx)
  }
  logx = sort(logx,decreasing=TRUE)
  tmp=1
  for (i in 2:length(logx)) tmp=tmp+exp(logx[i]-logx[1])
  return(logx[1]+log(tmp))
}
