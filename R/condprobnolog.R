condprobnolog <-
function(f,sumf){
  # ----------------------------------------------------------------------------
  # Conditional probabilities of belonging to clusters (version without log).
  # input:
  #   f       components density [n,k]
  #	  sumf 	  mixture density [n]
  # output:
  # 	t	conditional probabilities [n,k]
  # ----------------------------------------------------------------------------
  n = nrow(f)
  k = ncol(f)  
  t = f / (rep(sumf,k,bycol=TRUE)) # t [n,k]
  return(t)
}
