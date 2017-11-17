verif <-
function(x,V,W,nb.col.cluster,nbindmini){
  # fonction qui verifie le nombre d'observation differentes minimale par block
  check=1 # verifier si dans chaque case il y a au moins nbindmini observations differentes
  
  
  
  if(nb.col.cluster!=1){
    if(is.vector(V)){
      if(is.vector(W)){
        check=0
      }
      else{
        check=0;
      }
    }
    else{
      for (ir in 1:ncol(V)){
        if (check==0) break
        if(is.vector(W)){
          check=0
        }
        else{
          for (ic in 1:ncol(W)){
            if (length(unique(as.vector(x[V[,ir]==1,W[,ic]==1])))<nbindmini){
              check=0;
              break
            } 
          }
        }
        
      }
    }
  }
  else{
    for (ir in 1:ncol(V)){
        
            if (length(unique(as.vector(x[V[,ir]==1,])))<nbindmini){
              check=0;
              break
            } 
        
      }
    
  }
  
  
  
  return(check)
}
