bosclassif.no.kc <-
  function (x,y,to.predict,kr,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
            init='kmeans',disp=TRUE,iterordiEM=10) {
    
    # setting progress bar
    if(disp) 
    {
      pb <- progress_bar$new(
        format="1/2 [:bar] :percent",
        total=nbSEM, clear=FALSE, width=60
      )
    }

   # constant for polynomial probability (BOS) 
    tab_pej=tabpej(m)
    
    # ---
    n=nrow(x)
    d=ncol(x)
    missing=FALSE
    if (sum(x==0)>0) {
      missing=TRUE
      # aleatory imputation for missing values
      miss=which(x==0)
      x[miss]=sample(1:m,sum(x==0),replace=TRUE)
    }

##### if missing values, an SEM gibbs is required #####

if(missing){
  V=array(0,c(n,kr))       
  gamma=rep(0, kr)
  mu=array(0,c(kr,d,nbSEM+1))
  p=array(0,c(kr,d,nbSEM+1))
  res_mu=array(0,c(kr,d))
  res_p=array(0,c(kr,d))
  res_gamma=rep(0, kr)
  res_V=array(0,c(n,kr))  
  res_zr=y


  # estimate the gammas 
  for(k in 1:kr){
    gamma[k] <- length(which(y==k))/length(y)
    res_gamma[k] <- length(which(y==k))/length(y)
  }

  #creating the V (we know them)
  for(i in 1:n){
    V[i,y[i]] = 1
    res_V[i,y[i]] = 1
  }


    # ==== init ==== 
    #No need to initialize the V and W because we know them

    # ---- parameters initialization from partitions ----
    for (h in 1:d){
      for (k in 1:kr){
        res <- ordiemCpp(m,tab_pej,as.vector(x[which(V[,k]==1),h]),
                         tabmu0=1:m,tabp0=seq(0,1,0.2),
                         iter_max=iterordiEM)
        mu[k,h,1]=res[[1]]
        if(res[[2]]==1) {
          p[k,h,1] = 0.999
        }
        else{
          p[k,h,1]=res[[2]]
        }
      }
    }
    # ---- missing values initialization ----
    if (missing){
      x[miss]=0 
      for (h in 1:d){
        for (k in 1:kr){

          tmp=which(x[which(V[,k]==1),h]==0)

          if (length(tmp)>0){
            probaBOS=rep(0,m)
            for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,mu[k,h,1],]*(rep(p[k,h,1],m)^(0:(m-1)))))
            x[which(V[,k]==1),h][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
          }
        }
      }
    }
    # ============  SEM ============
    for (iter in 1:nbSEM){
      if (disp) pb$tick()
      # ==== SE step ==== 
      # --- no need to compute proba for column and row partition because we know it
      # we just need this part for the missing values
      
      
      # --- imputation des donnees manquantes ----
      #if (missing){
        x[miss]=0 
        for (h in 1:d){
          for (k in 1:kr){
            
            tmp=which(x[which(V[,k]==1),h]==0)
            
            if (length(tmp)>0){
              probaBOS=rep(0,m)
              for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,mu[k,h,iter],]*(rep(p[k,h,iter],m)^(0:(m-1)))))
              x[which(V[,k]==1),h][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
            }
          }
        }
      #}
      
      # ==== M step ==== 
      
        for (h in 1:d){
          for (k in 1:kr){
            xtmp=x
            xtmp[miss]=0
            tmp=as.vector(x[which(V[,k]==1),h])
            datablock_kl=tmp[tmp>0]
            res <- ordiemCpp(m,tab_pej,datablock_kl,
                      tabmu0=1:m,tabp0=p[k,h,iter],
                      iter_max=iterordiEM)
            mu[k,h,iter+1]=res[[1]]
            if(res[[2]]==1){
              p[k,h,iter+1] = 0.999  
            }
            else{
              p[k,h,iter+1]=res[[2]]  
            }
          }
        }
      
    }# for iter
    # ===== parameters computaton (mode and median after burn-in) =====
    for (h in 1:d){
      for (k in 1:kr){
        res_mu[k,h]=mode(mu[k,h,nbSEMburn:(nbSEM+1)])
        res_p[k,h]=median(p[k,h,nbSEMburn:(nbSEM+1)])
      }
    }
}
else{ #### case no missing value: NO SEM neeeded

  V=array(0,c(n,kr))       
  gamma=rep(0, kr)
  mu=array(0,c(kr,d))
  p=array(0,c(kr,d))
  res_mu=array(0,c(kr,d))
  res_p=array(0,c(kr,d))
  res_gamma=rep(0, kr)
  res_V=array(0,c(n,kr))  
  res_zr=y

  # estimate the gammas 
  for(k in 1:kr){
    gamma[k] <- length(which(y==k))/length(y)
    res_gamma[k] <- length(which(y==k))/length(y)
  }

  #creating the V (we know them)
  for(i in 1:n){
    V[i,y[i]] = 1
    res_V[i,y[i]] = 1
  }


  for (h in 1:d){
    for (k in 1:kr){

      res <- ordiemCpp(m,tab_pej,as.vector(x[which(V[,k]==1),h]),
                           tabmu0=1:m,tabp0=seq(0,1,0.1),
                           iter_max=iterordiEM)
      mu[k,h]=res[[1]]
      if(res[[2]]==1) {
        p[k,h] = 0.999
      }
      else{
        p[k,h]=res[[2]]
      }
    }
  }

  res_mu=mu
  res_p=p
}
    
    
     # ===== estimation of partitions  =====
      if(disp) 
      {
        pb2 <- progress_bar$new(
          format="2/2 [:bar] :percent",
          total=nbSEM, clear=FALSE, width=60
        )
      }
      Q=nbSEM
      Xhat=array(0,c(n,d,Q+1))
      Xhat[,,1:(Q+1)]=x
      for (iterQ in 1:Q){
        # --- no need for simulation about column  or row partitions
        if(disp) pb2$tick()
       # --- missing values simulation ---
        if (missing){
          tmpx=x
          tmpx[miss]=0 # on remet des 0 la ou il y avait des data manquantes
          for (h in 1:d){
            for (k in 1:kr){
              # recherche des cases manquantes
              tmp=which(tmpx[which(V[,k]==1),h]==0)
              # simulation des data manquantes
              if (length(tmp)>0){
                probaBOS=rep(0,m)
                for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,res_mu[k,h],]*(rep(res_p[k,h],m)^(0:(m-1)))))
                tmpx[which(V[,k]==1),h][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
              }
            }
          }
          Xhat[,,iterQ+1]=tmpx
        }
      }#iterQ
      # --- final partition estimation  ---

       # --- missing vlues final estimtion ---
      for (i in 1:n){
        for (h in 1:d){
          if (x[i,h]==0) x[i,h]=mode(Xhat[i,h,])
        }
      }
      
      xhat=x
     


  
    zr=res_zr



    # Now that every thing is estimated, we can estimate the predictions for to.predict

    # First, impute missing values of validation dataset:
    missing.validation = FALSE
    if (sum(to.predict==0)>0) {
      missing.validation=TRUE
      miss.val=which(to.predict==0)
    }
    if(missing.validation){

        #for (h in 1:d){
        #  for (k in 1:kr){
        #    # recherche des cases manquantes
        #    tmp=which(to.predict[which(V[,k]==1),h]==0)
        #    # simulation des data manquantes
        #    if (length(tmp)>0){
        #      probaBOS=rep(0,m)
        #      for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,mu[k,h],]*(rep(p[k,h],m)^(0:(m-1)))))
        #      to.predict[which(V[,k]==1),h][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
        #    }
        #  }
        #}
        to.predict[miss.val]=sample(1:m,sum(to.predict==0),replace=TRUE)
    }

    #Then predict:
    n.to.predict= nrow(to.predict)
    logprobas <- matrix(0,nrow=n.to.predict,ncol=kr)
    probas <- matrix(0,nrow=n.to.predict,ncol=kr)

    logprobas=matrix(log(gamma),nrow=n.to.predict,ncol=kr,byrow = T)

    for (k in 1:kr){
      for(h in 1:d){
        for (i in 1:n.to.predict){ 

          #logprobas[i,k] = logprobas[i,k] + 
          #        log(sum(tab_pej[to.predict[i,h], res_mu[k,h],] *
          #                  (rep(res_p[k,h],m)^(0:(m-1)))))
          
          sum.neg <- sum(tab_pej[to.predict[i,h], res_mu[k,h],] *
                            (rep(res_p[k,h],m)^(0:(m-1))))
          
          
          logprobas[i,k] = logprobas[i,k] + log(sum.neg)
              
   
        }
      }
    }
    


    V.to.predict <- array(0,c(n.to.predict,kr))
    for (i in 1:n.to.predict){
      tmp=logsum (logprobas[i,])
      probas[i,]=exp ( logprobas[i,] - tmp )
      #V.to.predict[i,]=rmultinom(1,1,probas[i,])
      idx = which(probas[i,]==max(probas[i,]))
      V.to.predict[i,idx] = 1
    }

    
    zr.to.predict=apply(apply(V.to.predict,c(1,2),sum),1,which.max)



    # --- return the result ---
    res=list(probas.to.predict=probas,
             zr.to.predict=zr.to.predict,
             V.to.predict=V.to.predict,
             xhat=xhat,
             res_mus=res_mu,
             res_ps=res_p,
             zr=zr,
             mus=mu,
             ps=p,
             gamma=gamma,
             V=V,
             string="classif.no.kc",
             m=m)
    
    return(res) 
  }
