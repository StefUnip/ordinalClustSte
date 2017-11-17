bosclassif.kc <-
  function (x,y,to.predict,kr,kc=1,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
            init='kmeans',disp=TRUE,iterordiEM=10) {
    # ----------------------------------------------------------------------------
    # Estimation of the latent BOS coclustering model via SEM algoritm
    # input
    #   x  : matrice n x p de donnees ordinales (individu en ligne, variable en colonne)
    #   kc : nb de classes en colonne
    #   kr : nb de classes en ligne  
    #   m  : nombre de modalites (identique pour toutes les variables) 
    #   nbSEM : nombre d'iterations de l'algo. SEM
    #   nbburn : taille de la periode de burn
    # ----------------------------------------------------------------------------
    
    # setting progress bar
    if(disp) 
    {
      pb <- progress_bar$new(
        format="1/2 [:bar] :percent",
        total=nbSEM, clear=FALSE, width=60
      )
    }

    # on charge en memoire les exposants des probas BOS sous forme polynomiales
    tab_pej=tabpej(m)
    

    # ---
    n=nrow(x)
    d=ncol(x)


    missing=FALSE
    if (sum(x==0)>0) {
      missing=TRUE
      # imputation aleatoire lors de l'init
      miss=which(x==0)
      x[miss]=sample(1:m,sum(x==0),replace=TRUE)
    }

  
    V=array(0,c(n,kr))
    W=array(0,c(d,kc,nbSEM+1))        
    gamma=rep(0, kr)
    rho=array(0,c(kc,nbSEM+1))
    mu=array(0,c(kr,kc,nbSEM+1))
    p=array(0,c(kr,kc,nbSEM+1))
    res_mu=array(0,c(kr,kc))
    res_p=array(0,c(kr,kc))
    res_gamma=rep(0, kr)
    res_rho=matrix(0,kc)
    res_V=array(0,c(n,kr)) 
    res_W=array(0,c(d,kc)) 
    res_zr=y
    res_zc=matrix(0,d)


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
      # --- init aleatoire des partitions ----
      if (init=='random'){
        W[,,1]=t(rmultinom(d,1,rep(1/kc,kc)))
        while (! verif(x,V,W[,,1],kc,nbindmini))
        {
          print('reload random init')
          W[,,1]=t(rmultinom(d,1,rep(1/kc,kc)))
        }
      }
      
      # --- init kmeans des partitions ----
      if (init=='kmeans'){
        tmpW=kmeans(t(x),kc,nstart=10)  
        for (i in 1:d) W[i,tmpW$cluster[i],1]=1
        while (! verif(x,V,W[,,1],kc,nbindmini))
        {
          print('reload kmeans init')
          W[,,1]=0
          tmpW=kmeans(t(x),kc,nstart=2)  
          for (i in 1:d) W[i,tmpW[i]$cluster,1]=1
        }  
      }
      # init des parametres a partir des partitions
      rho[,1]=getMeans(W[,,1])
      for (l in 1:kc){
        for (k in 1:kr){
          # init des parametres en fonction des partitions aleatoires
          # res=ordiem(as.vector(x[which(V[,k,1]==1),which(W[,l,1]==1)]),m,tabmu0=1:m,tabp0=seq(0,1,0.2),iter_max=iterordiEM)
          res <- ordiemCpp(m,tab_pej,as.vector(x[which(V[,k]==1),
                                                 which(W[,l,1]==1)]),
                           tabmu0=1:m,tabp0=seq(0,1,0.2),
                           iter_max=iterordiEM)
          mu[k,l,1]=res[[1]]
          if(res[[2]]==1) {
            p[k,l,1] = 0.999
          }
          else{
            p[k,l,1]=res[[2]]
          }
        }
      }
      # === init des valeurs manquantes ===
      if (missing){
        x[miss]=0 # on remet des 0 la ou il y avait des data manquantes
        for (l in 1:kc){
          for (k in 1:kr){
            # recherche des cases manquantes
            tmp=which(x[which(V[,k]==1),which(W[,l,1]==1)]==0)
            # simulation des data manquantes
            if (length(tmp)>0){
              probaBOS=rep(0,m)
              for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,mu[k,l,1],]*(rep(p[k,l,1],m)^(0:(m-1)))))
              x[which(V[,k]==1),which(W[,l,1]==1)][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
            }
          }
        }
      }



        # === debut du SEM ===
        for (iter in 1:nbSEM){
          #if (disp) cat('iteration ',iter,'/',nbSEM,'\n')
          if(disp) pb$tick()
          # ==== SE step ==== 
          # --- calcul des proba pour la simulation de la partition en colonne

          logprobaW=matrix(log(rho[,iter]),nrow=d,ncol=kc,byrow = T)
          for (k in 1:kr){
            for(l in 1:kc){
              for (h in 1:d){ 
                for (tmp in 1:m){
                  logprobaW[h,l]=logprobaW[h,l] + sum((x[,h]*V[,k])==tmp) * log(sum(tab_pej[tmp,mu[k,l,iter],]*(rep(p[k,l,iter],m)^(0:(m-1)))))     
                }
              }
            }
          }
          
          
          if (iter>1) probaWold=probaW
          probaW=matrix(0,d,kc)
          for (h in 1:d){
            tmp=logsum (logprobaW[h,])
            probaW[h,]=exp ( logprobaW[h,] -  tmp)
          }
          
          if (iter>1) {
            maxW=max(abs(probaWold-probaW))
          }
          casevide=TRUE
          restart=0
          while (casevide && (restart<1000)){
            for (h in 1:d){
              W[h,,iter+1]=rmultinom(1,1,probaW[h,])
            }
            
            casevide = (! verif(x,V,W[,,iter+1],kc,nbindmini))
            if (casevide){
              restart=restart+1
              #cat('restart number ',restart,'\n')
            }
          }
          # --- imputation des donnees manquantes ----
          if (missing){
            x[miss]=0 # on remet des 0 la ou il y avait des data manquantes
            for (l in 1:kc){
              for (k in 1:kr){
                # recherche des cases manquantes
                tmp=which(x[which(V[,k]==1),which(W[,l,iter+1]==1)]==0)
                # simulation des data manquantes
                if (length(tmp)>0){
                  probaBOS=rep(0,m)
                  for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,mu[k,l,iter],]*(rep(p[k,l,iter],m)^(0:(m-1)))))
                  x[which(V[,k]==1),which(W[,l,iter+1]==1)][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
                }
              }
            }
          }
          if ((restart==1000) && casevide){
            print('The algorithm is stopped for degenerancy reason')
            return(NULL)
          }
          # ==== M step ==== 
          
            rho[,iter+1]=getMeans(W[,,iter+1])
            for (l in 1:kc){
              for (k in 1:kr){
                tmp=as.vector(x[which(V[,k]==1),which(W[,l,iter+1]==1)])
                datablock_kl=tmp[tmp>0]
                # res=ordiem(datablock_kl,m,tabmu0=1:m,tabp0=p[k,l,iter],iter_max=iterordiEM)
                res <- ordiemCpp(m,tab_pej,datablock_kl,
                          tabmu0=1:m,tabp0=p[k,l,iter],
                          iter_max=iterordiEM)
                mu[k,l,iter+1]=res[[1]]
                if(res[[2]]==1) {
                  p[k,l,iter+1] = 0.999
                }
                else{
                  p[k,l,iter+1]=res[[2]]
                }
              }
            }
          
        }# for iter
        # ===== calcul des parametres (mode et median hors burn) =====
        for (l in 1:kc){
          res_rho[l]=median(rho[l,nbSEMburn:(nbSEM+1)])
          for (k in 1:kr){
            res_mu[k,l]=mode(mu[k,l,nbSEMburn:(nbSEM+1)])
            res_p[k,l]=median(p[k,l,nbSEMburn:(nbSEM+1)])
          }
        }
        res_rho=res_rho/sum(res_rho)



      # ===== estimation des partitions et des valeurs manquantes  =====
      if(disp) 
      {
        pb2 <- progress_bar$new(
          format="2/2 [:bar] :percent",
          total=nbSEM, clear=FALSE, width=60
        )
      }
      Q=nbSEM
      Xhat=array(0,c(n,d,Q+1))
      Wfinal=array(0,c(d,kc,Q+1)) 
      Wfinal[,,1]=W[,,iter+1]
      Xhat[,,1:(Q+1)]=x
      for (iterQ in 1:Q){
        if(disp) pb2$tick()
        # --- simulation des partitions en ligne et en colonne ---
        logprobaW=matrix(log(res_rho),nrow=d,ncol=kc,byrow = T)
        for (k in 1:kr){
          for(l in 1:kc){
            for (h in 1:d){ 
              for (tmp in 1:m){
                logprobaW[h,l]=logprobaW[h,l] + sum((Xhat[,h,iterQ]* V[,k])==tmp) * log(sum(tab_pej[tmp,res_mu[k,l],]*(rep(res_p[k,l],m)^(0:(m-1)))))     
              }
            }
          }
        }
        probaW=matrix(0,d,kc)
        for (h in 1:d){
          tmp=logsum (logprobaW[h,])
          probaW[h,]=exp ( logprobaW[h,] -  tmp)
          Wfinal[h,,iterQ+1]=rmultinom(1,1,probaW[h,])
        }
        # --- simulation des donnees manquantes ---
        if (missing){
          tmpx=x
          tmpx[miss]=0 # on remet des 0 la ou il y avait des data manquantes
          for (l in 1:kc){
            for (k in 1:kr){
              # recherche des cases manquantes
              tmp=which(tmpx[which(V[,k]==1),which(Wfinal[,l,iterQ]==1)]==0)
              # simulation des data manquantes
              if (length(tmp)>0){
                probaBOS=rep(0,m)
                for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,res_mu[k,l],]*(rep(res_p[k,l],m)^(0:(m-1)))))
                tmpx[which(V[,k]==1),which(Wfinal[,l,iterQ]==1)][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
              }
            }
          }
          Xhat[,,iterQ+1]=tmpx
        }
      }#iterQ
      # --- estimation de la partition finale par mode marginal ---

      res_zc=apply(apply(Wfinal,c(1,2),sum),1,which.max)

      for (h in 1:d) res_W[h,res_zc[h]]=1
      # --- estimation des valeurs manquantes par le mode marginal ---
      for (i in 1:n){
        for (h in 1:d){
          if (x[i,h]==0) x[i,h]=mode(Xhat[i,h,])
        }
      }
      # --- sauvegarde de la matrice completee ---
      xhat=x
      # --- approximation (plus faite)


    
    # estimation des partitions
    zr=res_zr
    zc=res_zc
    # calcul ICL Brault
    # formule adaptee a notre cas ( 1 au lieu de m-1 param par case), mais a verifier theoriquement
    if (!missing){
      icl=- (kr-1)/2 *log(n)- (kc-1)/2 *log(d)- kc*kr/2 *log(n*d)
      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log (res_gamma[k])
        }}
      for (h in 1:d){
        for (l in 1:kc){
          icl = icl + res_W[h,l] * log (res_rho[l])
        }}
      for (h in 1:d){
        for (l in 1:kc){
          for(i in 1:n){
            for(k in 1:kr){
              icl = icl + res_W[h,l] * res_V[i,k] * log ( (sum(tab_pej[x[i,h],res_mu[k,l],]*(rep(res_p[k,l],m)^(0:(m-1))))) )
            }}}}
    }
    else{
      
      # --  fin approximation par moyenne harmonique --
      icl=- (kr-1)/2 *log(n)- (kc-1)/2 *log(d)- kc*kr/2 *log(n*d)
      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log (res_gamma[k])
        }}
      for (h in 1:d){
        for (l in 1:kc){
          icl = icl + res_W[h,l] * log (res_rho[l])
        }}
      x[miss]=0
      for (h in 1:d){
        for (l in 1:kc){
          for(i in 1:n){
            for(k in 1:kr){
              if (x[i,h]>0){
                icl = icl + res_W[h,l] * res_V[i,k] * log ( (sum(tab_pej[x[i,h],res_mu[k,l],]*(rep(res_p[k,l],m)^(0:(m-1))))) )
              }
            }
          }
        }
      }
    }#if(!missing)

    # Now that every thing is estimated, we estimate the predictions for to.predict
    # First, impute missing values of validation dataset:
    missing.validation = FALSE
    if (sum(to.predict==0)>0) {
      missing.validation=TRUE
      miss.val=which(to.predict==0)
    }
    if(missing.validation){

        #for (l in 1:kc){
        #  for (k in 1:kr){
        #    # recherche des cases manquantes
        #    tmp=which(to.predict[which(V[,k]==1),l]==0)
        #    # simulation des data manquantes
        #    if (length(tmp)>0){
        #      probaBOS=rep(0,m)
        #      for (im in 1:m) probaBOS[im]=(sum(tab_pej[im,mu[k,l],]*(rep(p[k,l],m)^(0:(m-1)))))
        #      to.predict[which(V[,k]==1),l][tmp]=sample(1:m,length(tmp),prob=probaBOS,replace=TRUE) 
        #    }
        #  }
        #}
        to.predict[miss.val]=sample(1:m,sum(to.predict==0),replace=TRUE)
    }

    n.to.predict= nrow(to.predict)
    logprobas <- matrix(0,nrow=n.to.predict,ncol=kr)
    probas <- matrix(0,nrow=n.to.predict,ncol=kr)

    logprobas=matrix(log(gamma),nrow=n.to.predict,ncol=kr,byrow = T)
    for (k in 1:kr){
      for(l in 1:kc){
        for (i in 1:n.to.predict){ 
         for (tmp in 1:m){
           logprobas[i,k]=logprobas[i,k] +  
             sum((to.predict[i,]*res_W[,l])==tmp)* 
             log(sum(tab_pej[tmp,res_mu[k,l],]*
             (rep(res_p[k,l],m)^(0:(m-1)))))     
         }    
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
             mus=res_mu,
             ps=res_p,
             res_rho=res_rho,
             res_W=res_W,
             icl=icl,
             zr=zr,
             zc=zc,
             mus=mu,
             ps=p,
             gamma=gamma,
             rho=rho,
             V=V,
             W=W,
             probaW=probaW,
             string="classif.kc",
             m=m)

    return(res) 
  }
