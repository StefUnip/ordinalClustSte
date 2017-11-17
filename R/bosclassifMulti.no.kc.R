bosclassifMulti.no.kc <-
  function (x,y,to.predict,d.list,kr,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
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

    D<- length(m)

    # on charge en memoire les exposants des probas BOS sous forme polynomiales
    tab_pejs <- list()
    for(id in 1:D){
      tab_pejs[[id]] <- tabpej(m[id])
    }
    
    # ---
    n=nrow(x)
    nd <- c()
    for(id in 1:D){
      nd[id] <- length(d.list[[id]])
    }


    #separation des deux modalites
    xsep <- list()
    for(id in 1:D){
      xsep[[id]] <- x[,d.list[[id]]]
    }

    to.predict.sep <- list()
    for(id in 1:D){
      to.predict.sep[[id]] <- to.predict[,d.list[[id]]]
    }


    missing=FALSE
    if (sum(x==0)>0) {
      missing=TRUE
      # imputation aleatoire lors de l'init
      
      miss <- list()
      for(id in 1:D){
        miss[[id]] <- which(xsep[[id]]==0)
        xsep[[id]][miss[[id]]] <- sample(1:m[id], sum(xsep[[id]]==0), replace=TRUE)
        
      }
      
    }


    V=array(0,c(n,kr))       
    gamma=rep(0, kr)

    mus <- list()
    ps <- list() 
    res_mus <- list()
    res_ps <- list()
    for(id in 1:D){
      mus[[id]]=array(0,c(kr,nd[id],nbSEM+1))
      ps[[id]]=array(0,c(kr,nd[id],nbSEM+1))
      res_mus[[id]]=array(0,c(kr,nd[id]))
      res_ps[[id]]=array(0,c(kr,nd[id]))
    }
    

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

      # init des parametres a partir des partitions

      for(id in 1:D){
        for (h in 1:nd[id]){
          for (k in 1:kr){
            # init des parametres en fonction des partitions aleatoires
            # res=ordiem(as.vector(x[which(V[,k,1]==1),which(W[,l,1]==1)]),m,tabmu0=1:m,tabp0=seq(0,1,0.2),iter_max=iterordiEM)
            res <- ordiemCpp(m[id],tab_pejs[[id]],as.vector(xsep[[id]][which(V[,k]==1),h]),
                             tabmu0=1:m[id],tabp0=seq(0,1,0.2),
                             iter_max=iterordiEM)
            mus[[id]][k,h,1]=res[[1]]
            if(res[[2]]==1) {
              ps[[id]][k,h,1] = 0.999
            }
            else{
              ps[[id]][k,h,1]=res[[2]]
            }
          }
        }
      }
      


      # === init des valeurs manquantes ===
      if (missing){

        for(id in 1:D){
          xsep[[id]][miss[[id]]]=0 # on remet des 0 la ou il y avait des data manquantes
          for (h in 1:nd[id]){
            for (k in 1:kr){
              # recherche des cases manquantes
              tmp=which(xsep[[id]][which(V[,k]==1),h]==0)
              # simulation des data manquantes
              if (length(tmp)>0){
                probaBOS=rep(0,m[id])
                for (im in 1:m[id]) probaBOS[im]=(sum(tab_pejs[[id]][im,mus[[id]][k,h,1],]*(rep(ps[[id]][k,h,1],m[id])^(0:(m[id]-1)))))
                xsep[[id]][which(V[,k]==1),h][tmp]=sample(1:m[id],length(tmp),prob=probaBOS,replace=TRUE) 
              }
            }
          }
        }
        
      }
      # === debut du SEM ===
      for (iter in 1:nbSEM){
        if (disp) pb$tick()
        # ==== SE step ==== 
        # --- no need to compute proba for column and row partition because we know it
        # we just need this part for the missing values
        
        
        # --- imputation des donnees manquantes ----
        if (missing){
          for(id in 1:D){
            xsep[[id]][miss[[id]]]=0 # on remet des 0 la ou il y avait des data manquantes
            for (h in 1:nd[id]){
              for (k in 1:kr){
                # recherche des cases manquantes
                tmp=which(xsep[[id]][which(V[,k]==1),h]==0)
                # simulation des data manquantes
                if (length(tmp)>0){
                  probaBOS=rep(0,m[id])
                  for (im in 1:m[id]) probaBOS[im]=(sum(tab_pejs[[id]][im,mus[[id]][k,h,iter],]*(rep(ps[[id]][k,h,iter],m[id])^(0:(m[id]-1)))))
                  xsep[[id]][which(V[,k]==1),h][tmp]=sample(1:m[id],length(tmp),prob=probaBOS,replace=TRUE) 
                }
              }
            }
          }
        }
        
        # ==== M step ==== 
      
          for(id in 1:D){
            for (h in 1:nd[id]){
              for (k in 1:kr){
                xtmp=xsep[[id]]
                if (missing) xtmp[miss[[id]]]=0
                tmp=as.vector(xsep[[id]][which(V[,k]==1),h])
                datablock_kl=tmp[tmp>0]
                res <- ordiemCpp(m[id],tab_pejs[[id]],datablock_kl,
                          tabmu0=1:m[id],tabp0=ps[[id]][k,h,iter],
                          iter_max=iterordiEM)
                mus[[id]][k,h,iter+1]=res[[1]]
                if(res[[2]]==1){
                  ps[[id]][k,h,iter+1] = 0.999  
                }
                else{
                  ps[[id]][k,h,iter+1]=res[[2]]  
                }
              }
            }
          }
        
      }# for iter
      # ===== calcul des parametres (mode et median hors burn) =====
      for(id in 1:D){
        for (h in 1:nd[id]){
          for (k in 1:kr){
            res_mus[[id]][k,h]=mode(mus[[id]][k,h,nbSEMburn:(nbSEM+1)])
            res_ps[[id]][k,h]=median(ps[[id]][k,h,nbSEMburn:(nbSEM+1)])
          }
        }
      }
      
    
      # ===== estimation des partitions et des valeurs manquantes  =====
      if(disp) 
      {
        pb2 <- progress_bar$new(
          format="2/2 [:bar] :percent",
          total=nbSEM, clear=FALSE, width=60
        )
      }
      Q=nbSEM

      Xhat <- list()
      for(id in 1:D){
        Xhat[[id]]=array(0,c(n,nd[id],Q+1))
        Xhat[[id]][,,1:(Q+1)]=xsep[[id]]
      }
      
      for (iterQ in 1:Q){
        # --- no need for simulation about column  or row partitions
        if(disp) pb2$tick()
        # --- simulation des donnees manquantes ---
        if (missing){

          for(id in 1:D){
              tmpx=xsep[[id]]
              tmpx[miss[[id]]]=0 # on remet des 0 la ou il y avait des data manquantes
              for (h in 1:nd[id]){
                for (k in 1:kr){
                  # recherche des cases manquantes
                  tmp=which(tmpx[which(V[,k]==1),h]==0)
                  # simulation des data manquantes
                  if (length(tmp)>0){
                    probaBOS=rep(0,m[id])
                    for (im in 1:m[id]) probaBOS[im]=(sum(tab_pejs[[id]][im,res_mus[[id]][k,h],]*(rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1)))))
                    tmpx[which(V[,k]==1),h][tmp]=sample(1:m[id],length(tmp),prob=probaBOS,replace=TRUE) 
                  }
                }
              }
              Xhat[[id]][,,iterQ+1]=tmpx
            }
          }
          
      }#iterQ
      # --- estimation de la partition finale par mode marginal ---

      # --- estimation des valeurs manquantes par le mode marginal ---
      xhat <- list()
      for(id in 1:D){
        for (i in 1:n){
          for (h in 1:nd[id]){
            if (xsep[[id]][i,h]==0) xsep[[id]][i,h]=mode(Xhat[[id]][i,h,])
          }
        }

        # --- sauvegarde de la matrice completee ---
        xhat[[id]]=xsep[[id]]
      }
      

      # --- approximation de la loglik (non faite) ---
     
    
    # estimation des partitions
    zr=res_zr


    ####################### Don't think ICL is useful here #######################


    # Now that every thing is estimated, we can estimate the predictions for to.predict

    # First, impute missing values of validation dataset:
    missing.validation = FALSE
    if (sum(to.predict==0)>0) {
      missing.validation=TRUE
      miss.val <- list()
      for(id in 1:D){
        miss.val[[id]]=which(to.predict.sep[[id]]==0)
      }
      
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
        for(id in 1:D){
          to.predict.sep[[id]][miss.val[[id]]]=sample(1:m[id],length(miss.val[[id]]),replace=TRUE)
        }
        
    }

    #Then predict:
    n.to.predict= nrow(to.predict)
    logprobas <- matrix(0,nrow=n.to.predict,ncol=kr)
    probas <- matrix(0,nrow=n.to.predict,ncol=kr)

    logprobas=matrix(log(gamma),nrow=n.to.predict,ncol=kr,byrow = T)

    for(id in 1:D){
      for (k in 1:kr){
        for(h in 1:nd[id]){
          for (i in 1:n.to.predict){ 

            #logprobas[i,k] = logprobas[i,k] + 
            #        log(sum(tab_pejs[[id]][to.predict[i,h], res_mus[[id]][k,h],] *
            #                  (rep(res_ps[[id]][k,h],m)^(0:(m-1)))))
            
            sum.neg <- sum(tab_pejs[[id]][to.predict.sep[[id]][i,h], res_mus[[id]][k,h],] *
                              (rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1))))
            
            
            logprobas[i,k] = logprobas[i,k] + log(sum.neg)
                
      
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

    for(id in 1:D){
      xhat[[id]]=xhat[[id]]
      res_mus[[id]]=res_mus[[id]]
      res_ps[[id]]=res_ps[[id]]
      mus[[id]]=mus[[id]]
      ps[[id]]=ps[[id]]
    }

    # --- return the result ---
    res=list(probas.to.predict=probas,
             zr.to.predict=zr.to.predict,
             V.to.predict=V.to.predict,
             xhat=xhat,
             res_mus=res_mus,
             res_ps=res_ps,
             zr=zr,
             mus=mus,
             ps=ps,
             gamma=gamma,
             V=V,
             string="classifM.no.kc",
             m=m)
   
    return(res) 
  }
