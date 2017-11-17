bosclust <-
function (x,kr, m, nbSEM=50,nbSEMburn=20,nbindmini=4,init='kmeans',disp=TRUE,iterordiEM=10) {
    # ----------------------------------------------------------------------------
    # Estimation of the latent BOS clustering model via SEM algoritm
    # input
    #   x  : matrice n x p de donnees ordinales (individu en ligne, variable en colonne)
    #   kr : nb de classes en ligne  
    #   m  : vecteur de longueur D qui definit la modalite de chaque groupe
    #   d  : liste de longueur D qui defnit pour chaque modalite, quelles sont 
    #        les variables du dataset qui ont cette modalite
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

    # definir le nombre de modalite differentes :
    # on charge en memoire les exposants des probas BOS sous forme polynomiales
    tab_pej <- tabpej(m)
    

    
    # ---
    n=nrow(x)
    d = ncol(x)
    

   
    missing=FALSE
    if (sum(x==0)>0) {
      missing=TRUE
      # imputation aleatoire lors de l'init
      
        miss <- which(x==0)
        x[miss] <- sample(1:m, sum(x==0), replace=TRUE)
      
    }

    # initialisation des partitions 0 : 
    # en ligne :
    V=array(0,c(n,kr,nbSEM+1))  
    gamma <- array(0,c(kr,nbSEM+1))
    
  
    # initialisation des modes et precisions (mu et pi) 
    mu <- array(0,c(kr, d, nbSEM+1))
    p <- array(0,c(kr, d, nbSEM+1))
      
    res_mu <- array(0, c(kr, d))
    res_p <- array(0, c(kr, d))
    


    res_gamma=matrix(0,kr)
    res_V=array(0,c(n,kr))
    res_zr=matrix(0,n)
    

      # ==== init ==== 
      # --- init aleatoire des partitions ----
      if (init=='random'){
        V[,,1]=t(rmultinom(n,1,rep(1/kr,kr)))
        
        
        
        verif <- TRUE
        
        verif <- verif(x,V[,,1],rep(1,d),1,nbindmini)
        
        
        while (! verif)
        {
          print('reload random init')
          V[,,1]=t(rmultinom(n,1,rep(1/kr,kr)))
          
          verif <- verif(x,V[,,1],rep(1,d),1,nbindmini)

    
        }
      }

      # --- init kmeans des partitions ----
      if (init=='kmeans'){
        #x[miss]=NA
        #tmpV=kmeans(na.omit(x),kr,nstart=10)
        tmpV=kmeans(x,kr,nstart=10)
        for (i in 1:n) V[i,tmpV$cluster[i],1]=1
    
        verif <- TRUE
        verif <- verif(x,V[,,1],rep(1,d),1,nbindmini)
        
        while (! verif)
        {
          print('reload kmeans init')
          V[,,1]=0
          
          
          tmpV=kmeans(x,kr,nstart=10)
          for (i in 1:n) V[i,tmpV[i]$cluster,1]=1
          
          
          
          verif <- TRUE
          verif <- verif(x,V[,,1],rep(1,d),1,nbindmini)
          
        } 

      }
      # init des parametres a partir des partitions
      gamma[,1]=getMeans(V[,,1])
      
      

      for (k in 1:kr){

        
          for(h in 1:d){
            res <- ordiemCpp(m,tab_pej,as.vector(x[which(V[,k,1]==1),h]),
                                     tabmu0=1:m,tabp0=seq(0,1,0.2),
                                    iter_max=iterordiEM)
            mu[k,h,1] <- res[[1]]
            if(res[[2]]==1){
              p[k,h,1] <- 0.999
            }
            else{
              p[k,h,1] <- res[[2]]
            }
          }  
          
        
        
      }
      
      # === init des valeurs manquantes ===
      if (missing){
        
        # on remet des 0 la ou il y avait des data manquantes :
        x[miss] <- 0
        
        

        for (k in 1:kr){
          
          

            for(h in 1:d){
              # recherche des cases manquantes
              tmp = which(x[which(V[,k,1]==1),h]==0)
              # simulation des data manquantes
              if(length(tmp)>0){
                probaBOS = rep(0,m)
                for(im in 1:m) {
                    probaBOS[im]=
                    (sum(tab_pej[im, mu[k,h,1],]*
                           (rep(p[k,h,1],m)^(0:(m-1)))))                    
                }
                x[which(V[,k,1]==1),h][tmp]= sample(1:m, length(tmp), prob=probaBOS, replace=T)
              }
            }
              
            
          
        }
      }
      
      # === debut du SEM ===
      for (iter in 1:nbSEM){
        if (disp) pb$tick()
        # ==== SE step ==== 
        # --- calcul des probas pour la simulation de la partition en ligne
        logprobaV=matrix(0,n,kr)
        for (i in 1:n){
          for (k in 1:kr){
            logprobaV[i,k]=log(gamma[k,iter])
            
            
              for(h in 1:d){
                
                  logprobaV[i,k] = logprobaV[i,k] + 
                    log(sum(tab_pej[x[i,h], mu[k,h,iter],]*
                              (rep(p[k,h,iter],m)^(0:(m-1)))))
                
              }
            

          } 
        }
        if (iter>1) probaVold=probaV
        probaV=matrix(0,n,kr)
        for (i in 1:n){
          for (k in 1:kr){
            probaV[i,k]=exp ( logprobaV[i,k] - logsum (logprobaV[i,]) )
          }
        }
                
        
        
        
        casevide=TRUE
        restart=0
        while (casevide && (restart<5000)){
          
          
          for (i in 1:n){
            V[i,,iter+1]=rmultinom(1,1,probaV[i,])
          }
          casevide=FALSE
          verif <- TRUE
          verif <-  verif(x,V[,,iter+1],rep(1,d),1,nbindmini)
          
          if(!verif) casevide <- TRUE

          if (casevide){
            restart=restart+1
            # if(disp) cat('restart number ',restart,'\n')
          }
        }
        
        
        # --- imputation des donnees manquantes ----
        if (missing){
          # on remet des 0 la ou il y avait des data manquantes
          
           x[miss] <- 0
          
          
          
          for (k in 1:kr){
            
            
              for(h in 1:d){
                tmp <- which(x[which(V[,k,iter+1]==1),h]==0)
                if(length(tmp)>0){
                  probaBOS <- rep(0,m)
                  for(im in 1:m){
                    probaBOS[im] <- (sum(tab_pej[im, mu[k,h,iter],]*
                                           (rep(p[k,h,iter],m)^(0:(m-1)))))
                    x[which(V[,k,iter+1]==1),h][tmp]=sample(1:m, length(tmp), prob=probaBOS, replace=T)
                  }
                }
              }
                  
              
            
          }
          
        }
        if ((restart==1000) && casevide){
          print('The algorithm is stopped for degenerancy reason')
          return(NULL)
        }
        # ==== M step ==== 
        
          gamma[,iter+1]=getMeans(V[,,iter+1])
          
          for (k in 1:kr){
            
              
              for(h in 1:d){
                xtmp <- x
                if (missing) xtmp[miss] <- 0
                tmp <- as.vector(x[which(V[,k,iter+1]==1),h])
                datablock_kh <- tmp[tmp>0]
                res <- ordiemCpp(m, tab_pej, datablock_kh, tabmu0 = 1:m,
                              tabp0 = p[k,h,iter], iter_max = iterordiEM)
                
                mu[k,h,iter+1] <- res[[1]]
                
                if(res[[2]]==1){
                  p[k,h,iter+1] <- 0.999
                }
                else{
                  p[k,h,iter+1] <- res[[2]]
                }
                
              }
                
            }
          
        
        
      }# for iter
      # ===== calcul des parametres (mode et median hors burn) =====
      for (k in 1:kr){
        
           res_gamma[k] <- median(gamma[k,nbSEMburn:(nbSEM+1)])
          for(h in 1:d){
            res_mu[k,h] <- mode(mu[k,h,nbSEMburn:(nbSEM+1)])
            res_p[k,h] <- median(p[k,h,nbSEMburn:(nbSEM+1)])
          }
            
          
        
      }
      
      
      res_gamma=res_gamma/sum(res_gamma)
      
      # ===== estimation des partitions et des valeurs manquantes  =====

      if(disp) 
      {
        pb2 <- progress_bar$new(
          format="2/2 [:bar] :percent",
          total=nbSEM, clear=FALSE, width=60
        )
      }

      Q=nbSEM
      Xhat <- array(0, c(n,d,Q+1))
      
      Vfinal <- array(0,c(n,kr,Q+1)) 

      
      Vfinal[,,1]=V[,,iter+1]
      
      Xhat[,,1:(Q+1)] <- x

      for (iterQ in 1:Q){
        if(disp) pb2$tick()
        # --- simulation des partitions en ligne ---
        for (i in 1:n){
          for (k in 1:kr){
            logprobaV[i,k]=log(res_gamma[k])
            
              for(h in 1:d){
                
                  sum.neg <- sum(tab_pej[Xhat[i,h,iterQ], res_mu[k,h],]*
                              (rep(res_p[k,h], m)^(0:(m-1))))

                  if(sum.neg == 0){
                    sum.neg = 1e-300
                  }
                  logprobaV[i,k] <- logprobaV[i,k] + 
                    log(sum.neg)
                
              }
            
          } 
        }

        for (i in 1:n){
          for (k in 1:kr){
            probaV[i,k]=exp ( logprobaV[i,k] - logsum (logprobaV[i,]) )
          }

          Vfinal[i,,iterQ+1]=rmultinom(1,1,probaV[i,])

        }
        
      
        
        # --- simulation des donnees manquantes ---
        if (missing){
          tmpx <- list()
          
          
          tmpx <- x
          tmpx[miss] <- 0
          
          
          for (k in 1:kr){
            
              
              for(h in 1:d){
                  # recherche des cases manquantes
                  tmp <- which(tmpx[which(Vfinal[,k,iterQ]==1),h]==0)
                  # simulation des donnees manquantes 
                  if(length(tmp)>0){
                    probaBOS= rep(0,m)
                    for(im in 1:m){
                      probaBOS[im] <- 
                        (sum(tab_pej[im, res_mu[k,h],]*
                               (rep(res_p[k,h],m)^(0:(m-1)))))
                    }
                    tmpx[which(Vfinal[,k,iterQ]==1),][tmp] <-
                      sample(1:m, length(tmp), prob=probaBOS, replace=TRUE)
                  }
                
                Xhat[,,iterQ+1]=tmpx
              }
                
            

          }
          
        }# end if missing
      }#iterQ
      
      # --- estimation de la partition finale par mode marginal ---
      res_zr=apply(apply(Vfinal,c(1,2),sum),1,which.max)
      
      
      for (i in 1:n) res_V[i,res_zr[i]]=1
      
      
      
      # --- estimation des valeurs manquantes par le mode marginal ---
      for (i in 1:n){
        
          for(h in 1:d){
            if(x[i,h]==0) x[i,h] = mode(Xhat[i,h])
          }
          # --- sauvegarde de la matrice completee ---
          Xhat <- x
        
      }
      
    # estimation des partitions
    zr=res_zr
 

    
    #if(disp) print("computing ICL")
    # calcul ICL Brault
    # formule adaptee a notre cas ( 1 au lieu de m-1 param par case), mais a verifier theoriquement
    if (!missing){

      
      icl=- (kr-1)/2 *log(n)

      


      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log (res_gamma[k])
        }
      }
      
      
        for(h in 1:d){

            for(i in 1:n){
              for(k in 1:kr){
                icl <- icl +  res_V[i,k] *
                  log ( (sum(tab_pej[x[i,h],res_mu[k,h],]*
                               (rep(res_p[k,h],m)^(0:(m-1))))))
              }
            }
          
        }
      
    }
    
    else{

      icl=- (kr-1)/2 *log(n)

      


      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log (res_gamma[k])
        }
      }
      


        for(h in 1:d){

            for(i in 1:n){
              for(k in 1:kr){
                icl <- icl + res_V[i,k] *
                  log ( (sum(tab_pej[x[i,h],res_mu[k,h],]*
                               (rep(res_p[k,h],m)^(0:(m-1))))))
              }
            }
          
        }
      


      
    }#if(!missing)
    # --- return the result ---


      res_mu <- res_mu
      res_p <- res_p

      mu <- mu
      p <- p
    
    res=list(xhat=Xhat, 
             mu=mu, 
             p=p, 
             gamma=gamma,
             V=V,
             res_mu=res_mu,
             res_p=res_p, 
             res_gamma=res_gamma,
             res_V=res_V,
             icl=icl,
             zr=zr,
             probaV=probaV,
             string="clust",
             m=m)

    return(res)
  }
