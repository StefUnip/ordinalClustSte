bosclassifMulti.kc <-
  function (x,y,to.predict,d.list,kr,kc=1,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
    init='kmeans',disp=TRUE,iterordiEM=10) {
    # ----------------------------------------------------------------------------
    # Estimation of the latent BOS coclustering model via SEM algoritm
    # input
    #   x  : matrice n x p de donnees ordinales (individu en ligne, variable en colonne)
    #   kc : nb de classes en colonne
    #   kr : nb de classes en ligne  
    #   m  : nombre de modalites (vector) 
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

    # to know hoy many super blocks we have:
    D <- length(m)

    # on charge en memoire les exposants des probas BOS sous forme polynomiales
    tab_pejs <-list()
    for(id in 1:D){
      tab_pejs[[id]]=tabpej(m[id])
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
    res_gamma=rep(0, kr)
    res_V=array(0,c(n,kr)) 
    res_zr=y


    W <- list()
    rho <- list()
    mus <- list()
    ps <- list()
    res_mus <- list()
    res_ps <- list()
    res_rho <- list()
    res_W <-list()
    res_zc <- list()
    for(id in 1:D){
      W[[id]]=array(0,c(nd[id],kc[id],nbSEM+1))  
      rho[[id]]=array(0,c(kc[id],nbSEM+1))
      mus[[id]]=array(0,c(kr,kc[id],nbSEM+1))
      ps[[id]]=array(0,c(kr,kc[id],nbSEM+1))
      res_mus[[id]]=array(0,c(kr,kc[id]))
      res_ps[[id]]=array(0,c(kr,kc[id]))
      res_rho[[id]]=matrix(0,kc[id])
      res_W[[id]]=array(0,c(nd[id],kc[id])) 
      res_zc[[id]]=matrix(0,nd[id])
    }
          

    

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
        for(id in 1:D){
          W[[id]][,,1]=t(rmultinom(nd[id],1,rep(1/kc[id],kc[id])))
          while (! verif(xsep[[id]],V,W[[id]][,,1],kc[id],nbindmini))
          {
            print('reload random init')
            W[[id]][,,1]=t(rmultinom(nd[id],1,rep(1/kc[id],kc[id])))
          }
        }
        
      }
      # --- init kmeans des partitions ----
      if (init=='kmeans'){
        for(id in 1:D){
          tmpW=kmeans(t(xsep[[id]]),kc[id],nstart=10)  
          for (i in 1:nd[id]) W[[id]][i,tmpW$cluster[i],1]=1
          while (! verif(xsep[[id]],V,W[[id]][,,1],kc[id],nbindmini))
          {
            print('reload kmeans init')
            W[[id]][,,1]=0
            tmpW=kmeans(t(xsep[[id]]),kc[id],nstart=2)  
            for (i in 1:nd[id]) W[[id]][i,tmpW[i]$cluster,1]=1
          }
        }
          
      }
      # init des parametres a partir des partitions
      for(id in 1:D){
        rho[[id]][,1]=getMeans(W[[id]][,,1])
        for (l in 1:kc[id]){
          for (k in 1:kr){
            # init des parametres en fonction des partitions aleatoires
            # res=ordiem(as.vector(x[which(V[,k,1]==1),which(W[,l,1]==1)]),m,tabmu0=1:m,tabp0=seq(0,1,0.2),iter_max=iterordiEM)
            res <- ordiemCpp(m[id],tab_pejs[[id]],as.vector(xsep[[id]][which(V[,k]==1),
                                                   which(W[[id]][,l,1]==1)]),
                             tabmu0=1:m[id],tabp0=seq(0,1,0.2),
                             iter_max=iterordiEM)
            mus[[id]][k,l,1]=res[[1]]
            if(res[[2]]==1) {
              ps[[id]][k,l,1] = 0.999
            }
            else{
              ps[[id]][k,l,1]=res[[2]]
            }
          }
        }
      }
      
      # === init des valeurs manquantes ===
      if (missing){
        for(id in 1:D){
          xsep[[id]][miss[[id]]]=0 # on remet des 0 la ou il y avait des data manquantes
          for (l in 1:kc[id]){
            for (k in 1:kr){
              # recherche des cases manquantes
              tmp=which(xsep[[id]][which(V[,k]==1),which(W[[id]][,l,1]==1)]==0)
              # simulation des data manquantes
              if (length(tmp)>0){
                probaBOS=rep(0,m[id])
                for (im in 1:m[id]) probaBOS[im]=(sum(tab_pejs[[id]][im,mus[[id]][k,l,1],]*(rep(ps[[id]][k,l,1],m[id])^(0:(m[id]-1)))))
                xsep[[id]][which(V[,k]==1),which(W[[id]][,l,1]==1)][tmp]=sample(1:m[id],length(tmp),prob=probaBOS,replace=TRUE) 
              }
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

        logprobaW <- list()
        for(id in 1:D){
         logprobaW[[id]]=matrix(log(rho[[id]][,iter]),nrow=nd[id],ncol=kc[id],byrow = T)
         for (k in 1:kr){
           for(l in 1:kc[id]){
             for (h in 1:nd[id]){ 
               for (tmp in 1:m[id]){
                 logprobaW[[id]][h,l] =
                 logprobaW[[id]][h,l] + sum((xsep[[id]][,h]*V[,k])==tmp) * 
                 log(sum(tab_pejs[[id]][tmp,mus[[id]][k,l,iter],]*(rep(ps[[id]][k,l,iter],m[id])^(0:(m[id]-1)))))     
               }
             }
           }
         }
        }

        
        
        
        if (iter>1) probaWold=probaW
        probaW <- list()
        for(id in 1:D){
          probaW[[id]]=matrix(0,nd[id],kc[id])
          for (h in 1:nd[id]){
            tmp=logsum (logprobaW[[id]][h,])
            probaW[[id]][h,]=exp ( logprobaW[[id]][h,] -  tmp)
          }
        }
        
        
        
        if (iter>1) {
            maxW = c();
            for(id in 1:D){
              maxW[id]=max(abs(probaWold[[id]]-probaW[[id]]))
            } 
        }


        casevide=TRUE
        restart=0
        while (casevide && (restart<1000)){
          for(id in 1:D){
            for (h in 1:nd[id]){
              W[[id]][h,,iter+1]=rmultinom(1,1,probaW[[id]][h,])
            }
            
            casevide = (! verif(xsep[[id]],V,W[[id]][,,iter+1],kc[id],nbindmini))
            if (casevide){
              restart=restart+1
              cat('restart number ',restart,'\n')
            }
          }
          
        }


        # --- imputation des donnees manquantes ----
        if (missing){
          for(id in 1:D){
            xsep[[id]][miss[[id]]]=0 # on remet des 0 la ou il y avait des data manquantes
            for (l in 1:kc[id]){
              for (k in 1:kr){
                # recherche des cases manquantes
                tmp=which(xsep[[id]][which(V[,k]==1),which(W[[id]][,l,iter+1]==1)]==0)
                # simulation des data manquantes
                if (length(tmp)>0){
                  probaBOS=rep(0,m[id])
                  for (im in 1:m[id]) probaBOS[im]=(sum(tab_pejs[[id]][im,mus[[id]][k,l,iter],]*(rep(ps[[id]][k,l,iter],m[id])^(0:(m[id]-1)))))
                  xsep[[id]][which(V[,k]==1),which(W[[id]][,l,iter+1]==1)][tmp]=sample(1:m[id],length(tmp),prob=probaBOS,replace=TRUE) 
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
        
          for(id in 1:D){
            rho[[id]][,iter+1]=getMeans(W[[id]][,,iter+1])
            for (l in 1:kc[id]){
              for (k in 1:kr){
                xtmp=xsep[[id]]
                if (missing) xtmp[miss[[id]]]=0
                tmp=as.vector(xsep[[id]][which(V[,k]==1),which(W[[id]][,l,iter+1]==1)])
                datablock_kl=tmp[tmp>0]
                # res=ordiem(datablock_kl,m,tabmu0=1:m,tabp0=p[k,l,iter],iter_max=iterordiEM)
                res <- ordiemCpp(m[id],tab_pejs[[id]],datablock_kl,
                          tabmu0=1:m[id],tabp0=ps[[id]][k,l,iter],
                          iter_max=iterordiEM)
                mus[[id]][k,l,iter+1]=res[[1]]
                if(res[[2]]==1) {
                  ps[[id]][k,l,iter+1] = 0.999
                }
                else{
                  ps[[id]][k,l,iter+1]=res[[2]]
                }
              }
            }
          }
          
        
      }# for iter
      # ===== calcul des parametres (mode et median hors burn) =====
      for(id in 1:D){
        for (l in 1:kc[id]){
          res_rho[[id]][l]=median(rho[[id]][l,nbSEMburn:(nbSEM+1)])
          for (k in 1:kr){
            res_mus[[id]][k,l]=mode(mus[[id]][k,l,nbSEMburn:(nbSEM+1)])
            res_ps[[id]][k,l]=median(ps[[id]][k,l,nbSEMburn:(nbSEM+1)])
          }
        }
        res_rho[[id]]=res_rho[[id]]/sum(res_rho[[id]])
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
      Wfinal <- list()
      for(id in 1:D){
        Xhat[[id]]=array(0,c(n,nd[id],Q+1))
        Wfinal[[id]]=array(0,c(nd[id],kc[id],Q+1)) 
        Wfinal[[id]][,,1]=W[[id]][,,iter+1]
        Xhat[[id]][,,1:(Q+1)]=xsep[[id]]
      }
      
      for (iterQ in 1:Q){
        if(disp) pb2$tick()
        # --- simulation des partitions en ligne et en colonne ---
        for(id in 1:D){
          logprobaW[[id]]=matrix(log(res_rho[[id]]),nrow=nd[id],ncol=kc[id],byrow = T)
          for (k in 1:kr){
            for(l in 1:kc[id]){
              for (h in 1:nd[id]){ 
                for (tmp in 1:m[id]){
                  logprobaW[[id]][h,l] =
                  logprobaW[[id]][h,l] + sum((Xhat[[id]][,h,iterQ]* V[,k])==tmp) * 
                  log(sum(tab_pejs[[id]][tmp,res_mus[[id]][k,l],]*(rep(res_ps[[id]][k,l],m[id])^(0:(m[id]-1)))))     
                }
              }
            }
          }
        }

        

        for(id in 1:D){

          for (h in 1:nd[id]){
            tmp=logsum (logprobaW[[id]][h,])
            probaW[[id]][h,]=exp ( logprobaW[[id]][h,] -  tmp)
            Wfinal[[id]][h,,iterQ+1]=rmultinom(1,1,probaW[[id]][h,])
          }
        }
        
        # --- simulation des donnees manquantes ---
        if (missing){
          for(id in 1:D){
            tmpx=xsep[[id]]
            tmpx[miss[[id]]]=0 # on remet des 0 la ou il y avait des data manquantes
            for (l in 1:kc[id]){
              for (k in 1:kr){
                # recherche des cases manquantes
                tmp=which(tmpx[which(V[,k]==1),which(Wfinal[[id]][,l,iterQ]==1)]==0)
                # simulation des data manquantes
                if (length(tmp)>0){
                  probaBOS=rep(0,m[id])
                  for (im in 1:m[id]) probaBOS[im]=(sum(tab_pejs[[id]][im,res_mus[[id]][k,l],]*(rep(res_ps[[id]][k,l],m[id])^(0:(m[id]-1)))))
                  tmpx[which(V[,k]==1),which(Wfinal[[id]][,l,iterQ]==1)][tmp]=sample(1:m[id],length(tmp),prob=probaBOS,replace=TRUE) 
                }
              }
            }

            Xhat[[id]][,,iterQ+1]=tmpx
          }
          
        }
      }#iterQ
      # --- estimation de la partition finale par mode marginal ---


      for(id in 1:D){
        res_zc[[id]]=apply(apply(Wfinal[[id]],c(1,2),sum),1,which.max)
        for (h in 1:nd[id]) res_W[[id]][h,res_zc[[id]][h]]=1
      }
      

      
      # --- estimation des valeurs manquantes par le mode marginal ---
      for(id in 1:D){
        for (i in 1:n){
          for (h in 1:nd[id]){
            if (xsep[[id]][i,h]==0) xsep[[id]][i,h]=mode(Xhat[[id]][i,h,])
          }
        }
      }
      

      # --- sauvegarde de la matrice completee ---
      xhat <- list()
      for(id in 1:D){
        xhat[[id]]=xsep[[id]]
      }
      
      # --- approximation de la loglik (non faite) ---
      
    # estimation des partitions
    zr=res_zr
    zc <- list()
    for(id in 1:D){
      zc[[id]]=res_zc[[id]]
    }
    
    #if(disp) print("computing ICL")
    # calcul ICL Brault
    # formule adaptee a notre cas ( 1 au lieu de m-1 param par case), mais a verifier theoriquement
    if (!missing){

      
      icl=- (kr-1)/2 *log(n)

      for(id in 1:D){
        d = nd[id]
        # print(icl)
        icl <- icl - (kc[id]-1)/2 *log(d) - kc[id]*kr/2 *log(n*d)
      }


      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log (res_gamma[k])
        }
      }
      
      
      for(id in 1:D){
        for (d in 1:nd[id]){
          for (h in 1:kc[id]){
            icl = icl + res_W[[id]][d,h] * log (res_rho[[id]][h])
          }
        }
      }

      
      for(id in 1:D){
        for(d in 1:nd[id]){
          for(h in 1:kc[id]){
            for(i in 1:n){
              for(k in 1:kr){
                icl <- icl + res_W[[id]][d,h] * res_V[i,k] *
                  log ( (sum(tab_pejs[[id]][xsep[[id]][i,d],res_mus[[id]][k,h],]*
                               (rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1))))))
              }
            }
          }
        }
      } 
    }
    
    else{

      icl=- (kr-1)/2 *log(n)


      for(id in 1:D){
        d = nd[id]
        # print(icl)
        icl <- icl - (kc[id]-1)/2 *log(d) - kc[id]*kr/2 *log(n*d)
      }
      


      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log (res_gamma[k])
        }
      }
      

      for(id in 1:D){
        for (d in 1:nd[id]){
          for (h in 1:kc[id]){
            icl = icl + res_W[[id]][d,h] * log (res_rho[[id]][h])
          }
        }
      }
      
      
      for(id in 1:D){
        # xsep[[id]][miss[[id]]] <- 0
        for(d in 1:nd[id]){
          for(h in 1:kc[id]){
            for(i in 1:n){
              for(k in 1:kr){
                icl <- icl + res_W[[id]][d,h] * res_V[i,k] *
                  log ( (sum(tab_pejs[[id]][xsep[[id]][i,d],res_mus[[id]][k,h],]*
                               (rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1))))))
              }
            }
          }
        }
      }

    }#if(!missing)

    # Now that every thing is estimated, we estimate the predictions for to.predict

    n.to.predict= nrow(to.predict)
    logprobas <- matrix(0,nrow=n.to.predict,ncol=kr)
    probas <- matrix(0,nrow=n.to.predict,ncol=kr)

    logprobas=matrix(log(gamma),nrow=n.to.predict,ncol=kr,byrow = T)
    for(id in 1:D){
      for (k in 1:kr){
        for(l in 1:kc[id]){
          for (i in 1:n.to.predict){ 
           for (tmp in 1:m[id]){
             logprobas[i,k]=logprobas[i,k] +  
               sum((to.predict.sep[[id]][i,]*res_W[[id]][,l])==tmp)* 
               log(sum(tab_pejs[[id]][tmp,res_mus[[id]][k,l],]*
               (rep(res_ps[[id]][k,l],m[id])^(0:(m[id]-1)))))     
           }    
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
      res_rho[[id]]=res_rho[[id]]
      res_W[[id]]=res_W[[id]]
      mus[[id]]=mus[[id]]
      ps[[id]]=ps[[id]]
      rho[[id]]=rho[[id]]
      W[[id]]=W[[id]]
    }



    # --- return the result ---
    res=list(probas.to.predict=probas,
             zr.to.predict=zr.to.predict,
             V.to.predict=V.to.predict,
             xhat=xhat,
             res_mus=res_mus,
             res_ps=res_ps,
             res_rho=res_rho,
             res_W=res_W,
             icl=icl,
             zr=zr,
             zc=zc,
             mus=mus,
             ps=ps,
             gamma=gamma,
             rho=rho,
             V=V,
             W=W,
             probaW=probaW,
             string="classifM.kc",
             m=m)

    return(res) 
  }
