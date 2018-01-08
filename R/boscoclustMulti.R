boscoclustMulti <-
function (x,kr,kc, m, d.list, nbSEM=50,nbSEMburn=20,nbindmini=4,init='kmeans',
  disp=TRUE,iterordiEM=10) {
    # setting progress bar
    if(disp) 
    {
      pb <- progress_bar$new(
        format="1/2 [:bar] :percent",
        total=nbSEM, clear=FALSE, width=60
      )
    }

    # to know how many super blocks we have:
    D <- length(m)
    # constant for polynomial probability (BOS) 
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


    # different levels separation
    xsep <- list()
    for(id in 1:D){
      xsep[[id]] <- x[,d.list[[id]]]
    }

    missing=FALSE
    if (sum(x==0)>0) {
      missing=TRUE
      
      
      miss <- list()
      for(id in 1:D){
        miss[[id]] <- which(xsep[[id]]==0)
        xsep[[id]][miss[[id]]] <- sample(1:m[id], sum(xsep[[id]]==0), replace=TRUE)
        
      }
      
    }
    loglik=rep(-Inf)
    # --- aleatory initialization for partitions ----
    # lines :
    V=array(0,c(n,kr,nbSEM+1))
    W <- list()
    for(id in 1:D){
      W[[id]] <- array(0, c(nd[id], kc[id],nbSEM+1))
    }
    # columns :  
    gamma <- array(0,c(kr,nbSEM+1))
    rho <- list()
    for(id in 1:D){
      rho[[id]] <- array(0, c(kc[id], nbSEM+1))
    }
  
    # params initialisation  
    mus <- list()
    ps <- list()
    res_mus <- list()
    res_ps <- list()
    for(id in 1:D){
      mus[[id]] <- array(0,c(kr, kc[id], nbSEM+1))
      ps[[id]] <- array(0,c(kr, kc[id], nbSEM+1))
      
      res_mus[[id]] <- array(0, c(kr, kc[id]))
      res_ps[[id]] <- array(0, c(kr, kc[id]))
    }


    res_gamma=matrix(0,kr)
    res_rho  <- list()
    for(id in 1:D){
      res_rho[[id]]=matrix(0,kc[id])
    }
   
    
    res_V=array(0,c(n,kr))
    res_W <- list()
    for(id in 1:D){
      res_W[[id]] <- array(0, c(nd[id], kc[id])) 
    }

    res_zr=matrix(0,n)
    res_zc <- list()
    for(id in 1:D){
      res_zc[[id]]<- matrix(0, nd[id])
    }


      # ==== init ==== 
      # --- aleatory initialization for partitions ----
      if (init=='random'){
        V[,,1]=t(rmultinom(n,1,rep(1/kr,kr)))
        
        for(id in 1:D){
          W[[id]][,,1]=t(rmultinom(nd[id],1,rep(1/kc[id],kc[id])))
        }

        
        verif <- TRUE
        verif.vec <- c()
        for(id in 1:D){
          verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],W[[id]][,,1],kc[id],nbindmini))
        }
        if(sum(verif.vec==0)>0){
          verif <- FALSE
        }
        
        while (! verif)
        {
          if(disp)print('reload random init')
          V[,,1]=t(rmultinom(n,1,rep(1/kr,kr)))
          for(id in 1:D){
            W[[id]][,,1]=t(rmultinom(nd[id],1,rep(1/kc[id],kc[id])))
          }
          
          
          verif.vec <- c()
          for(id in 1:D){
            verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],W[[id]][,,1],kc[id],nbindmini))
          }
          if(sum(verif.vec==0)>0){
            verif <- FALSE
          }
        }
      }

       # --- kmeans initialization for partitions ----
      if (init=='kmeans'){
        #x[miss]=NA
        #tmpV=kmeans(na.omit(x),kr,nstart=10)
        tmpV=kmeans(x,kr,nstart=10)
        for (i in 1:n) V[i,tmpV$cluster[i],1]=1
        #tmpW=kmeans(na.omit(t(x)),kc,nstart=10)
        
        tmpW <- list()
        for(id in 1:D){
          tmpW[[id]] <- kmeans(t(xsep[[id]]),kc[id],nstart=10)
          for (i in 1:nd[id]) W[[id]][i,tmpW[[id]]$cluster[i],1]=1
        }
          
        verif <- TRUE
        verif.vec <- c()
        for(id in 1:D){
          verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],W[[id]][,,1],kc[id],nbindmini))
        }
        if(sum(verif.vec==0)>0){
          verif <- FALSE
        }
        
        while (! verif)
        {
          if(disp) print('reload kmeans init')
          V[,,1]=0
          for(id in 1:D){
            W[[id]][,,1]=0
          }
          
          tmpV=kmeans(x,kr,nstart=10)
          for (i in 1:n) V[i,tmpV[i]$cluster,1]=1
          
          for(id in 1:D){
            tmpW[[id]]=kmeans(t(xsep[[id]]),kc[id],nstart=2)  
            for (i in 1:nd[id]) W[[id]][i,tmpW[[id]]$cluster[i],1]=1
          }
          
          verif <- TRUE
          verif.vec <- c()
          
          for(id in 1:D){
            verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],W[[id]][,,1],kc[id],nbindmini))
          }
          if(sum(verif.vec==0)>0){
            verif <- FALSE
          }
        }  
      }
      # ---- parameters initialization from partitions ----
      gamma[,1]=getMeans(V[,,1])
      for(id in 1:D){
        rho[[id]][,1]=getMeans(W[[id]][,,1])
      }
      

      for (k in 1:kr){

        for(id in 1:D){
          for(h in 1:kc[id]){
          
            res <- ordiemCpp(m[id],tab_pejs[[id]],as.vector(xsep[[id]][which(V[,k,1]==1),
                                               which(W[[id]][,h,1]==1)]),
                                     tabmu0=1:m[id],tabp0=seq(0,1,0.2),
                                    iter_max=iterordiEM)

            mus[[id]][k,h,1] <- res[[1]]

            if(res[[2]]==1) {
              ps[[id]][k,h,1] = 0.999
            }
            else{
              ps[[id]][k,h,1]=res[[2]]
            }
            
          }
        }
        
      }
      
     
      if (missing){
        
        
        for(id in 1:D){
          xsep[[id]][miss[[id]]] <- 0
        }
        

        for (k in 1:kr){
          
          for(id in 1:D){
            for(h in 1:kc[id]){
              
              tmp = which(xsep[[id]][which(V[,k,1]==1),which(W[[id]][,h,1]==1)]==0)
              
              if(length(tmp)>0){
                probaBOS = rep(0,m[id])
                for(im in 1:m[id]) {
                    probaBOS[im]=
                    (sum(tab_pejs[[id]][im, mus[[id]][k,h,1],]*
                           (rep(ps[[id]][k,h,1],m[id])^(0:(m[id]-1)))))                    
                }
                xsep[[id]][which(V[,k,1]==1), which(W[[id]][,h,1]==1)][tmp]= sample(1:m[id], length(tmp), prob=probaBOS, replace=T)
              }
            }
          }
        }
      }
      
      # ============  SEM ============
      for (iter in 1:nbSEM){
        if (disp) pb$tick()
        # ==== SE step ==== 
        
      logprobaV=matrix(log(gamma[,iter]),nrow=n,ncol=kr,byrow = T)
      logprobaW = list()
      for(id in 1:D){
        logprobaW[[id]]=matrix(log(rho[[id]][,iter]),nrow=nd[id],ncol=kc[id],byrow = T)
      }

      for(id in 1:D){
        for (k in 1:kr){
          for(h in 1:kc[id]){
            for (i in 1:n){ 
             for (tmp in 1:m[id]){
               logprobaV[i,k]=logprobaV[i,k] +  sum((xsep[[id]][i,]*W[[id]][,h,iter])==tmp)* log(sum(tab_pejs[[id]][tmp,mus[[id]][k,h,iter],]*
                                                                                                    (rep(ps[[id]][k,h,iter],m[id])^(0:(m[id]-1)))))     
             }    
            }
            for (d in 1:nd[id]){ 
              for (tmp in 1:m[id]){
                logprobaW[[id]][d,h]=logprobaW[[id]][d,h] + sum((xsep[[id]][,d]*V[,k,iter])==tmp) * log(sum(tab_pejs[[id]][tmp,mus[[id]][k,h,iter],]*(rep(ps[[id]][k,h,iter],m[id])^(0:(m[id]-1)))))     
              }
            }
          }
        }
      }
      
      
      if (iter>1) probaVold=probaV
      probaV=matrix(0,n,kr)

      
      for (i in 1:n){
        tmp=logsum (logprobaV[i,])
        probaV[i,]=exp ( logprobaV[i,] - tmp )
      }
      
      
      
      if(iter == 1){
        probaWold = list();
        probaW = list();
      }
      for(id in 1:D){
        if (iter>1) {probaWold[[id]]=probaW[[id]]}
        probaW[[id]]=matrix(0,nd[id],kc[id])
      }
      

      for(id in 1:D){
        for (d in 1:nd[id]){
          tmp=logsum (logprobaW[[id]][d,])
          probaW[[id]][d,]=exp ( logprobaW[[id]][d,] -  tmp)
        }
      }
      
      if (iter>1) {
        maxV=max(abs(probaVold-probaV))
        maxW = c();
        for(id in 1:D){
          maxW[id]=max(abs(probaWold[[id]]-probaW[[id]]))
        }
        
      }
    
        
        casevide=TRUE
        restart=0
        while (casevide && (restart<5000)){
          
          for(id in 1:D){
            for(d in 1:nd[id]){
                W[[id]][d,,iter+1] <- rmultinom(1,1,probaW[[id]][d,])

            }
          }
          
          for (i in 1:n){
            V[i,,iter+1]=rmultinom(1,1,probaV[i,])
          }
          casevide=FALSE
          verif.vec <- c()
          for(id in 1:D){
            verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,iter+1],W[[id]][,,iter+1],kc[id],nbindmini))
          }
          if(sum(verif.vec==0)>0) casevide <- TRUE

          if (casevide){
            restart=restart+1
            # cat('restart number ',restart,'\n')
          }
        }
        
        
        # --- missing vaues imputation ----
        if (missing){
          
          for(id in 1:D){
            xsep[[id]][miss[[id]]] <- 0
          }
          
          
          for (k in 1:kr){
            
            for(id in 1:D){
              for(h in 1:kc[id]){
                tmp <- which(xsep[[id]][which(V[,k,iter+1]==1), 
                                         which(W[[id]][,h,iter+1]==1)]==0)
                if(length(tmp)>0){
                  probaBOS <- rep(0,m[id])
                  for(im in 1:m[id]){
                    probaBOS[im] <- (sum(tab_pejs[[id]][im, mus[[id]][k,h,iter],]*
                                           (rep(ps[[id]][k,h,iter],m[id])^(0:(m[id]-1)))))
                    xsep[[id]][which(V[,k,iter+1]==1), 
                               which(W[[id]][,h,iter+1]==1)][tmp]=sample(1:m[id], length(tmp), prob=probaBOS, replace=T)
                  }
                }
              }
            }
          }
          
        }
        if ((restart==1000) && casevide){
          if(disp) print('The algorithm is stopped for degenerancy reason')
          return(NULL)
        }
        # ==== M step ==== 
      
          gamma[,iter+1]=getMeans(V[,,iter+1])
          
          for(id in 1:D){
            rho[[id]][,iter+1]=getMeans(W[[id]][,,iter+1])
          }

          for (k in 1:kr){
            
            for(id in 1:D){
              for(h in 1:kc[id]){
                res <- ordiemCpp(m[id], tab_pejs[[id]],as.vector(xsep[[id]][which(V[,k,iter+1]==1),
                                                   which(W[[id]][,h,iter+1]==1)]),
                              tabmu0 = 1:m[id], tabp0 = ps[[id]][k,h,iter],
                              iter_max = iterordiEM)
                mus[[id]][k,h,iter+1] <- res[[1]]
                if(res[[2]]==1) {
                  ps[[id]][k,h,iter+1] = 0.999
                }
                else{
                  ps[[id]][k,h,iter+1]=res[[2]]
                }
              }
            }
            
          }
        
      }# for iter
       #===== parameters computaton (mode and median after burn-in) =====
      for (k in 1:kr){
        
        for(id in 1:D){
          for(h in 1:kc[id]){
            if(k==1) res_rho[[id]][h] <- median(rho[[id]][h, nbSEMburn:(nbSEM+1)])
            if(h==1) res_gamma[k] <- median(gamma[k,nbSEMburn:(nbSEM+1)])
            res_mus[[id]][k,h] <- mode(mus[[id]][k,h,nbSEMburn:(nbSEM+1)])
            res_ps[[id]][k,h] <- median(ps[[id]][k,h,nbSEMburn:(nbSEM+1)])
          }
        }
      }
      
      for(id in 1:D){
        res_rho[[id]]=res_rho[[id]]/sum(res_rho[[id]])
      }
      
      res_gamma=res_gamma/sum(res_gamma)
      
      # --- partition simulation  ---

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
        Xhat[[id]] <- array(0, c(n,nd[id],Q+1))
        Wfinal[[id]] <- array(0,c(nd[id],kc[id],Q+1)) 
      }
      Vfinal <- array(0,c(n,kr,Q+1)) 
      
      
      
      Vfinal[,,1]=V[,,iter+1]
      for(id in 1:D){
        Wfinal[[id]][,,1]=W[[id]][,,iter+1]
      }
      
      for(id in 1:D){
        Xhat[[id]][,,1:(Q+1)] <- xsep[[id]]
      }

      for (iterQ in 1:Q){
        if(disp) pb2$tick()
        
        
        logprobaV=matrix(log(res_gamma),nrow=n,ncol=kr,byrow = T)
        logprobaW = list()
        for(id in 1:D){
          logprobaW[[id]]=matrix(log(res_rho[[id]]),nrow=nd[id],ncol=kc[id],byrow = T)
        }
        
        for (k in 1:kr){
          for(id in 1:D){
            for(h in 1:kc[id]){
              for (i in 1:n){ 
                for (tmp in 1:m[id]){
                  logprobaV[i,k]=logprobaV[i,k] +  sum((Xhat[[id]][i,,iterQ]* Wfinal[[id]][,h,iterQ])==tmp)  * 
                    log(sum(tab_pejs[[id]][tmp,res_mus[[id]][k,h],]*(rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1)))))     
                }    
              }
              for (d in 1:nd[id]){ 
                for (tmp in 1:m[id]){
                  logprobaW[[id]][d,h]=logprobaW[[id]][d,h] + sum((Xhat[[id]][,d,iterQ]* Vfinal[,k,iterQ])==tmp) * 
                    log(sum(tab_pejs[[id]][tmp,res_mus[[id]][k,h],]*(rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1)))))     
                }
              }
            }
          }
        }
        probaV=matrix(0,n,kr)
        for (i in 1:n){
          tmp=logsum (logprobaV[i,])
          probaV[i,]=exp ( logprobaV[i,] - tmp )
          Vfinal[i,,iterQ+1]=rmultinom(1,1,probaV[i,])
        }
        probaW <- list()
        for(id in 1:D){
          probaW[[id]]=matrix(0,nd[id],kc[id])
          for (d in 1:nd[id]){
            tmp=logsum (logprobaW[[id]][d,])
            probaW[[id]][d,]=exp ( logprobaW[[id]][d,] -  tmp)
            Wfinal[[id]][d,,iterQ+1]=rmultinom(1,1,probaW[[id]][d,])
          }
        }
        
      
       
        if (missing){
          tmpx <- list()
          
          for(id in 1:D){
            tmpx[[id]] <- xsep[[id]]
            tmpx[[id]][miss[[id]]] <- 0
          }
          
          for (k in 1:kr){
            
            for(id in 1:D){
              for(h in 1:kc[id]){
                
                tmp <- which(tmpx[[id]][which(Vfinal[,k,iterQ]==1), which(Wfinal[[id]][,h,iterQ]==1)]==0)
                # simulation for missing values
                if(length(tmp)>0){
                  probaBOS= rep(0,m[id])
                  for(im in 1:m[id]){
                    probaBOS[im] <- 
                      (sum(tab_pejs[[id]][im, res_mus[[id]][k,h],]*
                             (rep(res_ps[[id]][k,h],m[id])^(0:(m[id]-1)))))
                  }
                  tmpx[[id]][which(Vfinal[,k,iterQ]==1), which(Wfinal[[id]][,h,iterQ]==1)][tmp] <-
                    sample(1:m[id], length(tmp), prob=probaBOS, replace=TRUE)
                }
              }
              Xhat[[id]][,,iterQ+1]=tmpx[[id]]
            }

          }
          
        }# end if missing
      }#iterQ
      
      # --- final partition estimation  ---
      res_zr=apply(apply(Vfinal,c(1,2),sum),1,which.max)
      for(id in 1:D){
        res_zc[[id]] <- 
          apply(apply(Wfinal[[id]], c(1,2),sum), 1, which.max)
      }
      

      for (i in 1:n) res_V[i,res_zr[i]]=1
      
      for(id in 1:D){
        for(d in 1:nd[id]) res_W[[id]][d,res_zc[[id]][d]] <- 1
      }
      
      
     # --- missing values final estimation ---
      for (i in 1:n){
        
        for(id in 1:D){
          for(d in 1:nd[id]){
            if(xsep[[id]][i,d]==0) xsep[[id]][i,d] = mode(Xhat[[id]][i,d,])
          }
         
          Xhat[[id]] <- xsep[[id]]
        }

      }
      

    
    
    zr=res_zr
    zc <- list()
    for(id in 1:D){
      zc[[id]] <- res_zc[[id]]
    }

    #if(disp) print("computing ICL")
    # computing ICL
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
    # --- return the result ---
    
  
    
    
    for(id in 1:D){
      res_mus[[id]] <- res_mus[[id]]
      res_ps[[id]] <- res_ps[[id]]
      res_rho[[id]] <- res_rho[[id]]
      res_W[[id]] <- res_W[[id]]


      mus[[id]] <- mus[[id]]
      ps[[id]] <- ps[[id]]
      rho[[id]] <- rho[[id]]
      W[[id]] <- W[[id]]

    }
    
   # if(disp) print("return result")
    
    res=list(xhat=Xhat,
             mus=mus,
             ps=ps, 
             gamma=gamma,
             rho=rho,
             V=V,
             W=W,
             res_mus=res_mus,
             res_ps=res_ps, 
             res_gamma=res_gamma,
             res_rho=res_rho,
             res_V=res_V,
             res_W=res_W,
             icl=icl,
             zr=zr,
             zc=zc,
             probaV=probaV,
             probaW=probaW,
             string="coclustM",
             m=m)
   
    return(res)
  }
