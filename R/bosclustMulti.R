bosclustMulti <-
function (x,kr, m, d.list, nbSEM=50,nbSEMburn=20,nbindmini=4,
        init='kmeans',disp=TRUE,iterordiEM=10) {
    
    
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

  
    V=array(0,c(n,kr,nbSEM+1))  
    gamma <- array(0,c(kr,nbSEM+1))
    
  
   
    mus <- list()
    ps <- list()
    res_mus <- list()
    res_ps <- list()
    for(id in 1:D){
      mus[[id]] <- array(0,c(kr, nd[id], nbSEM+1))
      ps[[id]] <- array(0,c(kr, nd[id], nbSEM+1))
      
      res_mus[[id]] <- array(0, c(kr, nd[id]))
      res_ps[[id]] <- array(0, c(kr, nd[id]))
    }


    res_gamma=matrix(0,kr)
    res_V=array(0,c(n,kr))
    res_zr=matrix(0,n)
    

      # ==== init ==== 
       # --- aleatory initialization for partitions ----
      if (init=='random'){
        V[,,1]=t(rmultinom(n,1,rep(1/kr,kr)))
        
        
        
        verif <- TRUE
        verif.vec <- c()
        for(id in 1:D){
          verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],rep(1,length(d.list[[id]])),1,nbindmini))
        }
        if(sum(verif.vec==0)>0){
          verif <- FALSE
        }
        
        while (! verif)
        {
          #print('reload random init')
          V[,,1]=t(rmultinom(n,1,rep(1/kr,kr)))
          
          
          verif.vec <- c()
          for(id in 1:D){
            verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],rep(1,length(d.list[[id]])),1,nbindmini))
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
    
        verif <- TRUE
        verif.vec <- c()
        for(id in 1:D){
          verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],rep(1,length(d.list[[id]])),1,nbindmini))
        }
        if(sum(verif.vec==0)>0){
          verif <- FALSE
        }
        
        while (! verif)
        {
          #print('reload kmeans init')
          V[,,1]=0
          
          
          tmpV=kmeans(x,kr,nstart=10)
          for (i in 1:n) V[i,tmpV[i]$cluster,1]=1
          
          
          
          verif <- TRUE
          verif.vec <- c()
          
          for(id in 1:D){
            verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,1],rep(1,length(d.list[[id]])),1,nbindmini))
          }
          if(sum(verif.vec==0)>0){
            verif <- FALSE
          }
        } 

      }
       # ---- parameters initialization from partitions ----
      gamma[,1]=getMeans(V[,,1])
      
      

      for (k in 1:kr){

        for(id in 1:D){
          for(d in 1:nd[id]){
            res <- ordiemCpp(m[id],tab_pejs[[id]],as.vector(xsep[[id]][which(V[,k,1]==1),d]),
                                     tabmu0=1:m[id],tabp0=seq(0,1,0.2),
                                    iter_max=iterordiEM)
            mus[[id]][k,d,1] <- res[[1]]
            if(res[[2]]==1){
              ps[[id]][k,d,1] <- 0.999
            }
            else{
              ps[[id]][k,d,1] <- res[[2]]
            }
          }  
          
        }
        
      }
      
     # ---- missing values initialization ----
      if (missing){
        
        
        for(id in 1:D){
          xsep[[id]][miss[[id]]] <- 0
        }
        

        for (k in 1:kr){
          
          for(id in 1:D){

            for(d in 1:nd[id]){
              
              tmp = which(xsep[[id]][which(V[,k,1]==1),d]==0)
              # simulation
              if(length(tmp)>0){
                probaBOS = rep(0,m[id])
                for(im in 1:m[id]) {
                    probaBOS[im]=
                    (sum(tab_pejs[[id]][im, mus[[id]][k,d,1],]*
                           (rep(ps[[id]][k,d,1],m[id])^(0:(m[id]-1)))))                    
                }
                xsep[[id]][which(V[,k,1]==1),d][tmp]= sample(1:m[id], length(tmp), prob=probaBOS, replace=T)
              }
            }
              
            
          }
        }
      }
      
      # ============  SEM ============
      for (iter in 1:nbSEM){
        if (disp) pb$tick()
        # ==== SE step ==== 
       
        logprobaV=matrix(0,n,kr)
        for (i in 1:n){
          for (k in 1:kr){
            logprobaV[i,k]=log(gamma[k,iter])
            
            for(id in 1:D){
              for(d in 1:nd[id]){
                
                  logprobaV[i,k] = logprobaV[i,k] + 
                    log(sum(tab_pejs[[id]][xsep[[id]][i,d], mus[[id]][k,d,iter],]*
                              (rep(ps[[id]][k,d,iter],m[id])^(0:(m[id]-1)))))
                
              }
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
          verif.vec <- c()
          for(id in 1:D){
            verif.vec <- c(verif.vec, verif(xsep[[id]],V[,,iter+1],rep(1,length(d.list[[id]])),1,nbindmini))
          }
          if(sum(verif.vec==0)>0) casevide <- TRUE

          if (casevide){
            restart=restart+1
            # if(disp) cat('restart number ',restart,'\n')
          }
        }
        
        
       # --- missing values imputation ----
        if (missing){
          
          for(id in 1:D){
            xsep[[id]][miss[[id]]] <- 0
          }
          
          
          for (k in 1:kr){
            
            for(id in 1:D){
              for(d in 1:nd[id]){
                tmp <- which(xsep[[id]][which(V[,k,iter+1]==1),d]==0)
                if(length(tmp)>0){
                  probaBOS <- rep(0,m[id])
                  for(im in 1:m[id]){
                    probaBOS[im] <- (sum(tab_pejs[[id]][im, mus[[id]][k,d,iter],]*
                                           (rep(ps[[id]][k,d,iter],m[id])^(0:(m[id]-1)))))
                    xsep[[id]][which(V[,k,iter+1]==1),d][tmp]=sample(1:m[id], length(tmp), prob=probaBOS, replace=T)
                  }
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
            
            for(id in 1:D){
              
              for(d in 1:nd[id]){
                xtmp <- xsep[[id]]
                if (missing) xtmp[miss[[id]]] <- 0
                tmp <- as.vector(xsep[[id]][which(V[,k,iter+1]==1),d])
                datablock_kh <- tmp[tmp>0]
                res <- ordiemCpp(m[id], tab_pejs[[id]], datablock_kh, tabmu0 = 1:m[id],
                              tabp0 = ps[[id]][k,d,iter], iter_max = iterordiEM)
                
                mus[[id]][k,d,iter+1] <- res[[1]]
                
                if(res[[2]]==1){
                  ps[[id]][k,d,iter+1] <- 0.999
                }
                else{
                  ps[[id]][k,d,iter+1] <- res[[2]]
                }
                
              }
                
            }
          }
        
        
      }# for iter
       #===== parameters computaton (mode and median after burn-in) =====
      for (k in 1:kr){
        
        for(id in 1:D){
          if(id==1) res_gamma[k] <- median(gamma[k,nbSEMburn:(nbSEM+1)])
          for(d in 1:nd[id]){
            res_mus[[id]][k,d] <- mode(mus[[id]][k,d,nbSEMburn:(nbSEM+1)])
            res_ps[[id]][k,d] <- median(ps[[id]][k,d,nbSEMburn:(nbSEM+1)])
          }
            
          
        }
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
      for(id in 1:D){
        Xhat[[id]] <- array(0, c(n,nd[id],Q+1))
      }
      Vfinal <- array(0,c(n,kr,Q+1)) 
      
      
      
      Vfinal[,,1]=V[,,iter+1]
      
      
      for(id in 1:D){
        Xhat[[id]][,,1:(Q+1)] <- xsep[[id]]
      }

      for (iterQ in 1:Q){
        if(disp) pb2$tick()
        for (i in 1:n){
          for (k in 1:kr){
            logprobaV[i,k]=log(res_gamma[k])
            
            for(id in 1:D){
              for(d in 1:nd[id]){
                
                  sum.neg <- sum(tab_pejs[[id]][Xhat[[id]][i,d,iterQ], res_mus[[id]][k,d],]*
                              (rep(res_ps[[id]][k,d], m[id])^(0:(m[id]-1))))

                  if(sum.neg == 0){
                    sum.neg = 1e-300
                  }
                  logprobaV[i,k] <- logprobaV[i,k] + 
                    log(sum.neg)
                
              }
            }
          } 
        }

        for (i in 1:n){
          for (k in 1:kr){
            probaV[i,k]=exp ( logprobaV[i,k] - logsum (logprobaV[i,]) )
          }

          Vfinal[i,,iterQ+1]=rmultinom(1,1,probaV[i,])

        }
        
      
        
         # --- missing values imputation ----
        if (missing){
          tmpx <- list()
          
          for(id in 1:D){
            tmpx[[id]] <- xsep[[id]]
            tmpx[[id]][miss[[id]]] <- 0
          }
          
          for (k in 1:kr){
            
            for(id in 1:D){
              
              for(d in 1:nd[id]){
                  
                  tmp <- which(tmpx[[id]][which(Vfinal[,k,iterQ]==1),d]==0)
                  
                  if(length(tmp)>0){
                    probaBOS= rep(0,m[id])
                    for(im in 1:m[id]){
                      probaBOS[im] <- 
                        (sum(tab_pejs[[id]][im, res_mus[[id]][k,d],]*
                               (rep(res_ps[[id]][k,d],m[id])^(0:(m[id]-1)))))
                    }
                    tmpx[[id]][which(Vfinal[,k,iterQ]==1),][tmp] <-
                      sample(1:m[id], length(tmp), prob=probaBOS, replace=TRUE)
                  }
                
                Xhat[[id]][,,iterQ+1]=tmpx[[id]]
              }
                
            }

          }
          
        }# end if missing
      }#iterQ
      
       # --- final partition estimation  ---

      res_zr=apply(apply(Vfinal,c(1,2),sum),1,which.max)
      
      
      for (i in 1:n) res_V[i,res_zr[i]]=1
      
      
      
      # --- missing values final estimation ---
      for (i in 1:n){
        
        for(id in 1:D){
          for(d in 1:nd[id]){
            if(xsep[[id]][i,d]==0) xsep[[id]][i,d] = mode(Xhat[[id]][i,d])
          }
          # --- sauvegarde de la matrice completee ---
          Xhat[[id]] <- xsep[[id]]
        }
      }
      

    # estimation des partitions
    zr=res_zr
 

    
    #if(disp) print("computing ICL")
    # computing ICL
    if (!missing){

      
      icl=- (kr-1)/2 *log(n)

      #icl=- (kr-1)/2 *log(n) - (d-1)/2 *log(d)- d*kr/2 *log(n*d)

      for(id in 1:D){
        icl = icl - (nd[id]-1)/2 *log(nd[id])- nd[id]*kr/2 *log(n*nd[id])
      }


      for(i in 1:n){
        for(k in 1:kr){
          icl = icl + res_V[i,k] * log(res_gamma[k])
        }
      }

      for(id in 1:D){
        icl = icl +  nd[id]*log(1/nd[id])  
      }
      
      

      
      for(id in 1:D){
        for(d in 1:nd[id]){

            for(i in 1:n){
              for(k in 1:kr){
                icl <- icl +  res_V[i,k] *
                  log ( (sum(tab_pejs[[id]][xsep[[id]][i,d],res_mus[[id]][k,d],]*
                               (rep(res_ps[[id]][k,d],m[id])^(0:(m[id]-1))))))
              }
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
      


      
      
      for(id in 1:D){
        # xsep[[id]][miss[[id]]] <- 0
        for(d in 1:nd[id]){

            for(i in 1:n){
              for(k in 1:kr){
                icl <- icl + res_V[i,k] *
                  log ( (sum(tab_pejs[[id]][xsep[[id]][i,d],res_mus[[id]][k,d],]*
                               (rep(res_ps[[id]][k,d],m[id])^(0:(m[id]-1))))))
              }
            }
          
        }
      }


      
    }#if(!missing)
    # --- return the result ---


    for(id in 1:D){
      res_mus[[id]] <- res_mus[[id]]
      res_ps[[id]] <- res_ps[[id]]

      mus[[id]] <- mus[[id]]
      ps[[id]] <- ps[[id]]
    }
    res=list(xhat=Xhat, 
             mus=mus, 
             ps=ps, 
             gamma=gamma,
             V=V,
             res_mus=res_mus,
             res_ps=res_ps, 
             res_gamma=res_gamma,
             res_V=res_V,
             icl=icl,
             zr=zr,
             probaV=probaV,
             string="clustM",
             m=m)
    
    return(res)
  }
