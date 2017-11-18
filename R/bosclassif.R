bosclassif <-
  function (x,y,to.predict,kr,kc=0,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
            init='kmeans',disp=TRUE,iterordiEM=10) {
    
    # Defining if parsimonious or not
    parsimonious = FALSE
    if(kc==0){
      res = bosclassif.no.kc(x=x,y=y,to.predict=to.predict,kr=kr,m=m,nbSEM=nbSEM,nbSEMburn=nbSEMburn,nbindmini=nbindmini,
            init=init,disp=disp,iterordiEM=iterordiEM)
    }
    else{
      res = bosclassif.kc(x=x,y=y,to.predict=to.predict,kr=kr,kc=kc,m=m,nbSEM=nbSEM,nbSEMburn=nbSEMburn,nbindmini=nbindmini,
            init=init,disp=disp,iterordiEM=iterordiEM)
    }


    return(res) 
  }
