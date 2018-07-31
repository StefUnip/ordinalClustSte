setClass (
  "xhat",
  representation=representation(
    xhat="matrix"
  )
)
setClass (
  "W",
  representation=representation(
    W="matrix"
  )
)
setClass (
  "zc",
  representation=representation(
    zc="vector"
  )
)
setClass (
  "zcchain",
  representation=representation(
    zcchain="list"
  )
)
setClass (
  "zrchain",
  representation=representation(
    zrchain="vector"
  )
)
setClass (
  "rho",
  representation=representation(
    rho="vector"
  )
)
setClass (
  "rhochain",
  representation=representation(
    rhochain="list"
  )
)
setClass (
  "pichain",
  representation=representation(
    pichain="vector"
  )
)
setClass (
  "params",
  representation=representation(
    params="list"
  )
)
setClass (
  "paramschain",
  representation=representation(
    paramschain="list"
  )
)
setClass (
  "dlist",
  representation=representation(
    dlist="vector"
  )
)



# Result for co-clustering
setClass (
  "ResultCoclustOrdinal",
  
  # Defining slot type
  representation = representation (
    V = "matrix",
    zr = "vector",
    pi = "vector",
    m = "vector",
    icl = "numeric",
    name="character"
  ),
  contains=c("W","zc","rho","params","paramschain","xhat","pichain","rhochain","zrchain","zcchain")
)


# Result for clustering
setClass (
  "ResultClustOrdinal",
  
  # Defining slot type
  representation = representation (
    V = "matrix",
    zr = "vector",
    pi = "vector",
    m = "vector",
    icl = "numeric",
    name="character"
  ),
  contains=c("params","paramschain","xhat","zrchain","pichain")
)



# Result for classification
setClass (
  "ResultClassifOrdinal",
  
  # Defining slot type
  representation = representation (
    V = "matrix",
    zr = "vector",
    pi = "vector",
    icl = "numeric",
    kr = "integer",
    kc = "vector",
    J = "vector",
    number_distrib = "integer",
    m = "vector",
    nbSEM = "integer",
    name="character"
  ),
  contains=c("W","zc","rho","params","dlist","xhat")
)

setClass (
  "ResultPredictionOrdinal",
  
  # Defining slot type
  representation = representation (
    zr_topredict = "vector",
    V_topredict = "matrix"
  )
)