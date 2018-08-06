methods::setMethod(
  f="predict",
  signature = "ResultClassifOrdinal",
  definition = function(object, x=matrix(0,nrow=1,ncol=1)) {
    res <- predictions(object, x)
    return(res)
  }
)


methods::setGeneric("plot",function(object){standardGeneric("plot")}) 

methods::setMethod(
  f="plot",
  signature = c("ResultClassifOrdinal"),
  definition = function(object) {
    res <- bosplot(object)
    return(res)
  }
)

methods::setMethod(
  f="plot",
  signature = c("ResultCoclustOrdinal"),
  definition = function(object) {
    res <- bosplot(object)
    return(res)
  }
)

methods::setMethod(
  f="plot",
  signature = c("ResultClustOrdinal"),
  definition = function(object) {
    res <- bosplot(object)
    return(res)
  }
)