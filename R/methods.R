setMethod(
  f="predict",
  signature = "ResultClassifOrdinal",
  definition = function(object, x=matrix(0,nrow=1,ncol=1)) {
    res <- predictions(object, x)
    return(res)
  }
)


setGeneric("plot",function(object){standardGeneric("plot")}) 

setMethod(
  f="plot",
  signature = c("ResultClassifOrdinal"),
  definition = function(object) {
    res <- bosplot(object)
    return(res)
  }
)

setMethod(
  f="plot",
  signature = c("ResultCoclustOrdinal"),
  definition = function(object) {
    res <- bosplot(object)
    return(res)
  }
)

setMethod(
  f="plot",
  signature = c("ResultClustOrdinal"),
  definition = function(object) {
    res <- bosplot(object)
    return(res)
  }
)