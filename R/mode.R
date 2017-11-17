mode <-
function (x){
  # mode of a vector x
  return(as.integer(names(sort(-table(x)))[1]))
}
