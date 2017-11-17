rdirichlet <-
function (n, alpha){
  # fonction rdirichlet du packahe gtools
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}
