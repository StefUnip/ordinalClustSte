\name{plot}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plots the result of a classification, clustering or co-clustering that were performed from the following functions:
bosclassif,bosclust,boscoclust.

}
\usage{
plot(object) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{x} here~~
  Result of one of the fuction listed above.
}
}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Margot Selosse, Julien Jacques, Christophe Biernacki.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
\dontshow{

    data("Msimulated")

    nbSEM=15
    nbSEMburn=10
    nbindmini=1
    init="kmeans"

    kr=2
    kc=2
    m=3

    res <- boscoclust(x=Msimulated,kr = kr,kc=kc,m=m, init=init,
                      nbSEM=nbSEM,nbSEMburn=nbSEMburn,nbindmini=nbindmini)

    plot(res)

  }
\donttest{

  data("dataqol")
  set.seed(5)

  # loading the ordinal data 
  M <- as.matrix(dataqol[,2:29])

  m=4

  krow = 4

  nbSEM=50
  nbSEMburn=40
  nbindmini=1

  object <- bosclust(x=M, kr=krow, m=m, nbSEM=nbSEM,
      nbindmini=nbindmini, nbSEMburn=nbSEMburn,
      init='kmeans',disp=TRUE)

  plot(object)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"co-clustering"}% __ONLY ONE__ keyword per line
