\name{boscoclust}
\alias{boscoclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
boscoclust
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs a co-clustering on ordinal data by using the latent block model (cf references for further  details). A BOS distribution is used, and the parameters inference is realized with an SEM-Gbbs algorithm.
}
\usage{
boscoclust(x,kr,kc,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
  init='kmeans',disp=TRUE,iterordiEM=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Matrix of ordinal data. The missing values should be equal to 0.
}
  \item{kr}{
%%     ~~Describe \code{kr} here~~
Number of row clusters.
}
  \item{kc}{
%%     ~~Describe \code{kc} here~~
Number of column clusters.
}
  \item{m}{
%%     ~~Describe \code{m} here~~
Number of categories for ordinal variables.
}
  \item{nbSEM}{
%%     ~~Describe \code{nbSEM} here~~
Number of SEM-Gibbs iterations realized to estimate parameters.
}
  \item{nbSEMburn}{
%%     ~~Describe \code{nbSEMburn} here~~
Number of SEM-Gibbs burning iterations for estimating parameters. This parameter must be inferior to nbSEM.
}
  \item{nbindmini}{
%%     ~~Describe \code{nbindmini} here~~
Minimum number of cells belonging to a block.
}
  \item{init}{
%%     ~~Describe \code{init} here~~
String that indicates the kind of initialisation. Must one of these two words :  "kmeans" or "random".
}
  \item{disp}{
%%     ~~Describe \code{disp} here~~
Boolean that should be set to TRUE if the user wants the progress bars to be displayed in the console.
}
  \item{iterordiEM}{
%%     ~~Describe \code{iterordiEM} here~~
Number of iterations for the internal EM algorithm that estimates the parameters of BOS distribution.
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{xhat}{dataset completed with missing values completed.}
  \item{mus}{For each iteration, an array of dimension kr*kc that represents the mus that were estimated at each iteration.}
  \item{ps}{For each iteration, an array of dimension kr*kc that represents the pis that were estimated at each iteration.} 
  \item{gamma}{Array of gammas (row mixing proportions) that were estimated at each iteration.}
  \item{rho}{Array of rhos (columns mixing proportions) that were estimated at each iteration.}
  \item{V}{For each iteration, array of dimension N*G such that V[i,g]=1 if i belongs to cluster g}
  \item{W}{For each iteration, array of dimension J*H such that W[j,h]=1 if j belongs to cluster h.}
  \item{res_mus}{Array of dimension kr*kc that represents the resulting BOS position parameters mus.}
  \item{res_ps}{Array of dimension kr*kc that represents the resulting BOS precision parameters pis.}
  \item{res_gamma}{Vector with the resulting gamma.}
  \item{res_rho}{Vector with the resulting rho.}
  \item{res_V}{Array with the resulting V.}
  \item{res_W}{Array with the resulting W.}
  \item{icl}{ICL-BIC result value. Must be maximized.}
  \item{zr}{Vector with resulting row partitions.}
  \item{zc}{Vector with resulting column partitions.}
  \item{probaV}{Array representing the probability for each row to belong to each row-cluster.}
  \item{probaW}{Array representing the probability for each column to belong to each column-cluster.}
}
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

    kr=2
    kc=2
    m=3

    res <- boscoclust(x=Msimulated,kr = kr,kc=kc,m=m, nbSEM=nbSEM,
                      nbindmini=nbindmini,nbSEMburn=nbSEMburn, iterordiEM=5)
  }
  \donttest{
    # loading the real dataset
    data("dataqol")
    set.seed(5)

    # loading the ordinal data 
    M <- as.matrix(dataqol[,2:29])


    # defining different number of categories:
    m=4


    # defining number of row and column clusters
    krow = 5
    kcol = 4

    # configuration for the inference
    nbSEM=50
    nbSEMburn=40
    nbindmini=2

    # Co-clustering execution
    object <- boscoclust(x=M,kr=krow,kc=kcol,m=m,nbSEM=nbSEM,nbSEMburn=nbSEMburn,
                         nbindmini=nbindmini, init='kmeans',disp=TRUE)

  }
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"co-clustering"}% __ONLY ONE__ keyword per line
