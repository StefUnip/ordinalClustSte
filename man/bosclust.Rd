\name{bosclust}
\alias{bosclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
bosclust
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs a clustering on ordinal data by using the multiple latent block model (cf references for further  details). It allows the user to define D groups of variables that have different number of levels. A BOS distribution is used, and the parameters inference is realized with an SEM-Gbbs algorithm.
}
\usage{
bosclust(x, kr, m, nbSEM = 50, nbSEMburn = 20, nbindmini = 4, 
  init = "kmeans", disp = TRUE, iterordiEM = 10)
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
  \item{m}{
%%     ~~Describe \code{m} here~~
Defines the number of levels for the variables.
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
  \item{xhat}{Represents the dataset, with missing values completed.}
  \item{mu}{Represents the BOS position parameter mus that were estimated at each iteration of the SEM algorithm.}
  \item{p}{Represents the BOS precision parameter pis that were estimated at each iteration of the SE algorithm.} 
  \item{gamma}{Vector of length kr indicating the row mixing proportions that were estimated at each iteration of the SE algorithm.}
  \item{V}{For each iteration of the SEM algorithm, array of dimension N*G such that V[i,g]=1 if i belongs to cluster g.}
  \item{res_mu}{Represents the resulting position parameter mu.}
  \item{res_p}{Represents the resulting precision parameter pi.}
  \item{res_gamma}{Vector with the resulting gammas.}
  \item{res_V}{Array with the resulting V.}
  \item{icl}{ICL-BIC result value. Must be maximized.}
  \item{zr}{Vector with resulting row partitions.}
  \item{probaV}{Array representing the probability for each row to belong to each row-cluster.}
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
  #\dontshow{
  #  data("Msimulated")

  #  nbSEM=15
  #  nbSEMburn=10
  #  nbindmini=1

  #  kr=2
  #  m=3

  #  res <- bosclust(x=Msimulated,kr = kr,m=m,
  #                  nbSEM=nbSEM,nbSEMburn=nbSEMburn,
  #                  nbindmini=#nbindmini,
  #                  disp=FALSE,iterordiEM=5)
  #}
  \donttest{
    library(ordinalClust)
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
  }
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"clustering"}% __ONLY ONE__ keyword per line
