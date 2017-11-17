\name{bosclustMulti}
\alias{bosclustMulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
bosclustMulti
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs a clustering on ordinal data by using the multiple latent block model (cf references for further  details). It allows the user to define D groups of variables that have different number of levels. A BOS distribution is used, and the parameters inference is realized with an SEM-Gbbs algorithm.
}
\usage{
bosclustMulti(x, kr, m, d.list, nbSEM = 50, nbSEMburn = 20, nbindmini = 4,
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
Vector of length D. Element d defines the number of levels for the variables of group d.
}
  \item{d.list}{
%%     ~~Describe \code{d.list} here~~
List of length D. Item d is a vector with the column indices of the variables of group d.
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
  \item{xhat}{List of length D. Item d represents the group d in dataset, with missing values completed.}
  \item{mus}{List of length D. Item d represents the mus that were estimated at each iteration for group d.}
  \item{ps}{List of length D. Item d represents the pis that were estimated at each iteration for group d.} 
  \item{gamma}{Vector of length kr indicating the row mixing proportions that were estimated at each iteration.}
  \item{V}{For each iteration, array of dimension N*G such that V[i,g]=1 if i belongs to cluster g.}
  \item{res_mus}{List of length D. Item d represents the resulting mus for group d.}
  \item{res_ps}{List of length D. Item d represents the resulting pis for group d.}
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
  \dontrun{
    # loading the real dataset
    data("dataqol")
    set.seed(5)

    # loading the ordinal data 
    M <- as.matrix(dataqol[,2:31])


    # defining different number of categories:
    m=c(4,7)


    # defining number of row and column clusters
    krow = 5

    # configuration for the inference
    nbSEM=70
    nbSEMburn=50
    nbindmini=1

    d.list <- list(1:28,29:30)

    # Co-clustering execution
    object <- bosclustMulti(x=M,kr=krow,m=m, d.list=d.list,
                        nbSEM=nbSEM,nbSEMburn=nbSEMburn,
                         nbindmini=nbindmini, init='kmeans',disp=TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"clustering"}% __ONLY ONE__ keyword per line
