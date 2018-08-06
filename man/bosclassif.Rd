\name{bosclassif}
\alias{bosclassif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
bosclassif
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs a classification on a dataset with features of the ordinal kind, and a label variable of the integer type (1,2,...,kr). The classification function proposes two classification models. The first one, (chosen by the argument kc=0), is a multivariate BOS model assuming that, conditionally on the class of the observations, the feature are independent. The second model is a parsimonious version of the first model. Parcimony is introduced by grouping the features into clusters (as in co-clustering) and assuming that the features of a cluster have a common distribution. If the data contains ordinal features with D different numbers of levels, the data is going to be seen as D matrices of ordinal data.
}
\usage{
bosclassif(x, y, idx_list=c(0), kr, kc=0, init, nbSEM, nbSEMburn, nbindmini, m=0) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  %%     ~~Describe \code{x} here~~
  Matrix made of ordinal data, of dimension N*Jtot. The features with same numbers of levels must be placed side by side. The missing values should be coded as NA.
  }
  \item{y}{
  %%     ~~Describe \code{x} here~~
  Vector of length N. It should represent the classes corresponding to each row of x. Must be labelled with integers (1,2,...,kr).
  }
  \item{idx_list}{
  %%     ~~Describe \code{x} here~~
  Vector of length D. This argument is useful when variables have different numbers of levels. Element d should indicate where the variables with number of levels m[d] begins in matrix x.
  }
  \item{kr}{
  %%     ~~Describe \code{kr} here~~
  Number of row classes.
  }
  \item{kc}{
  %%     ~~Describe \code{kc} here~~
  Vector of length D. d^{th} element indicates the number of column clusters. Set to 0 to choose a classical multivariate BOS model.
  }
  \item{m}{
  %%     ~~Describe \code{m} here~~
  Vector of length D. d^{th} element defines the ordinal data's number of levels.
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
  String that indicates the kind of initialisation. Must be one of th following words :  "kmeans" or "random".
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%\item{comp1 }{Description of 'comp1'}
%%\item{comp2 }{Description of 'comp2'}
%% ...
  Return an object. The slots are:
  \item{@zr}{Vector of length N with resulting row partitions.}
  \item{@zc}{List of length D. d^{th} item is a vector of length J[d] representing the columns partitions for the group of variables d.}
  \item{@J}{Vector of length D. d^{th} item represents the number of columns for d^{th} group of variables.}
  \item{@W}{List of length D. Item d is a matrix of dimension J*kc[d] such that W[j,h]=1 if j belongs to cluster h.}
  \item{@V}{Matrix of dimension N*kr such that V[i,g]=1 if i belongs to cluster g.}
  \item{@icl}{ICL value for co-clustering.}
  \item{@kr}{Number of row classes.}
  \item{@name}{Name of the result.}
  \item{@number_distrib}{Number of groups of variables.}
  \item{@pi}{Vector of length kr. Row mixing proportions.}
  \item{@rho}{List of length D. d^{th} item represents the column mixing proportion for d^{th} group of variables.}
  \item{@dlist}{List of length d. d^{th} item represents the indexes of group of variables d.}
  \item{@kc}{Vector of length D. d^{th} element represents the number of clusters column H for d^{th} group of variables.}
  \item{@m}{Vector of length D. d^{th} element represents the number of levels of d^{th} group of variables.}
  \item{@nbSEM}{Number of SEM-Gibbs algorithm iteration.}
  \item{@params}{List of length D. d^{th} item represents the blocks paramaters for group of variables d.}
  \item{@xhat}{List of length D. d^{th} item represents the d^{th} group of variables dataset, with missing values completed.}
  
  





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
# loading the real dataset
data("dataqol.classif")

set.seed(5)

# loading the ordinal data
M <- as.matrix(dataqol.classif[,2:29])


# creating the classes values
y <- as.vector(dataqol.classif$death)


# sampling datasets for training and to predict
nb.sample <- ceiling(nrow(M)*2/3)
sample.train <- sample(1:nrow(M), nb.sample, replace=FALSE)

M.train <- M[sample.train,]
M.validation <- M[-sample.train,]
nb.missing.validation <- length(which(M.validation==0))
m <- c(4)
M.validation[which(M.validation==0)] <- sample(1:m, nb.missing.validation,replace=TRUE)


y.train <- y[sample.train]
y.validation <- y[-sample.train]



# configuration for SEM algorithm
nbSEM=50
nbSEMburn=40
nbindmini=1
init="kmeans"

# number of classes to predict
kr <- 2
# different kc to test with cross-validation
kcol <- 1


res <- bosclassif(x=M.train,y=y.train,kr=kr,kc=kcol,m=m,
                  nbSEM=nbSEM,nbSEMburn=nbSEMburn,
                  nbindmini=nbindmini,init=init)

predictions <- predict(res, M.validation)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"co-clustering"}% __ONLY ONE__ keyword per line
