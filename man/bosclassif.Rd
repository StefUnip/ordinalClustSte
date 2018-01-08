\name{bosclassif}
\alias{bosclassif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
bosclassif
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs a classification on a dataset with features of the ordinal kind, and a label variable of the integer type (1,2,...,kr). The classification function proposes two classification models. The first one, (chosen by the option kc=0), is a multivariate BOS model assuming that, conditionally on the class of the observations, the feature are independent. The second model is a parsimonious version of the first model. Parcimony is introduced by grouping the features into clusters (as in co-clustering) and assuming that the features of a cluster have a common distribution. 
}
\usage{
bosclassif(x,y,to.predict,kr,kc=0,m,nbSEM=50,nbSEMburn=20,nbindmini=4,
  init='kmeans',disp=TRUE,iterordiEM=10) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Matrix of ordinal data. The missing values should be equal to 0.
}
\item{y}{
%%     ~~Describe \code{x} here~~
Vector of classes for each row of x. Must be labelled with integers (1,2,...,kr).
}
\item{to.predict}{
%%     ~~Describe \code{x} here~~
Matrix of ordinal data with same number of features than x. Represents the observation the user wants to give a label (1,...,kr).
}
  \item{kr}{
%%     ~~Describe \code{kr} here~~
Number of classes.
}
  \item{kc}{
%%     ~~Describe \code{kc} here~~
Set to 0 to choose a classical multivariate BOS model. Otherwise, o choose a parsimonious model, set an integer that indicates the number of column clusters.
}
  \item{m}{
%%     ~~Describe \code{m} here~~
Integer that defines the ordinal data's number of levels.
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

  \item{probas.to.predict}{Matrix with nrow=nrow(to.predict) and ncol=kr. For each observation of argument to.predict, indicates the probability to  belong to each class.}
  \item{zr.to.predict}{Vector of length nrow(to.predict) which indicates at which class the observations of to predict has been associated.}
  \item{V.to.predict}{Matrix with nrow=nrow(to.predict) and ncol=kr. Indicates at which class the observations of to.predict has been associated.}
  \item{xhat}{dataset with missing values completed.}
  \item{mus}{For each iteration, an array of dimension kr*kc that represents the mus that were estimated at each iteration.}
  \item{ps}{For each iteration, an array of dimension kr*kc that represents the pis that were estimated at each iteration.} 
  \item{gamma}{Array of gammas (row mixing proportions) that were estimated at each iteration.}
  \item{rho}{Array of rhos (columns mixing proportions) that were estimated at each iteration.}
  \item{W}{For each iteration, array of dimension J*H such that W[j,h]=1 if j belongs to cluster h.}
  \item{res_mus}{Array of dimension kr*kc that represents the resulting BOS position parameters mus.}
  \item{res_ps}{Array of dimension kr*kc that represents the resulting BOS precision parameters pis.}
  \item{res_gamma}{Vector with the resulting gamma.}
  \item{res_rho}{Vector with the resulting rho.}
  \item{res_V}{Array with the resulting V.}
  \item{res_W}{Array with the resulting W.}
  \item{icl}{ICL-BIC result value. Must be maximized.}
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

# number of classes to predict
kr <- 2
# different kc to test with cross-validation
kcol <- 1


res <- bosclassif(x=M.train,y=y.train,to.predict=M.validation,kr,kc=kcol,m=m,
                  nbSEM=nbSEM,nbSEMburn=nbSEMburn,
                  nbindmini=nbindmini, disp=TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"ordinal-data"}% use one of  RShowDoc("KEYWORDS")
\keyword{"co-clustering"}% __ONLY ONE__ keyword per line
