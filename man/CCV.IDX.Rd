\name{CCV.IDX}
\alias{CCV.IDX}
\title{Correlation Cluster Validity (CCV) index
}
\description{
Computes the CCVP and CCVS (M. Popescu et al., 2013) indexes for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax}.
}
\usage{
CCV.IDX(x, cmax, cmin = 2, indexlist = "all", method = 'FCM', fzm = 2,
  iter = 100, nstart = 20)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{cmax}{a maximum number of clusters to be considered.
}
  \item{cmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{indexlist}{a character string indicating which The generalized C index be computed ("\code{all}","\code{CCVP}","\code{CCVS}"). More than one indexes can be selected.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"FCM"} or \code{"EM"}). The default is \code{"FCM"}.
}
  \item{fzm}{a number greater than 1 giving the degree of fuzzification for \code{method = "FCM"}. The default is \code{2}.
}
  \item{iter}{a maximum number of iterations for \code{method = "FCM"}. The default is \code{100}.
}
  \item{nstart}{a maximum number of initial random sets for FCM for \code{method = "FCM"}. The default is \code{20}.
}
}
\details{
A new cluster validity framework that compares the structure in the data to the structure of dissimilarity matrices induced by a matrix transformation of the partition being tested.  The largest value of \eqn{CCV(c)} indicates a valid optimal partition.
}
\value{
  Each of the followings shows the values of each index for \code{c} from \code{cmin} to \code{cmax} in a data frame.
  \item{CCVP}{the Pearson Correlation Cluster Validity index.}
  \item{CCVS}{the Spearman’s (rho) Correlation Cluster Validity index.}
}
\references{
M. Popescu, J. C. Bezdek, T. C. Havens and J. M. Keller (2013). "A Cluster Validity Framework Based on Induced Partition Dissimilarity." \url{https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6246717&isnumber=6340245}
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}
\seealso{
\code{\link{R1_data}, \link{TANG.IDX}, \link{FzzyCVIs}, \link{WP.IDX}, \link{Hvalid}}
}
\examples{

library(UniversalCVI)

# Iris data
x = iris[,1:4]

# ---- FCM algorithm ----


# Compute all the indices by CCV.IDX
FCM.ALL.CCV = CCV.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "all",
  method =  'FCM', fzm = 2, iter = 100, nstart = 20)
print(FCM.ALL.CCV)

# Compute CCVP index
FCM.CCVP = CCV.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "CCVP",
  method =  'FCM', fzm = 2, iter = 100, nstart = 20)
print(FCM.CCVP)


# ---- EM algorithm ----

# Compute all the indices by CCV.IDX
EM.ALL.CCV = CCV.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "all",
  method =  'EM', iter = 100, nstart = 20)
print(EM.ALL.CCV)

# Compute CCVP index
EM.CCVP = CCV.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "CCVP",
  method =  'EM', iter = 100, nstart = 20)
print(EM.CCVP)
}
