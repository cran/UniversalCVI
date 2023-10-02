\name{GC.IDX}
\alias{GC.IDX}
\title{The generalized C index
}
\description{
Computes the GC1 GC2 GC3 and GC4 (J. C. Bezdek et al., 2016) indexes for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax}.
}
\usage{
GC.IDX(x, cmax, cmin = 2, indexlist = "all", method = 'FCM', fzm = 2,
  iter = 100, nstart = 20)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{cmax}{a maximum number of clusters to be considered.
}
  \item{cmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{indexlist}{a character string indicating which The generalized C index be computed ("\code{all}","\code{GC1}","\code{GC2}","\code{GC3}","\code{GC4}"). More than one indexes can be selected.
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
The GC index is a soft version of the C-index, formulated based on relational transformations of the membership degree matrix \eqn{\mu}. It comprises four distinct variants, each with its own definition.\cr
The smallest value of \eqn{GC(c)} indicates a valid optimal partition.
}
\value{
  Each of the followings shows the values of each index for \code{c} from \code{cmin} to \code{cmax} in a data frame.
  \item{GC1}{the generalized C index (\eqn{\sum\cdot \sim} Sum-Product).}
  \item{GC2}{the generalized C index (\eqn{\sum\wedge \sim} Sum-Min).}
  \item{GC3}{the generalized C index (\eqn{\vee\cdot \sim} Max-Product).}
  \item{GC4}{the generalized C index (\eqn{\vee\wedge \sim} Max-Min).}
}
\references{
J. C. Bezdek, M. Moshtaghi, T. Runkler, and C. Leckie, “The generalized
c index for internal fuzzy cluster validity,” IEEE Transactions on Fuzzy
Systems, vol. 24, no. 6, pp. 1500–1512, 2016. \url{https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7429723&isnumber=7797168}
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

# Compute all the indices by GC.IDX
FCM.all.GC = GC.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "all",
  method =  'FCM', fzm = 2, iter = 100, nstart = 5)
print(FCM.all.GC)

# Compute GC2 index
FCM.GC2 = GC.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "GC2",
  method = 'FCM', fzm = 2, iter = 100, nstart = 5)
print(FCM.GC2)

# ---- EM algorithm ----

# Compute all the indices by GC.IDX
EM.all.GC = GC.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "all",
  method =  'EM', iter = 100, nstart = 5)
print(EM.all.GC)

# Compute GC2 index
EM.GC2 = GC.IDX(scale(x), cmax = 10, cmin = 2, indexlist = "GC2",
  method =  'EM', iter = 100, nstart = 5)
print(EM.GC2)
}
