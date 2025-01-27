\name{KPBM.IDX}
\alias{KPBM.IDX}
\title{
Modified Kernel form of Pakhira-Bandyopadhyay-Maulik (KPBM) index
}
\description{
Computes the KPBM (C. Alok, 2010) index for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax}.
}
\usage{
KPBM.IDX(x, cmax, cmin = 2, method = "FCM", fzm = 2, nstart = 20, iter = 100)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{cmax}{a maximum number of clusters to be considered.
}
  \item{cmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"FCM"} or \code{"EM"}). The default is \code{"FCM"}.
}
  \item{fzm}{a number greater than 1 giving the degree of fuzzification for \code{method = "FCM"}. The default is \code{2}.
}
  \item{nstart}{a maximum number of initial random sets for FCM for \code{method = "FCM"}. The default is \code{20}.
}
  \item{iter}{a maximum number of iterations for \code{method = "FCM"}. The default is \code{100}.
}
}
\details{
The KPBM index is defined as \cr
\deqn{KPBM(c) = \left(\frac{\max_{j \neq k}\| {v}_j-{v}_k\|}{c\sum_{j=1}^c\sum_{i=1}^n\mu_{ij}\| {x}_i-{v}_j\|}\right)^2.}
The largest value of \eqn{KPBM(c)} indicates a valid optimal partition.\cr
}
\value{
  \item{KPBM}{the KPBM index for \code{c} from \code{cmin} to \code{cmax} shown in a data frame where the first and the second columns are \code{c} and the KPBM index, respectively.}
}
\references{
C. Alok. (2010). "An investigation of clustering algorithms and soft computing approaches for pattern recognition",  Department of Computer Science, Assam University.
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul}
\seealso{\code{\link{R1_data}, \link{TANG.IDX}, \link{FzzyCVIs}, \link{WP.IDX}, \link{Hvalid}}}
\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- FCM algorithm ----

# Compute the KPBM index
FCM.KPBM = KPBM.IDX(scale(x), cmax = 15, cmin = 2, method = "FCM",
  fzm = 2, nstart = 20, iter = 100)
print(FCM.KPBM)

# The optimal number of cluster
FCM.KPBM[which.max(FCM.KPBM$KPBM),]

# ---- EM algorithm ----

# Compute the KPBM index
EM.KPBM = KPBM.IDX(scale(x), cmax = 15, cmin = 2, method = "EM",
  nstart = 20, iter = 100)
print(EM.KPBM)

# The optimal number of cluster
EM.KPBM[which.max(EM.KPBM$KPBM),]
}
