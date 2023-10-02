\name{KWON.IDX}
\alias{KWON.IDX}
\title{
KWON index
}
\description{
Computes the KWON (S. H. Kwon, 1998) index for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax}.
}
\usage{
KWON.IDX(x, cmax, cmin = 2, method = "FCM", fzm = 2, nstart = 20, iter = 100)
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
The KWON index is defined as \cr
\deqn{KWON(c) = \dfrac{\sum_{j=1}^c\sum_{i=1}^n \mu_{ij}^2  \|{x}_i-{v}_j\|^2 +\frac{1}{c}\sum_{j=1}^c\| {v}_j-{v}_0\|^2}{\min_{i \neq j} \| {v}_i-{v}_j\|^2}.
}
The smallest value of \eqn{KWON(c)} indicates a valid optimal partition.\cr

}
\value{
  \item{KWON}{the KWON index for \code{c} from \code{cmin} to \code{cmax} shown in a data frame where the first and the second columns are \code{c} and the KWON index, respectively.}
}
\references{
S. H. Kwon, “Cluster validity index for fuzzy clustering,” Electronics
letters, vol. 34, no. 22, pp. 2176–2177, 1998. \doi{10.1049/el:19981523}
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}
\seealso{\code{\link{R1_data}, \link{TANG.IDX}, \link{FzzyCVIs}, \link{WP.IDX}, \link{Hvalid}}}
\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- FCM algorithm ----

# Compute the KWON index
FCM.KWON = KWON.IDX(scale(x), cmax = 15, cmin = 2, method = "FCM",
  fzm = 2, nstart = 20, iter = 100)
print(FCM.KWON)
# The optimal number of cluster
FCM.KWON[which.min(FCM.KWON$KWON),]

# ---- EM algorithm ----

# Compute the KWON index
EM.KWON = KWON.IDX(scale(x), cmax = 15, cmin = 2, method = "EM",
  nstart = 20, iter = 100)
print(EM.KWON)
# The optimal number of cluster
EM.KWON[which.min(EM.KWON$KWON),]
}
