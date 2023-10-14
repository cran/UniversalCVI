\name{TANG.IDX}
\alias{TANG.IDX}
\title{
Tang index
}
\description{
Computes the TANG (Y. Tang et al., 2005) index for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax}.
}
\usage{
TANG.IDX(x, cmax, cmin = 2, method = "FCM", fzm = 2, nstart = 20, iter = 100)
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
The Tang index  is defined as \cr
\deqn{TANG(c) = \frac{\sum_{j=1}^c \sum_{i=1}^n\mu_{ij}^2\| {x}_i-{v}_j\|^2 + \frac{1}{c(c-1)}\sum_{j\neq k}\| {v}_j-{v}_k\|^2}{\min_{j\neq k} \{ \| {v}_j-{v}_k\|^2 \}+\frac{1}{c}}.}
The smallest value of \eqn{TANG(c)} indicates a valid optimal partition.\cr
}
\value{
  \item{TANG}{the TANG index for \code{c} from \code{cmin} to \code{cmax} shown in a data frame where the first and the second columns are \code{c} and the TANG index, respectively.}
}
\references{
Y. Tang, F. Sun, and Z. Sun, “Improved validation index for fuzzy clustering,” in Proceedings of the 2005, American Control Conference, 2005., pp. 1120–1125 vol. 2, 2005. \url{https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1470111&isnumber=31519}
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

# Compute the TANG index
FCM.TANG = TANG.IDX(scale(x), cmax = 15, cmin = 2, method = "FCM",
  fzm = 2, nstart = 20, iter = 100)
print(FCM.TANG)

# The optimal number of cluster
FCM.TANG[which.min(FCM.TANG$TANG),]

# ---- EM algorithm ----

# Compute the TANG index
EM.TANG = TANG.IDX(scale(x), cmax = 15, cmin = 2, method = "EM",
  nstart = 20, iter = 100)
print(EM.TANG)

# The optimal number of cluster
EM.TANG[which.min(EM.TANG$TANG),]
}
