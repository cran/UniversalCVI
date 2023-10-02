\name{PB.IDX}
\alias{PB.IDX}
\title{
Point biserial correlation (PB)
}
\description{
Computes the PB (G. W. Miligan, 1980) index for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}
\usage{
PB.IDX(x, kmax, kmin = 2, method = "kmeans", corr = "pearson", nstart = 100)
}
\arguments{
 \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{kmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"kmeans"}, \code{"hclust_complete"}, \code{"hclust_average"}, \code{"hclust_single"}). The default is \code{"kmeans"}.
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed (\code{"pearson"}, \code{"kendall"} or \code{"spearman"}). The default is \code{"pearson"}.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
}
\details{
The largest value of \eqn{PB(k)} indicates a valid optimal partition.
}
\value{
  \item{PB}{the PB index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the PB index, respectively.}
}
\references{
G. W. Miligan, "An examination of the effect of six types of error perturbation on fifteen clustering algorithms," \emph{Psychometrika}, 45, 325-342 (1980).
}

\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{
\code{\link{Hvalid}, \link{Wvalid}, \link{DI.IDX}, \link{FzzyCVIs},  \link{R1_data}}
}

\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- Kmeans ----

# Compute PB index
K.PB = PB.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans",
  corr = "pearson", nstart = 100)
print(K.PB)

# The optimal number of cluster
K.PB[which.max(K.PB$PB),]

# ---- Hierarchical ----

# Average linkage

# Compute PB index
H.PB = PB.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average",
  corr = "pearson")
print(H.PB)

# The optimal number of cluster
H.PB[which.max(H.PB$PB),]
}
