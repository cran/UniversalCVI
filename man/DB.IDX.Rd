\name{DB.IDX}
\alias{DB.IDX}
\title{
Davies–Bouldin (DB) and DB* (DBs) indexes
}
\description{
Computes the DB (D. L. Davies and D. W. Bouldin, 1979) and DBs (M. Kim and R. S. Ramakrishna, 2005) indexes for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}
\usage{
DB.IDX(x, kmax, kmin = 2, method = "kmeans",
  indexlist = "all", p = 2, q = 2, nstart = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{kmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"kmeans"}, \code{"hclust_ward.D"}, \code{"hclust_ward.D2"}, \code{"hclust_complete"}, \code{"hclust_average"}, \code{"hclust_single"}). The default is \code{"kmeans"}.
}
  \item{indexlist}{a character string indicating which cluster validity indexes to be computed (\code{"all"}, \code{"DB"}, \code{"DBs"}). More than one indexes can be selected.
}
  \item{p}{the power of the Minkowski distance between centroids of clusters. The default is \code{2}.
}
  \item{q}{the power of dispersion measure of a cluster. The default is \code{2}.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
}
\details{
The lowest value of \eqn{DB(k),DBs(k)} indicates a valid optimal partition.
}
\value{
  \item{DB}{the DB index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the DB index, respectively.}
  \item{DBs}{the DBs index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the DBs index, respectively.}
}
\references{
D. L. Davies, D. W. Bouldin, "A cluster separation measure," \emph{IEEE Trans Pattern Anal Machine Intell}, 1, 224-227 (1979).

M. Kim, R. S. Ramakrishna, "New indices for cluster validity assessment," \emph{Pattern Recognition Letters}, 26, 2353-2363 (2005).
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

# Compute all the indices by DB.IDX
K.ALL = DB.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans",
  indexlist = "all", p = 2, q = 2, nstart = 100)
print(K.ALL)

# Compute DB index
K.DB = DB.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans",
  indexlist = "DB", p = 2, q = 2, nstart = 100)
print(K.DB)

# ---- Hierarchical ----

# Average linkage

# Compute all the indices by DB.IDX
H.ALL = DB.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average",
  indexlist = "all", p = 2, q = 2)
print(H.ALL)

# Compute DB index
H.DB = DB.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average",
  indexlist = "DB", p = 2, q = 2)
print(H.DB)
}
