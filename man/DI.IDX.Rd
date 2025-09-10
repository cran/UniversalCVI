\name{DI.IDX}
\alias{DI.IDX}
\title{
Dunn index
}
\description{
Computes the DI (J. C. Dunn, 1973) index for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}
\usage{
DI.IDX(x, kmax, kmin = 2, method = "kmeans", nstart = 100)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{kmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"kmeans"}, \code{"hclust_ward.D"}, \code{"hclust_ward.D2"}, \code{"hclust_complete"}, \code{"hclust_average"}, \code{"hclust_single"}). The default is \code{"kmeans"}.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
}
\details{
The DI index is defined as \cr
\deqn{DI(k) = \min_{i \ne j \in [k]}\left\{\frac{\min\left\{d(x_u,x_v)|x_u\in C_i,x_v \in C_j\right\}}{\max_{l \in [k]}\max\left\{d(x_u,x_v)|x_u,x_v \in C_l\right\}}\right\}.}
The largest value of \eqn{DI(k)} indicates a valid optimal partition.
}
\value{
  \item{DI}{the DI index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the DI index, respectively.}
}
\references{
J. C. Dunn, "A fuzzy relative of the ISODATA process and its use in detecting compact well-separated clusters," \emph{J Cybern}, 3(3), 32-57 (1973).
}

\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{
\code{\link{Hvalid}, \link{Wvalid}, \link{DB.IDX}, \link{FzzyCVIs},  \link{R1_data}}
}

\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- Kmeans ----

# Compute the DI index
K.DI = DI.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans", nstart = 100)
print(K.DI)

# The optimal number of cluster
K.DI[which.max(K.DI$DI),]

# ---- Hierarchical ----

# Average linkage

# Compute the DI index
H.DI = DI.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average")
print(H.DI)

# The optimal number of cluster
H.DI[which.max(H.DI$DI),]
}
