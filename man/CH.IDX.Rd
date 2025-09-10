\name{CH.IDX}
\alias{CH.IDX}
\title{
Calinski–Harabasz (CH) index
}
\description{
Computes the CH (T. Calinski and J. Harabasz, 1974) index for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.

}
\usage{
CH.IDX(x, kmax, kmin = 2, method = "kmeans", nstart = 100)
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
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
}
\details{
The CH index is defined as \cr
\deqn{CH(k) = \frac{n-k}{k-1}\frac{\sum_{i=1}^k|C|_id(v_i,\bar{x})}{\sum_{i=1}^k\sum_{x_j\in C_i}d(x_j,v_i)}}

The largest value of \eqn{CH(k)} indicates a valid optimal partition.
}
\value{
  \item{CH}{the CH index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the CH index, respectively.}
}

\references{
T. Calinski, J. Harabasz, "A dendrite method for cluster analysis," \emph{Communications in Statistics}, 3, 1-27 (1974).}

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

# Compute the CH index
K.CH = CH.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans", nstart = 100)
print(K.CH)

# The optimal number of cluster
K.CH[which.max(K.CH$CH),]

# ---- Hierarchical ----

# Average linkage

# Compute the CH index
H.CH = CH.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average")
print(H.CH)

# The optimal number of cluster
H.CH[which.max(H.CH$CH),]
}
