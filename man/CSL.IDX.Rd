\name{CSL.IDX}
\alias{CSL.IDX}
\title{
Chou-Su-Lai (CSL) index
}
\description{
Computes the CSL (C. H. Chou et al., 2004) index for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}
\usage{
CSL.IDX(x, kmax, kmin = 2, method = "kmeans", nstart = 100)
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
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
}
\details{
The CSL index is defined as \cr
\deqn{CSL(k) = \frac{\sum_{i=1}^k \left\{\frac{1}{|C_i|}\sum_{x_j \in C_i} \max_{x_l \in C_i} d(x_j,x_l)\right\}}{\sum_{i=1}^k \left\{\min_{j:j \ne i}d(v_i,v_j)\right\}}.
   }
The smallest value of \eqn{CSL(k)} indicates a valid optimal partition.
}
\value{
 \item{CSL}{the CSL index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the CSL index, respectively.}
}
\references{
C. H. Chou, M. C. Su, E. Lai, "A new cluster validity measure and its application to image compression," \emph{Pattern Anal Applic}, 7, 205-220 (2004).
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

# Compute the CSL index
K.CSL = CSL.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans", nstart = 100)
print(K.CSL)

# The optimal number of cluster
K.CSL[which.min(K.CSL$CSL),]

# ---- Hierarchical ----

# Average linkage

# Compute the CSL index
H.CSL = CSL.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average")
print(H.CSL)

# The optimal number of cluster
H.CSL[which.min(H.CSL$CSL),]
}
