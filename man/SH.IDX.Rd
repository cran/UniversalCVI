\name{SH.IDX}
\alias{SH.IDX}
\title{
Silhouette index
}
\description{
Computes the SH (Rousseeuw, 1987; Kaufman and Rousseeuw, 2009) index for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}

\usage{
SH.IDX(x, kmax, kmin = 2, method = "kmeans", nstart = 100)
}
%- maybe also 'usage' for other objects documented here.
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

For \eqn{i \in [n]}, \eqn{l \in [k]}, and \eqn{x_i \in C_l}, let

\deqn{a(i) = \dfrac{1}{|C_l|-1}\sum_{y \in C_l} \left\|x_i-y\right\| and}
\deqn{b(i) = \min_{r \neq l} \dfrac{1}{|C_r|} \sum_{y \in C_r} \left\|x_i-y\right\|.}
The silhouette value of one data point \eqn{x_j} is defined as:

\deqn{s(j) =
\begin{cases}
                     \dfrac{b(j) - a(j)}{\max\{a(j),b(i)\}}  &\text{  \ \ if \ }  |C_j| > 1  \\
                     0   &\text{  \ \ if \ }  |C_j| = 1
        \end{cases}.
        }

The silhouette index is defined as

\eqn{SH(k) = \dfrac{1}{n} \sum_{i = 1}^n s(i).}

The largest value of \eqn{SH(k)} indicates a valid optimal partition.
}
\value{
  \item{SH}{the SH index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the SH index, respectively.}
}

\references{
Rousseeuw, P.J., 1987. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. J. Comput. Appl. Math. 20, 53–65.

Kaufman, L. and Rousseeuw, P.J., 2009. Finding groups in data: an introduction to cluster analysis. John Wiley & Sons.

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

# ---- Hierarchical ----

# Average linkage

# Compute the SH index
H.SH = SH.IDX(scale(x), kmax = 10, kmin = 2, method = "hclust_average", nstart = 1)
print(H.SH)

# The optimal number of cluster
H.SH[which.max(H.SH$SH),]
}
