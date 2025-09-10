\name{STRPBM.IDX}
\alias{STRPBM.IDX}
\title{
Starczewski and Pakhira-Bandyopadhyay-Maulik for crisp clustering indexes
}
\description{
Computes the STR (A. Starczewski, 2017) and PBM (M. K. Pakhira et al., 2004) indexes for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}
\usage{
STRPBM.IDX(x, kmax, kmin = 2, method = "kmeans", indexlist = "all", nstart = 100)
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
  \item{indexlist}{a character string indicating which cluster validity indexes to be computed (\code{"all"}, \code{"STR"}, \code{"PBM"}). More than one indexes can be selected.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}

}
\details{
PBM index can be used with both crisp and fuzzy clustering algorithms.\cr
The largest value of \eqn{STR(k)} indicates a valid optimal partition.\cr
The largest value of \eqn{PBM(k)} indicates a valid optimal partition.
}
\value{
  \item{STR}{the STR index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the STR index, respectively.}
    \item{PBM}{the PBM index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the PBM index, respectively.}
}
\references{
M. K. Pakhira, S. Bandyopadhyay and U. Maulik, "Validity index for crisp and fuzzy clusters," \emph{Pattern Recogn} 37(3):487–501 (2004).

A. Starczewski, "A new validity index for crisp clusters," \emph{Pattern Anal Applic} 20, 687–700 (2017).

}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{
\code{\link{Wvalid}, \link{FzzyCVIs}, \link{DI.IDX},  \link{R1_data}}
}

\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- Kmeans ----

# Compute all the indices by STRPBM.IDX
K.ALL = STRPBM.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans",
  indexlist = "all", nstart = 100)
print(K.ALL)

# Compute STR index
K.STR = STRPBM.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans",
  indexlist = "STR", nstart = 100)
print(K.STR)

# ---- Hierarchical ----

# Average linkage

# Compute all the indices by STRPBM.IDX
H.ALL = STRPBM.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average",
  indexlist = "all")
print(H.ALL)

# Compute STR index
H.STR = STRPBM.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average",
  indexlist = "STR")
print(H.STR)

}

