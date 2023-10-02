\name{Wvalid}
\alias{Wvalid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wiroonsri(2024) correlation-based cluster validity indices
}
\description{
Computes the NC correlation, NCI, NCI1 and NCI2 cluster validity indices for the number of clusters from user specified \code{kmin} to \code{kmax} obtained from either K-means or hierarchical clustering based on the recent paper by Wiroonsri(2024).
}
\usage{
Wvalid(x, kmax, kmin = 2, method = "kmeans",
  corr = "pearson", nstart = 100, NCstart = TRUE)
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
  \item{corr}{a character string indicating which correlation coefficient is to be computed (\code{"pearson"}, \code{"kendall"} or \code{"spearman"}). The default is \code{"pearson"}.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
  \item{NCstart}{logical for \code{indexlist} includes the \code{"NC"}, \code{"NCI"}, \code{"NCI1"}, and \code{"NCI2"}), if \code{TRUE}, the NC correlation at \code{k=1} is defined as the ratio introduced in the reference. Otherwise, it is assigned as \code{0}.
}
}
\details{
The NC correlation computes the correlation between an actual distance between a pair of data points and a centroid distance of clusters that the two points locate in. NCI1 and NCI2 are the proportion and the subtraction, respectively, of the same two ratios. The first ratio is the NC improvement from \code{k-1} clusters to \code{k} clusters over the entire room for improvement. The second ratio is the NC improvement from \code{k} clusters to \code{k+1} clusters over the entire room for improvement. NCI is a combination of NCI1 and NCI2.
}
\value{
  \item{NC}{the NC correlations for \code{k} from \code{kmin-1} to \code{kmax+1} shown in a data frame where the first and the second columns are \code{k} and the NC, respectively.}
  Each of the followings shows the values of each index for \code{k} from \code{kmin} to \code{kmax} in a data frame.
  \item{NCI}{the NCI index.}
  \item{NCI1}{the NCI1 index.}
  \item{NCI2}{the NCI2 index.}
  }
\references{
N. Wiroonsri, "Clustering performance analysis using a new correlation based cluster validity index," Pattern Recognition, 145, 109910, 2024. \doi{10.1016/j.patcog.2023.109910}
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}
\seealso{
\code{\link{Hvalid}, \link{FzzyCVIs}, \link{DB.IDX},  \link{R1_data}}
}
\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- Kmeans ----

# Compute all the indices by Wvalid
K.NC = Wvalid(scale(x), kmax = 15, kmin=2, method = 'kmeans',
  corr='pearson', nstart=100, NCstart = TRUE)
print(K.NC)

# The optimal number of cluster
K.NC$NCI[which.max(K.NC$NCI$NCI),]

# ---- Hierarchical ----

# Average linkage

# Compute all the indices by Wvalid
H.NC = Wvalid(scale(x), kmax = 15, kmin=2, method = 'hclust_average',
  corr='pearson', nstart=100, NCstart = TRUE)
print(H.NC)

# The optimal number of cluster
H.NC$NCI[which.max(H.NC$NCI$NCI),]
}

