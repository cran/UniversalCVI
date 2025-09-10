\name{SF.IDX}
\alias{SF.IDX}
\title{The score function
}
\description{
Computes the SF (S. Saitta et al., 2007) index for a result either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax}.
}
\usage{
SF.IDX(x, kmax, kmin = 2, method = "kmeans", nstart = 100)
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
The smallest value of \eqn{SF(k)} indicates a valid optimal partition.
}
\value{
\item{SF}{the Score function index for \code{k} from \code{kmin} to \code{kmax} shown in a data frame where the first and the second columns are \code{k} and the SF index, respectively.}
}
\references{
S. Saitta, B. Raphael, I. Smith, "A bounded index for cluster validity," \emph{In Perner, P.: Machine Learning and Data Mining in Pattern Recognition, Lecture Notes in Computer Science}, 4571, Springer (2007).
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

# Compute the SF index
K.SF = SF.IDX(scale(x), kmax = 15, kmin = 2, method = "kmeans", nstart = 100)
print(K.SF)

# The optimal number of cluster
K.SF[which.min(K.SF$SF),]

# ---- Hierarchical ----

# Average linkage

# Compute the SF index
H.SF = SF.IDX(scale(x), kmax = 15, kmin = 2, method = "hclust_average")
print(H.SF)

# The optimal number of cluster
H.SF[which.min(H.SF$SF),]
}
