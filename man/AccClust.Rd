\name{AccClust}
\alias{AccClust}
\title{Accuracy detection for a clustering result with known classes
}
\description{
Computes the accuracy of a clustering result of a dataset with known classes from the k-means, fuzzy c-means, or EM algorithm.
}
\usage{
AccClust(x, label.names = "label", algorithm = "FCM", fzm = 2,
  scale = TRUE, nstart = 100, iter = 100)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{label.names}{a character string indicating the true label column name. The default is \code{"label"}
}
  \item{algorithm}{a character string indicating which clustering methods to be used (\code{"FCM"}, \code{"EM"}, \code{"Kmeans"}). More than one methods may be selected. The default is \code{"FCM"}.
}
  \item{fzm}{a number greater than 1 giving the degree of fuzzification for \code{method = "FCM"}. The default is \code{2}.
}
  \item{scale}{logical, if \code{TRUE} (default), the dataset is normalized before clustering.}
  \item{nstart}{a maximum number of initial random sets for FCM for \code{method = "FCM" or "Kmeans" or c("Kmeans","FCM")}. The default is \code{20}.
}
  \item{iter}{a maximum number of iterations for \code{method = "FCM"}. The default is \code{100}.
}
}
\value{
  \item{kmeans}{Accuracy score from \code{0} to \code{1} of the k-means result}
  \item{FCM}{Accuracy score from \code{0} to \code{1} of the FCM result}
  \item{EM}{Accuracy score from \code{0} to \code{1} of the EM result}
}
\references{
N. Wiroonsri, O. Preedasawakul, A correlation-based fuzzy cluster validity index with secondary options detector, arXiv:2308.14785, 2023
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}
\seealso{\code{\link{R1_data}, \link{D1_data}, \link{FzzyCVIs}, \link{WP.IDX}, \link{XB.IDX}, \link{Hvalid}}}
\examples{

library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data

# Check accuracy of clustering results obtained by kmeans, FCM, and EM clustering
AccClust(x, label.names = "label",algorithm = c("Kmeans","FCM","EM"), fzm = 2,
  scale = TRUE, nstart = 20,iter = 100)

# Check accuracy of a clustering result obtained by the FCM algoritm
AccClust(x, label.names = "label",algorithm = "FCM", fzm = 2,
  scale = TRUE, nstart = 20,iter = 100)
}
