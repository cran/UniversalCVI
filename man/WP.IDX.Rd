\name{WP.IDX}
\alias{WP.IDX}
\title{
Wiroonsri and Preedasawakul (WP) index

}
\description{
Computes the WPC (WP correlation), WP, WPCI1 and WPCI2 (N. Wiroonsri and O. Preedasawakul, 2023) indexes for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax}.
}

\usage{
WP.IDX(x, cmax, cmin = 2, corr = 'pearson', method = 'FCM', fzm = 2,
  gamma = (fzm^2*7)/4, sampling = 1, iter = 100, nstart = 20, NCstart = TRUE)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{cmax}{a maximum number of clusters to be considered.
}
  \item{cmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed (\code{"pearson"}, \code{"kendall"} or \code{"spearman"}). The default is \code{"pearson"}.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"FCM"} or \code{"EM"}). The default is \code{"FCM"}.
}
  \item{fzm}{a number greater than 1 giving the degree of fuzzification for \code{method = "FCM"}. The default is \code{2}.
}
  \item{gamma}{adjusted fuzziness parameter for \code{indexlist} = (\code{"WP"}, \code{"WPC"}, \code{"WPCI1"}, \code{"WPCI2"}). The default is computed from \eqn{7fzm^2/4}.
}
  \item{sampling}{a number greater than 0 and less than or equal to 1 indicating the undersampling proportion of data to be used. This argument is intended for handling a large dataset. The default is \code{1}.}
  \item{iter}{a maximum number of iterations for \code{method = "FCM"}. The default is \code{100}.
}
  \item{nstart}{a maximum number of initial random sets for FCM for \code{method = "FCM"}. The default is \code{20}.
}
  \item{NCstart}{logical for \code{indexlist} = (\code{"WP"}, \code{"WPC"}, \code{"WPCI1"},\code{"WPCI2"}), if \code{TRUE}, the WP correlation at c=1 is defined as an adjusted sd of the distances between all data points and their mean. Otherwise, the WP correlation at c=1 is defined as 0.
}
}
\details{

  The newly introduced index was inspired by the recently introduced Wiroonsri index which is only compatible with hard clustering methods. \cr\cr

  The WPC computes the correlation between the actual distance between a pair of data points and the distance between adjusted centroids with respect to the pair. WPCI1 and WPCI2 are the proportion and the subtraction, respectively, of the same two ratios. The first ratio is the WPC improvement from \code{c-1} clusters to \code{c} clusters over the entire room for improvement. The second ratio is the WPC improvement from \code{c} clusters to \code{c+1} clusters over the entire room for improvement. \code{WP} is defined as a combination of \code{WPCI1} and \code{WPCI2}.\cr\cr
The largest value of \code{WP(c)} indicates a valid optimal partition.
}

\value{
  \item{WPC}{the WP correlations for \code{c} from \code{cmin-1} to \code{cmax+1} shown in a data frame where the first and the second columns are \code{c} and the WPC, respectively.}
  Each of the followings show the value of each index for \code{c} from \code{cmin} to \code{cmax} in a data frame.
  \item{WP}{the WP index.}
  \item{WPCI1}{the WPCI1 index.}
  \item{WPCI2}{the WPCI2 index.}
}
\references{
N. Wiroonsri, O. Preedasawakul, "A correlation-based fuzzy cluster validity index with secondary options detector," arXiv:2308.14785, 2023
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}
\seealso{
  \code{\link{R1_data}, \link{TANG.IDX}, \link{FzzyCVIs}, \link{WP.IDX}, \link{Hvalid}}
}
\examples{
library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]

# ---- FCM algorithm ----

# Compute all the indices by WP.IDX using default gamma
FCM.WP = WP.IDX(scale(x), cmax = 15, cmin = 2, corr = 'pearson', method = 'FCM', fzm = 2,
  iter = 100, nstart = 20, NCstart = TRUE)
print(FCM.WP$WP)

# The optimal number of cluster
FCM.WP$WP[which.max(FCM.WP$WP$WPI),]


# ---- EM algorithm ----

# Compute all the indices by WP.IDX using default gamma
EM.WP = WP.IDX(scale(x), cmax = 15, cmin = 2, corr = 'pearson', method = 'EM',
  iter = 100, nstart = 20, NCstart = TRUE)
print(EM.WP$WP)

# The optimal number of cluster
EM.WP$WP[which.max(EM.WP$WP$WPI),]
}
