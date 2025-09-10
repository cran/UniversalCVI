\name{Hvalid}
\alias{Hvalid}
\title{
Wiroonsri(2024) correlation-based cluster validity indices and other well-known cluster validity indices
}
\description{
Computes the cluster validity indexes for a result of either kmeans or hierarchical clustering from user specified \code{kmin} to \code{kmax} used in Wiroonsri(2024). It includes the DI (J. C. Dunn, 1973) index, CH (T. Calinski and J. Harabasz, 1974) index, DB (D. L. Davies and D. W. Bouldin, 1979) index, PB (G. W. Miligan, 1985) index, CSL (C. H. Chou et al., 2004) index, PBM (M. K. Pakhira et al., 2004) index, DBs (M. Kim and R. S. Ramakrishna, 2005), Score function (S. Saitta et al., 2007), STR (A. Starczewski, 2017) index, NC, NCI, NCI1, and, NCI2 (N. Wiroonsri, 2024) indexes.
}
\usage{
Hvalid(x, kmax, kmin = 2, indexlist = "all", method = "kmeans",
  p = 2, q = 2, corr = "pearson", nstart = 100, sampling = 1, NCstart = TRUE)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{kmax}{a maximum number of clusters to be considered.
}
  \item{kmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{indexlist}{a character string indicating which cluster validity indexes to be computed (\code{"all"}, \code{"NC"}, \code{"NCI"}, \code{"NCI1"}, \code{"NCI2"}, \code{"PB"}, \code{"CSL"}, \code{"CH"}, \code{"DB"}, \code{"DBs"}, \code{"SF"}, \code{"DI"}, \code{"STR"}, \code{"PBM"}). More than one indexes can be selected.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"kmeans"}, \code{"hclust_complete"}, \code{"hclust_ward.D"}, \code{"hclust_ward.D2"}, \code{"hclust_average"}, \code{"hclust_single"}). The default is \code{"kmeans"}.
}
  \item{p}{the power of the Minkowski distance between centroids of clusters for \code{indexlist = c("DB","DBs")}. The default is \code{2}.
}
  \item{q}{the power of dispersion measure of a cluster for \code{indexlist = c("DB","DBs")}. The default is \code{2}.
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed (\code{"pearson"}, \code{"kendall"} or \code{"spearman"}). The default is \code{"pearson"}.
}
  \item{nstart}{a maximum number of initial random sets for kmeans for \code{method = "kmeans"}. The default is \code{100}.
}
  \item{sampling}{a number greater than 0 and less than or equal to 1 indicating the undersampling proportion of data to be used. This argument is intended for handling a large dataset. The default is \code{1}.}
  \item{NCstart}{logical for \code{indexlist} includes the \code{"NC"}, \code{"NCI"}, \code{"NCI1"}, and \code{"NCI2"}), if \code{TRUE}, the NC correlation at \code{k=1} is defined as the ratio introduced in the reference. Otherwise, it is assigned as \code{0}.
}
}
\details{
The well-known cluster validity indices used in Wiroonsri(2024). It includes the DI (J. C. Dunn, 1973) index, CH (T. Calinski and J. Harabasz, 1974) index, DB (D. L. Davies and D. W. Bouldin, 1979) index, PB (G. W. Miligan, 1980) index, CSL (C. H. Chou et al., 2004) index, PBM (M. K. Pakhira et al., 2004) index, DBs (M. Kim and R. S. Ramakrishna, 2005), Score function (S. Saitta et al., 2007), STR (A. Starczewski, 2017), NC, NCI, NCI1, and, NCI2 (N. Wiroonsri, 2024) indexes.

The NC correlation computes the correlation between an actual distance between a pair of data points and a centroid distance of clusters that the two points locate in. NCI1 and NCI2 are the proportion and the subtraction, respectively, of the same two ratios. The first ratio is the NC improvement from \code{k-1} clusters to \code{k} clusters over the entire room for improvement. The second ratio is the NC improvement from \code{k} clusters to \code{k+1} clusters over the entire room for improvement. NCI is a combination of NCI1 and NCI2.
}
\value{
  \item{NC}{the NC correlations for \code{k} from \code{kmin-1} to \code{kmax+1} shown in a data frame where the first and the second columns are \code{k} and the NC, respectively.}

Each of the followings shows the values of each index for \code{k} from \code{kmin} to \code{kmax} in a data frame.
  \item{NCI}{the NCI index.}
  \item{NCI1}{the NCI1 index.}
  \item{NCI2}{the NCI2 index.}
  \item{PB}{the PB index.}
  \item{DI}{the DI index.}
  \item{DB}{the DB index.}
  \item{DBs}{the DBs index.}
  \item{CSL}{the CSL index.}
  \item{CH}{the CH index.}
  \item{SF}{the Score function.}
  \item{STR}{the STR index.}
  \item{PBM}{the PBM index.}
}
\references{
J. C. Bezdek, N. R. Pal, "Some new indexes of cluster validity," \emph{IEEE Transactions on Systems, Man, and Cybernetics}, Part B, 28, 301-315 (1998).

T. Calinski, J. Harabasz, "A dendrite method for cluster analysis," \emph{Communications in Statistics}, 3, 1-27 (1974).

C. H. Chou, M. C. Su, E. Lai, "A new cluster validity measure and its application to image compression," \emph{Pattern Anal Applic}, 7, 205-220 (2004).

D. L. Davies, D. W. Bouldin, "A cluster separation measure," \emph{IEEE Trans Pattern Anal Machine Intell}, 1, 224-227 (1979).

J. C. Dunn, "A fuzzy relative of the ISODATA process and its use in detecting compact well-separated clusters," \emph{J Cybern}, 3(3), 32-57 (1973).

M. Kim, R. S. Ramakrishna, "New indices for cluster validity assessment," \emph{Pattern Recognition Letters}, 26, 2353-2363 (2005).

G. W. Miligan, "An examination of the effect of six types of error perturbation on fifteen clustering algorithms," \emph{Psychometrika}, 45, 325-342 (1980).

M. K. Pakhira, S. Bandyopadhyay and U. Maulik, "Validity index for crisp and fuzzy clusters," \emph{Pattern Recogn} 37(3):487–501 (2004).

S. Saitta, B. Raphael, I. Smith, "A bounded index for cluster validity," \emph{In Perner, P.: Machine Learning and Data Mining in Pattern Recognition, Lecture Notes in Computer Science}, 4571, Springer (2007).

A. Starczewski, "A new validity index for crisp clusters," \emph{Pattern Anal Applic} 20, 687–700 (2017).

N. Wiroonsri, "Clustering performance analysis using a new correlation based cluster validity index," Pattern Recognition, 145, 109910, 2024.

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

# Compute all the indices by Hvalid
Hvalid(scale(x), kmax = 15, kmin = 2, indexlist = "all",
  method = "kmeans", p = 2, q = 2, corr = "pearson", nstart = 100, NCstart = TRUE)

# Compute selected a set of indices ("NC","NCI","DI","DB")
Hvalid(scale(x), kmax = 15, kmin = 2, indexlist = c("NC","NCI","DI","DB"),
  method = "kmeans", p = 2, q = 2, corr = "pearson", nstart = 100, NCstart = TRUE)

# ---- Hierarchical ----

# Average linkage

# Compute all the indices by Hvalid
Hvalid(scale(x), kmax = 15, kmin = 2, indexlist = "all",
  method = "hclust_average", p = 2, q = 2, corr = "pearson", nstart = 100, NCstart = TRUE)

# Compute selected a set of indices ("NC","NCI","DI","DB")
Hvalid(scale(x), kmax = 15, kmin = 2, indexlist = c("NC","NCI","DI","DB"),
  method = "hclust_average", p = 2, q = 2, corr = "pearson", nstart = 100, NCstart = TRUE)

#---Plot and compare the indexes---

# Compute six cluster validity indexes of a kmeans clustering result for k from 2 to 15
IDX.list = c("NCI", "DI", "DB", "DBs", "CSL", "CH")

Hvalid.result = Hvalid(scale(x), kmax = 15, kmin = 2, indexlist = IDX.list,
  method = "hclust_average", p = 2, q = 2, corr = "pearson", nstart = 100, NCstart = TRUE)

# Plot the computed indexes
plot_idx(Hvalid.result)
}
