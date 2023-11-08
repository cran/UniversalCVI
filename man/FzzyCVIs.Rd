\name{FzzyCVIs}
\alias{FzzyCVIs}
\title{
Fuzzy cluster validity indexes used in Wiroonsri and Preedasawakul (2023)
}
\description{Computes the cluster validity indexes for a result of either FCM or EM clustering from user specified \code{cmin} to \code{cmax} used in Wiroonsri and Preedasawakul (2023). It includes the XB (X. L. Xie and G. Beni, 1991) index, KWON (S. H. Kwon, 1998) index, KWON2 (S. H. Kwon et al., 2021) index, TANG (Y. Tang et al., 2005) index , HF (F. Haouas et al., 2017) index, WL (C. H. Wu et al., 2015) index, PBM (M. K. Pakhira et al., 2004) index, KPBM (C. Alok, 2010) index, CCVP and CCVS (M. Popescu et al., 2013) index, GC1, GC2, GC3, and GC4 (J. C. Bezdek et al., 2016) indexes , WPC, WP, WPCI1, and, WPCI2 (N. Wiroonsri and O. Preedasawakul, 2023) indexes.
}
\usage{
FzzyCVIs(x, cmax, cmin = 2, indexlist = 'all', corr = 'pearson',
  method = 'FCM', fzm = 2, gamma = (fzm^2*7)/4, sampling = 1,
  iter = 100, nstart = 20, NCstart = TRUE)
}
\arguments{
  \item{x}{a numeric data frame or matrix where each column is a variable to be used for cluster analysis and each row is a data point.
}
  \item{cmax}{a maximum number of clusters to be considered.
}
  \item{cmin}{a minimum number of clusters to be considered. The default is \code{2}.
}
  \item{indexlist}{a character string indicating which cluster validity indexes to be computed (\code{"all"}, \code{"WPC"}, \code{"WP"}, \code{"WPCI1"}, \code{"WPCI2"}, \code{"XB"}, "\code{"KWON"}", "\code{"KWON2"}", "\code{"TANG"}", "\code{"HF"}", \code{"WL"}, \code{"PBM"}, \code{"KPBM"}, \code{"CCVP"}, \code{"CCVS"}, \code{"GC1"}, \code{"GC2"}, \code{"GC3"}, \code{"GC4"}). More than one indexes can be selected.
}
  \item{corr}{a character string indicating which correlation coefficient is to be computed (\code{"pearson"}, \code{"kendall"} or \code{"spearman"}) for \code{indexlist} = (\code{"WP"}, \code{"WPC"}, \code{"WPCI1"},\code{"WPCI2"}, \code{"CCVP"}, \code{"CCVS"}, \code{"GC1"}, \code{"GC2"}, \code{"GC3"} or \code{"GC4"}). The default is \code{"pearson"}.
}
  \item{method}{a character string indicating which clustering method to be used (\code{"FCM"} or \code{"EM"}). The default is \code{"FCM"}.
}
  \item{fzm}{a number greater than 1 giving the degree of fuzzification for \code{method = "FCM"}. The default is \code{2}.
}
  \item{gamma}{adjusted fuzziness parameter for \code{indexlist} = (\code{"WP"}, \code{"WPC"}, \code{"WPCI1"}, \code{"WPCI2"}). The default is \eqn{7fzm^2/4}.
}
  \item{sampling}{a number greater than 0 and less than or equal to 1 indicating the undersampling proportion of data to be used. This argument is intended for handling a large dataset. The default is \code{1}.}
  \item{iter}{a maximum number of iterations for \code{method = "FCM"}. The default is \code{100}.
}
  \item{nstart}{a maximum number of initial random sets for FCM for \code{method = "FCM"}. The default is \code{20}.
}
  \item{NCstart}{logical for \code{indexlist} includes either of the \code{"WP"}, \code{"WPC"}, \code{"WPCI1"}, and \code{"WPCI2"}), if \code{TRUE}, the WP correlation at \code{c=1} is defined as the ratio introduced in the reference. Otherwise, it is assigned as \code{0}.
}
}
\details{

  The well-known cluster validity indexes for either FCM or EM clustering. It includes the XB (X. L. Xie and G. Beni., 1991) index, KWON (S. H. Kwon, 1998) index, KWON2 (S. H. Kwon et al., 2021) index, TANG (Y. Tang et al., 2005) index , HF (F. Haouas et al., 2017) index, WL (C. H. Wu et al., 2015) index, PBM (M. K. Pakhira et al., 2004) index, KPBM (C. Alok, 2010) index, CCVP and CCVS (M. Popescu et al., 2013) index, GC1, GC2, GC3, and GC4 (J. C. Bezdek et al., 2016) indexes , WPC, WP, WPCI1, and, WPCI2 (N. Wiroonsri and O. Preedasawakul, 2023) indexes.

The WPC computes the correlation between the actual distance between a pair of data points and the distance between adjusted centroids with respect to the pair. WPCI1 and WPCI2 are the proportion and the subtraction, respectively, of the same two ratios. The first ratio is the WPC improvement from \code{c-1} clusters to \code{c} clusters over the entire room for improvement. The second ratio is the WPC improvement from \code{c} clusters to \code{c+1} clusters over the entire room for improvement. \code{WP} is defined as a combination of \code{WPCI1} and \code{WPCI2}.
  }

\value{
  \item{WPC}{the WP correlation from \code{c} from \code{cmin-1} to \code{cmax+1} shown in a data frame.}
  Each of the followings shows the values of each index for \code{c} from \code{cmin} to \code{cmax} in a data frame.
  \item{WP}{the WP index.}
  \item{WPCI1}{the WPCI1 index.}
  \item{WPCI2}{the WPCI2 index.}
  \item{XB}{the XB index.}
  \item{KWON}{the KWON index.}
  \item{KWON2}{the KWON2 index.}
  \item{TANG}{the TANG index.}
  \item{HF}{the HF index.}
  \item{WL}{the WL index.}
  \item{PBM}{the PBM index}
  \item{KPBM}{the KPBM index}
  \item{CCVP}{the Pearson Correlation Cluster Validity index.}
  \item{CCVS}{the Spearman’s (rho) Correlation Cluster Validity index.}
  \item{GC1}{the generalized C index (\eqn{\sum\cdot \sim} Sum-Product).}
  \item{GC2}{the generalized C index (\eqn{\sum\wedge \sim} Sum-Min).}
  \item{GC3}{the generalized C index (\eqn{\vee\cdot \sim} Max-Product).}
  \item{GC4}{the generalized C index (\eqn{\vee\wedge \sim} Max-Min).}
}

\references{

C. Alok. (2010). "An investigation of clustering algorithms and soft computing approaches for pattern recognition," Department of Computer Science, Assam University.\cr\cr

J. C. Bezdek, M. Moshtaghi, T. Runkler, C. Leckie, “The generalized
c index for internal fuzzy cluster validity,” IEEE Transactions on Fuzzy
Systems, vol. 24, no. 6, pp. 1500–1512, 2016. \cr\cr

F. Haouas, Z. Ben Dhiaf, A. Hammouda, B. Solaiman, "A new efficient fuzzy cluster validity index: Application to images clustering," 2017 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE), Naples, Italy, 2017, pp. 1-6.\cr\cr

S. H. Kwon, “Cluster validity index for fuzzy clustering,” Electronics
letters, vol. 34, no. 22, pp. 2176–2177, 1998. \cr\cr

S. H. Kwon, J. Kim, S. H. Son, “Improved cluster validity index
for fuzzy clustering,” Electronics Letters, vol. 57, no. 21, pp. 792–794,
2021. \cr\cr

M. K. Pakhira, S. Bandyopadhyay, U. Maulik, “Validity index for crisp and fuzzy clusters,” Pattern recognition, vol. 37, no. 3, pp. 487–501, 2004.\cr\cr

M. Popescu, J. C. Bezdek, T. C. Havens, J. M. Keller, "A Cluster Validity Framework Based on Induced Partition Dissimilarity," in IEEE Transactions on Cybernetics, vol. 43, no. 1, pp. 308-320, Feb. 2013.\cr\cr

Y. Tang, F. Sun, Z. Sun, “Improved validation index for fuzzy clustering,” in Proceedings of the 2005, American Control Conference, 2005., pp. 1120–1125 vol. 2, 2005.\cr\cr

N. Wiroonsri, O. Preedasawakul, "A correlation-based fuzzy cluster validity index with secondary options detector," arXiv:2308.14785, 2023 \cr\cr

C. H. Wu, C. S. Ouyang, L. W. Chen, L. W. Lu, “A new
fuzzy clustering validity index with a median factor for centroid-based clustering,” IEEE Transactions on Fuzzy Systems, vol. 23, no. 3, pp. 701–718, 2015.\cr\cr

X. Xie, G. Beni, “A validity measure for fuzzy clustering,” IEEE
Transactions on Pattern Analysis and Machine Intelligence, vol. 13, no. 8,
pp. 841–847, 1991.\cr\cr

}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}

\seealso{
\code{\link{WP.IDX}, \link{GC.IDX}, \link{CCV.IDX},  \link{R1_data}}
}

\examples{

library(UniversalCVI)

# Iris data
x = iris[,1:4]

# ---- FCM algorithm ----


# Compute selected a set of indices ("WPC","WP","XB") using default gamma
F.s = FzzyCVIs(scale(x), cmax = 10, cmin = 2, indexlist = c("WPC","WP","XB"),
  corr = 'pearson', method = 'FCM', fzm = 2, iter = 100, nstart = 20, NCstart = TRUE)

# Plot the computed indexes
plot_idx(F.s)

# ---- EM algorithm ----

# Compute all the indices by FzzyCVIs using default gamma
E.all = FzzyCVIs(scale(x), cmax = 10, cmin = 2, indexlist = 'all', corr = 'pearson',
  method = 'EM', iter = 100, nstart = 20, NCstart = TRUE)

# Plot the computed indexes
plot_idx(E.all)

}
