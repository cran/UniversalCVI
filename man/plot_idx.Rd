\name{plot_idx}
\alias{plot_idx}
\title{
Plots for visualizing CVIs
}
\description{
Plot and compare upto 8 indices computed by the algorithms in this package.
}
\usage{
plot_idx(idxresult,selected.idx = NULL)
}
\arguments{
  \item{idxresult}{a result from one of the algorithms \code{FzzyCVIs, WP.IDX, GC.IDX, CCV.IDX, XB.IDX, WL.IDX, TANG.IDX, PBM.IDX, KWON.IDX, KWON2.IDX, KPBM.IDX, HF.IDX, Hvalid, Wvalid, SF.IDX, PB.IDX, DI.IDX, DB.IDX, CSL.IDX, CH.IDX or STRPBM.IDX}.
}
  \item{selected.idx}{a numeric vector indicates a part of the indexes from the \code{idxresult} in respective order selected by a user. For instance, \code{selected.idx = 3} or \code{selected.idx = c(1,3,5)} may be selected. If not specified, the full idxresult will be considered.
  }
}
\value{
Plots of upto 8 cluster validity indices computed from \code{FzzyCVIs, WP.IDX, GC.IDX, CCV.IDX, XB.IDX, WL.IDX, TANG.IDX, PBM.IDX, KWON.IDX, KWON2.IDX, KPBM.IDX, HF.IDX, Hvalid, Wvalid, SF.IDX, PB.IDX, DI.IDX, DB.IDX, CSL.IDX, CH.IDX or STRPBM.IDX}. When using the isolated index algorithm, all the plots computed by that algorithm will be shown. When using \code{FzzyCVIs or Hvalid} with more than 8 selected indices, the first 8 indices will be plotted.
}
\references{
N. Wiroonsri, O. Preedasawakul, "A correlation-based fuzzy cluster validity index with secondary options detector," arXiv:2308.14785, 2023
}
\author{
Nathakhun Wiroonsri and Onthada Preedasawakul
}
\seealso{
\code{\link{FzzyCVIs}, \link{WP.IDX}, \link{XB.IDX}, \link{Hvalid}}
}
\examples{

library(UniversalCVI)

# Iris data
x = iris[,1:4]

# ----Compute all the indices by FzzyCVIs ----
FCVIs = FzzyCVIs(scale(x), cmax = 10, cmin = 2, indexlist = 'all', corr = 'pearson',
                 method = 'FCM', fzm = 2, iter = 100, nstart = 20, NCstart = TRUE)

# plots of the eight indices by default
plot_idx(idxresult = FCVIs)

# plots of a specific selected.idx
plot_idx(idxresult = FCVIs, selected.idx = c(2,5,7))

# ----Compute all the indices by Wvalid ----
FCM.NC = Wvalid(scale(x), kmax = 10, kmin=2, method = 'kmeans',
  corr='pearson', nstart=100, NCstart = TRUE)

# plots of the four indices by default
plot_idx(idxresult = FCM.NC)

# ----Compute all the indices by XB.IDX ----

FCM.XB = XB.IDX(scale(x), cmax = 10, cmin = 2, method = "FCM",
  fzm = 2, nstart = 20, iter = 100)
plot_idx(idxresult = FCM.XB)
}
