# UniversalCVI
  **UnversalCVI** package is an effective tool for evaluating clustering results by several cluster validity indices. It has functions to compute several indices as listed below for a user specified range of numbers of clusters and compare them in grid plots. The package is compatible with six clustering methods including K-means, fuzzy C-means, EM clustering, and hierarchical clustering (single, average, and complete linkage). Moreover, the UniversalCVI package has a function that computes the accuracy of clustering results when the true classes are known.

UniversalCVI requires the use of the two packages:  [`mclust`](https://CRAN.R-project.org/package=mclust) and [`e1071`](https://CRAN.R-project.org/package=e1071) which are for performing the fuzzy C-means (FCM) and EM algorithms, respectively.

In addition to the evaluation tools, the UniversalCVI package also includes 17 simulated datasets intially used for testing and comparing cluster validity indices in several perspectives written in Wiroonsri(2024) and Wiroonsri and Preedasawakul(2023). 

The cluster validity indices available in this package are listed as follows:

**Hard clustering:**

Dunn's index, Calinski–Harabasz index, Davies–Bouldin’s index, Point biserial correlation index, Chou-Su-Lai measure, Davies–Bouldin*’s index, Score function, Starczewski index, Pakhira–Bandyopadhyay–Maulik (for crisp clustering) index, and Wiroonsri index.

**Fuzzy clustering:**

Xie–Beni index, KWON index, KWON2 index, TANG index , HF index, Wu–Li index, Pakhira–Bandyopadhyay–Maulik (for fuzzy clustering) index, KPBM index, Correlation Cluster Validity index, Generalized C index, Wiroonsri and Preedasawakul index.


## Installation

If you have not already installed `mclust` and `e1071` in your local system, install these package as following first: 

``` r
install.packages(c('e1071','mclust'))
```
Install `UniversalCVI` package

``` r
install.packages('UniversalCVI')
```
Load R package into R working space

``` r
 suppressPackageStartupMessages({
library(UniversalCVI)
library(e1071)
library(mclust)
})
```

## Example


### Check accuracy of clustering results when the true classes are known

``` r

### Use a dataset in this package
x = R1_data

# Check accuracy of clustering results obtained by kmeans, FCM, and EM clustering
AccClust(x, label.names = "label",algorithm = c("FCM","EM","Kmeans"), fzm = 2,
  scale = TRUE, nstart = 20,iter = 100)
  
x = D1_data

# Check accuracy of a clustering result obtained by the FCM algoritm
AccClust(x, label.names = "label",algorithm = "FCM", fzm = 2,
  scale = TRUE, nstart = 20,iter = 100)
```

### Compute hard cluster validity indices
Using Hvalid to compute all index in function for a clustering result from 2 to 15


``` r
library(UniversalCVI)

# The data is from Wiroonsri (2024).
x = R1_data[,1:2]


# Compute six cluster validity indices of a kmeans clustering result for k from 2 to 15
IDX.list = c("NCI", "DI", "DB", "STR", "CSL", "CH")

Hvalid.result = Hvalid(scale(x), kmax = 15, kmin = 2, indexlist = IDX.list,
  method = "hclust_average", p = 2, q = 2, corr = "pearson", nstart = 100, NCstart = TRUE)

# Plot the computed indices for k from 2 to 15
plot_idx(Hvalid.result)

```

### Soft clustering
Using FzzyCVIs to compute all the fuzzy cluster validity indices for a clustering result for c from 2 to 15


``` r
library(UniversalCVI)

x = R1_data[,1:2]

# Compute six cluster validity indices of a FCM clustering result for c from 2 to 15
IDX.list = c("WP", "PBM", "TANG", "XB", "GC2", "KWON2")
FCVIs = FzzyCVIs(scale(x), cmax = 15, cmin = 2, indexlist = IDX.list, corr = 'pearson',
         method = 'FCM', fzm = 2, iter = 100, nstart = 20, NCstart = TRUE)
# Plot the computed indices for c from 2 to 15
plot_idx(FCVIs)

```

## License

The UniversalCVI package as a whole is distributed under [GPL(>=3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
