# **impute-methods: Quantitative proteomics data imputation**

The imputation methods clustered by `DEP2` package

#### **man**

Imputes missing values in a proteomics dataset by random draws from a manually defined distribution (scale = 0.3, shift = 1.8)

#### **bcpa**

Bayesian missing value imputation are available, as implemented in the and `pcaMethods::pca` functions. See [pca](https://www.rdocumentation.org/packages/tofsims/versions/1.0.2/topics/PCA-class) for deta

#### **knn**

Nearest neighbour averaging, as implemented in the `impute::impute.knn` function. See [impute.knn](https://www.rdocumentation.org/packages/impute/versions/1.46.0/topics/impute.knn) for details and additional parameters.

#### **QRILC**

A missing data imputation method that performs the imputation of left-censored missing data using random draws from a truncated distribution with parameters estimated using quantile regression. Implemented in the `imputeLCMD::impute.QRILC` function. See [impute.QRILC](https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.QRILC) for details and additional parameters.

#### **MLE**

Maximum likelihood-based imputation method using the EM algorithm. Implemented in the `norm::imp.norm` function. See [imp.norm](https://www.rdocumentation.org/packages/norm/versions/1.0-9.5/topics/imp.norm) for details and additional parameters. 

#### **MinDet**

Performs the imputation of left-censored missing data using a deterministic minimal value approach. Considering a expression data with *n* samples and *p* features, for each sample, the missing entries are replaced with a minimal value observed in that sample. The minimal value observed is estimated as being the q-th quantile (default `q = 0.01`) of the observed values in that sample. Implemented in the `imputeLCMD::impute.MinDet` function. See [impute.MinDet](https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.MinDet) for details and additional parameters.

#### **MinProb**

Performs the imputation of left-censored missing data by random draws from a Gaussian distribution centred to a minimal value. Considering an expression data matrix with *n* samples and *p* features, for each sample, the mean value of the Gaussian distribution is set to a minimal observed value in that sample. The minimal value observed is estimated as being the q-th quantile (default `q = 0.01`) of the observed values in that sample. The standard deviation is estimated as the median of the feature standard deviations. Note that when estimating the standard deviation of the Gaussian distribution, only the peptides/proteins which present more than 50% recorded values are considered. Implemented in the `imputeLCMD::impute.MinProb` function. See [impute.MinProb](https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.MinProb) for details and additional parameters.

#### **min**

Replaces the missing values by the smallest non-missing value in the data.

#### **zero**

Replaces the missing values by 0.

#### **mixed on proteins**

It uses a MAR(missing values at random) and MNAR(missing values not at random) imputation method on different subsets of proteins. Here, we consider a protein to have missing values not at random (MNAR) if it has missing values in all replicates of at least one condition, and use MinDet imputation method for MANR, knn for MAR.

#### **mixed on samples**

It is performed on a subset of samples. we imputed the **Control** sample using the “MinProb” method and the samples using the “knn” method.

#### **RF**

It impute missing values by **missForest** approach (provided by `missForest` package), an iterative imputation method (missForest) based on a random forest using the built-in out-of-bag error estimates. See [MissForest](https://academic.oup.com/bioinformatics/article/28/1/112/219101) for more information.

