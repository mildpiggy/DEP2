### **Normalization methods: Normalize proteomics data** 

It is conventional to normalize the proteomics assay to minimize the influence of systematic fluctuations. Now, we integrate some alternative normalization methods. It's recommended to check 'Normalization' plot.

#### **Log2 only**

Only log2 transform quantity/intensity, without further normalization methods.


#### **vsn**

Variance stabilizing transformation using the vsn package. VSN eliminate the dependency between variances and mean abundances and scaling data from different samples into the same level through parametric transformations and maximum likelihood estimation.

#### **diff.median**

Median normalization. Centers all samples (columns) so that they all match the grand median by subtracting the respective columns medians differences to the grand median.


#### **quantiles**

Quantile normalization. If inter-class effect proportion and/or batch effects are strong, Quantile may lead to poor performance while also not sufficiently addressing batch effects (doi: 10.1038/s41598-020-72664-6).

