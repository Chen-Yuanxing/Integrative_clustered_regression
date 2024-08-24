# Integrative_clustered_regression(ICR)
The files implement the **Integrative clustered regression (ICR)** method proposed in **Heterogeneity-aware Clustered Distributed Learning for Multi-source Data Analysis**, which has been published in *Journal of Machine Learning Research*. ICR aims to achieve the goal of simultaneous estimation, variable selection, and clustering in the framework of distributed learning. 

`main_functions.cpp` includes the main functions about proximal ADMM algorithm, written by Rcpp, of this method. `auxiliary_function.R` implements tuning parameter selection and obtain initial estimates. 
`example5_ICR_slurm` includes the data generating function for Example 5 in this paper and the implementation of the ICR method. The implementation of both the IP and Oracle methods are summarized in `example5_IP_slurm` and `example5_Oracke_slurm`, respectively.
