# Random-projection emsemble dimension reduction
We introduce a new framework for dimension reduction in the context of high-dimensional regression. The main idea is to apply many low-dimensional random projections to the data, fit a base regression model after each projection and choose good projections according their empirical performance. We then aggregate the chosen projections using singular value decomposition. In general, our goals are:
* Estimate dimension reduction directions, more specifically, central mean subspace (CMS). The performance of our methods can be evaluated by simulation studies.
* Use the dimension reduction directions to further estimate the regression function and make predictions, these can be evaluated through real data analysis.
### Functions to Use ###
* `RPEDR` (Algorithm 1) is our main algorithm, it outputs p dimension reduction directions with decreasing order of importance.
* `RPEDR_dim_est` (Algorithm 2) uses information from Algorithm 1 to estimate the dimension of CMS.
* `RPEDR_double` (Algorithm 3) has the potential to improve the estimation capability compared with Algorithm 1 when the true dimension of CMS is known.
A example to use these functions to reproduce our simulation results can be found in script file `simulation_example.R`
The scripts contained in folder `Real_Data` can be used to reproduce results in our Real world application section.
