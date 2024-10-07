# Random-projection emsemble dimension reduction
We introduce a new framework for dimension reduction in the context of high-dimensional regression. The main idea is to apply many low-dimensional random projections to the data, fit a base regression model after each projection and choose good projections according their empirical performance. We then aggregate the chosen projections using singular value decomposition. In general, our goals are:
* Estimate dimension reduction directions, more specifically, central mean subspace (CMS). The performance of our methods can be evaluated by simulation studies, 
### Functions to Use ###
* `RPEDR`
