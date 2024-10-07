# Random-projection emsemble dimension reduction
We introduce a new framework for dimension reduction in the context of high-dimensional regression. The main idea is to apply many low-dimensional random projections to the data, fit a base regression model after each projection and choose good projections according their empirical performance. We then aggregate the chosen projections using singular value decomposition. 

## Objectives
- **Estimate Dimension Reduction Directions**: Specifically, we aim to estimate the central mean subspace (CMS). The performance of our methods can be evaluated through simulation studies.
- **Prediction via Dimension Reduction**: Using the estimated dimension reduction directions, we further estimate the regression function and make predictions, which can be evaluated through real data analysis.

## Key Functions
- `RPEDR` (Algorithm 1): The main algorithm, which outputs $p$ dimension reduction directions, ordered by their importance.
- `RPEDR_dim_est` (Algorithm 2): Utilizes the output of Algorithm 1 to estimate the dimension of the CMS.
- `RPEDR_double` (Algorithm 3): A variant of Algorithm 1, designed to improve estimation performance when the true dimension of the CMS is known.

## Usage Example
To reproduce our simulation results, refer to the script file `simulation_example.R`.

Scripts contained in the `Real_Data` folder can be used to replicate the results in the real-world application section of our work.
