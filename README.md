# GRridge
R package for better prediction by use of co-data: Adaptive group-regularized ridge regression

This is the DEVELOPMENT VERSION of GRridge. Additional functionality includes

- Default Method in grridge(): ExactStable, including regularization of coefficient matrix
- Allows to enter one partition (grouping) as a simple list
- Added checks on whether indices in the partition are valid
- Added post-weighting variable selection using elastic net
- Allow re-calibration during testing
- Allow survival response
