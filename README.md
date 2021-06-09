## OMP-variations
This repo contains variations of Orthogonal Matching Pursuit (OMP) Algorithm family.

## OMP_MMV (OMP with Multiple Measurement Vector):
A fast implementation of OMP with Multiple Measurement Vector (to recover jointly row-sparse matrix). However, it works for sparse vectors as well

## POMP (Piecewise Orthogonal Matching Pursuit): 
This variation of OMP algorithm finds piecewise sparse vector. Especially handy when you find indices of joint codebooks/sensing matrix. Details here: https://ieeexplore.ieee.org/document/7472550

## DGOMP (Detector based OMP): 
DGOMP is the variation of OMP which doesn't require sparsity (number of non zero elements in the vector) as a prior. It performs statistical hypothesis testing on each iteration on the residual to determine algorithm stopping criterion. Details here: https://link.springer.com/article/10.1186/1687-6180-2014-178
