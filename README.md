# Compressive recovery of a low-rank matrix
**[[Full Text][4]]**

## Introduction
In many practical problems of interest, one would like to estimate a matrix from a sampling of its entries. The matrix we wish to estimate is known to be structured in the sense that it is
low-rank or approximately low-rank. Given below is a practical scenario where one would
like to be able to recover a low-rank matrix from a sampling of its entries.

_The Netflix problem_ - Users are given the opportunity to rate movies, but users typically rate only very few movies so that there are very few scattered
observed entries of this data matrix. Yet one would like to complete this matrix so that Netflix
might recommend titles that any particular user is likely to be willing to order. In this case, the data matrix of all
user-ratings may be approximately low-rank because it is commonly believed that only a few factors contribute
to an individual’s tastes or preferences.

In this project, I present an Alternating Direction Method of Multipliers
(ADMM) algorithm that has been widely used for solving several convex and non-convex optimization problems
by breaking them into smaller sub-problems. Further, I even went beyond and improved the results presented in the paper by imposing a sparsity constraint
on the matrix in the DCT basis. This is a common condition for natural images, and I demonstrate a significant
improvement in reconstruction results.

## Results
Given below are the reconstructions in the case of imposing the sparsity
constraint and not imposing the sparsity constraint. As we can see, imposing the sparsity constraint significantly improves the reconstruction results.

<div class="row">
  <div class="column">
    <img src="/assets/low_rank/matrix.png">
  </div>
</div>

[4]: https://arunabh98.github.io/reports/matrix_completion.pdf
