# Linear Algebra Toolbox

A Python library implementing common matrix decomposition techniques (LU, Cholesky, QR) for solving systems of linear equations. This project provides educational implementations and does *not* aim to replace highly optimized libraries like NumPy or SciPy for production use.

## Features

*   **LU Decomposition:** `lu_decomp(A)` - Performs LU decomposition *without* pivoting.  Handles singularity detection.
*   **Cholesky Decomposition:** `cholesky_d(A)` - Performs Cholesky decomposition for symmetric, positive-definite matrices. Includes a positive-definiteness check.
*   **QR Decomposition:** `qr_decomp(A)` - Performs QR decomposition using the Gram-Schmidt process. Handles linearly dependent columns detection.
*   **Solvers:**
    *   `solve_lu_system(L, U, b)` - Solves `Ax = b` given the LU decomposition.
    *   `solve_cholesky_system(L, b)` - Solves `Ax = b` given the Cholesky decomposition.
    *   `solve_qr_system(Q, R, b)` - Solves `Ax = b` given the QR decomposition.

