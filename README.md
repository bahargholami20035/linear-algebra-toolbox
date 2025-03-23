# Linear Algebra Toolbox

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Python library implementing common matrix decomposition techniques (LU, Cholesky, QR) for solving systems of linear equations. This project provides educational implementations and does *not* aim to replace highly optimized libraries like NumPy or SciPy for production use.

## Features

*   **LU Decomposition:** `lu_decomp(A)` - Performs LU decomposition *without* pivoting.  Handles singularity detection.
*   **Cholesky Decomposition:** `cholesky_d(A)` - Performs Cholesky decomposition for symmetric, positive-definite matrices. Includes a positive-definiteness check.
*   **QR Decomposition:** `qr_decomp(A)` - Performs QR decomposition using the Gram-Schmidt process. Handles linearly dependent columns detection.
*   **Solvers:**
    *   `solve_lu_system(L, U, b)` - Solves `Ax = b` given the LU decomposition.
    *   `solve_cholesky_system(L, b)` - Solves `Ax = b` given the Cholesky decomposition.
    *   `solve_qr_system(Q, R, b)` - Solves `Ax = b` given the QR decomposition.

## Installation

This is a simple, self-contained library. You can either:

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/YOUR_GITHUB_USERNAME/linear-algebra-toolbox.git
    cd linear-algebra-toolbox
    ```

2.  **Copy the code:** Directly copy the Python code from `linear_algebra_toolbox.py` (or whatever you name the file) into your project.
