import numpy as np

def lu_decomp(a):
    n = len(a)
    l = np.zeros((n, n))
    u = np.zeros((n, n))

    for i in range(n):
        l[i][i] = 1.0
        for j in range(i, n):
            s = sum(l[i][k] * u[k][j] for k in range(i))
            u[i][j] = a[i][j] - s

        for j in range(i + 1, n):
            s = sum(l[j][k] * u[k][i] for k in range(i))
            if np.isclose(u[i][i], 0):  # Check for singularity
                raise ValueError("LU decomposition failed: Zero pivot detected.")
            l[j][i] = (a[j][i] - s) / u[i][i]

    return l, u

def solve_lu_system(l, u, b):
    n = len(l)
    y = np.zeros(n)
    x = np.zeros(n)
    
    for i in range(n):
        y[i] = b[i] - sum(l[i][j] * y[j] for j in range(i))

    for i in range(n - 1, -1, -1):
        if np.isclose(u[i][i], 0):  # Check for division by zero
            raise ValueError("Singular matrix detected in back-substitution.")
        x[i] = (y[i] - sum(u[i][j] * x[j] for j in range(i + 1, n))) / u[i][i]

    return x

def cholesky_d(a):
    n = len(a)
    l = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1):
            tmp_sum = sum(l[i][k] * l[j][k] for k in range(j))
            if i == j:
                val = a[i][i] - tmp_sum
                if val <= 0:  # Check for positive definiteness
                    raise ValueError("Matrix is not positive definite.")
                l[i][j] = np.sqrt(val)
            else:
                l[i][j] = (a[i][j] - tmp_sum) / l[j][j]

    return l

def solve_cholesky_system(l, b):
    n = len(l)
    y = np.zeros(n)
    x = np.zeros(n)

    for i in range(n):
        y[i] = (b[i] - sum(l[i, j] * y[j] for j in range(i))) / l[i, i]

    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(l.T[i, j] * x[j] for j in range(i + 1, n))) / l[i, i]

    return x

def qr_decomp(a_matrix):
    m, n = a_matrix.shape
    q = np.zeros((m, n))
    r = np.zeros((n, n))

    for j in range(n):
        v = a_matrix[:, j]
        for i in range(j):
            r[i, j] = np.dot(q[:, i], a_matrix[:, j])
            v = v - r[i, j] * q[:, i]

        r[j, j] = np.linalg.norm(v)
        if np.isclose(r[j, j], 0):  # Check for linearly dependent columns
            raise ValueError("QR decomposition failed: Linearly dependent column detected.")
        q[:, j] = v / r[j, j]

    return q, r

def solve_qr_system(q_mat, r_mat, b_vec):
    qtb = np.dot(q_mat.T, b_vec)
    n = len(r_mat)
    x_sol = np.zeros(n)

    for i in range(n - 1, -1, -1):
        if np.isclose(r_mat[i, i], 0):  # Check for division by zero
            raise ValueError("Singular matrix detected in QR back-substitution.")
        x_sol[i] = (qtb[i] - sum(r_mat[i, j] * x_sol[j] for j in range(i + 1, n))) / r_mat[i, i]

    return x_sol

# --- Examples ---
A1 = np.array([[2, 1, 1], [4, 3, 3], [8, 7, 9]], dtype=float)
b1 = np.array([1, 3, 2], dtype=float)

L, U = lu_decomp(A1)
x = solve_lu_system(L, U, b1)
print("LU Result:", x)
print("A1 @ x", A1 @ x)

A2 = np.array([[4, 2, -2], [2, 5, 1], [-2, 1, 10]], dtype=float)
b2 = np.array([4, 1, -3], dtype=float)

L2 = cholesky_d(A2)
x2 = solve_cholesky_system(L2, b2)
print("Cholesky Result:", x2)

A3 = np.array([[1, 1, 1], [1, 2, 4], [1, 3, 9]], dtype=float)
b3 = np.array([6, 15, 38], dtype=float)

Q3, R3 = qr_decomp(A3)
x3 = solve_qr_system(Q3, R3, b3)
print("QR Result:", x3)

# Example: A singular matrix. This will now properly detect singularity.
A4 = np.array([[1, 2, 3], [2, 4, 6], [1, 1, 1]], dtype=float)
b4 = np.array([1, 2, 3], dtype=float)

try:
    L4, U4 = lu_decomp(A4)
    x4 = solve_lu_system(L4, U4, b4)
    print("Singular Matrix Result (LU):", x4)
except ValueError as e:
    print(f"LU Error: {e}")

try:
    Q4, R4 = qr_decomp(A4)
    x4 = solve_qr_system(Q4, R4, b4)
    print("Singular Matrix Result (QR):", x4)
except ValueError as e:
    print(f"QR Error: {e}")
