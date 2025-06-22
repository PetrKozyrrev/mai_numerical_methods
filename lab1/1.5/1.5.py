import math

def mat_copy(A):
    return [row[:] for row in A]

def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    p = len(B)
    result = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(p):
                result[i][j] += A[i][k] * B[k][j]
    return result

def mat_identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def vector_norm(v):
    return math.sqrt(sum(x ** 2 for x in v))

def vector_add(u, v):
    return [a + b for a, b in zip(u, v)]

def vector_sub(u, v):
    return [a - b for a, b in zip(u, v)]

def vector_mul_scalar(v, s):
    return [x * s for x in v]

def qr_decomposition(A):
    n = len(A)
    Q = mat_identity(n)
    R = mat_copy(A)

    for k in range(n - 1):

        x = [R[i][k] for i in range(k, n)]

        e1 = [1.0] + [0.0] * (n - k - 1)

        sign = 1.0 if x[0] >= 0 else -1.0

        norm_x = vector_norm(x)
        v = vector_add(x, vector_mul_scalar(e1, sign * norm_x))
        norm_v = vector_norm(v)
        v = vector_mul_scalar(v, 1.0 / norm_v)

        H = mat_identity(n)
        for i in range(k, n):
            for j in range(k, n):
                H[i][j] -= 2.0 * v[i - k] * v[j - k]

        Q = mat_mul(Q, H)
        R = mat_mul(H, R)

    return Q, R

def qr_algorithm(A, max_iterations=100,eps = -1e9):
    Ak = mat_copy(A)
    n = len(Ak)
    for k in range(max_iterations):
        Q, R = qr_decomposition(Ak)
        Ak = mat_mul(R, Q)

        for i in R:
            print(*i)
        print("")

        for i in range(n):
            sm = 0
            for j in range(i,n):
                sm += A[i][j]**2

            if(math.sqrt(sm) < eps):
                break

    for i in Ak:
        print(*i)

    return [Ak[i][i] for i in range(n)]

A = [
    [-1, 2, 9],
    [9, 3, 4],
    [8, -4, -6]
]

eigenvalues = qr_algorithm(A)
print("Приближённые собственные значения:", eigenvalues)
