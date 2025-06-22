import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

def plot_2(Xi, Yi, X, Y, title=''):
    plt.figure(figsize=(8, 6))
    plt.scatter(Xi, Yi, c='orange', label='Численное решение', zorder=0.5)
    plt.plot(X, Y, c='blue', label='Точное решение', linewidth=2)
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_shooting_iterations(all_paths, X_exact, Y_exact):
    plt.figure(figsize=(8, 6))
    plt.plot(X_exact, Y_exact, c='blue', label='Точное решение', linewidth=1.5)
    for i, (Xk, Yk) in enumerate(all_paths):
        plt.plot(Xk, Yk, label=f'Итерация {i}', linewidth=2, linestyle='--')

    plt.title("Метод стрельбы – итерации")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def runge_kutta(x0, y0, z0, x_end, h, dz, dy):
    xk, yk, zk = x0, y0, z0
    Xk, Yk, Zk = [xk], [yk], [zk]
    while abs(xk - x_end) > abs(h) / 10:
        L1 = h * dz(xk, yk, zk)
        K1 = h * dy(xk, yk, zk)

        L2 = h * dz(xk + h / 2, yk + K1 / 2, zk + L1 / 2)
        K2 = h * dy(xk + h / 2, yk + K1 / 2, zk + L1 / 2)

        L3 = h * dz(xk + h / 2, yk + K2 / 2, zk + L2 / 2)
        K3 = h * dy(xk + h / 2, yk + K2 / 2, zk + L2 / 2)

        L4 = h * dz(xk + h, yk + K3, zk + L3)
        K4 = h * dy(xk + h, yk + K3, zk + L3)

        yk += (K1 + 2 * K2 + 2 * K3 + K4) / 6
        zk += (L1 + 2 * L2 + 2 * L3 + L4) / 6
        xk += h

        Xk.append(xk)
        Yk.append(yk)
        Zk.append(zk)

    return Xk, Yk, Zk


def shooting_method(mu1, mu2, h):
    a, b = 0, 1
    z_a = 0.75  # y'(0)
    z_b = (math.exp(2) * (math.e + 2)) / ((math.e + 1) ** 2)  # y'(1)
    eps = 1e-6

    all_paths = []

    def dy(x, y, z): return z

    def dz(x, y, z): return (2 * z + math.exp(x) * y) / (math.exp(x) + 1)

    def run(mu):
        return runge_kutta(a, mu, z_a, b, h, dz, dy)

    X1, Y1, Z1 = run(mu1)
    X2, Y2, Z2 = run(mu2)

    all_paths.append((X1, Y1))
    all_paths.append((X2, Y2))

    phi1 = Z1[-1] - z_b
    phi2 = Z2[-1] - z_b

    while abs(phi2) > eps:
        dphi = (phi2 - phi1) / (mu2 - mu1)
        mu1, phi1 = mu2, phi2
        mu2 = mu2 - phi2 / dphi
        Xk, Yk, Zk = run(mu2)
        phi2 = Zk[-1] - z_b
        all_paths.append((Xk, Yk))

    return mu2, Xk, Yk, all_paths


def runge_romberg_method(get_func, h, r, p):
    y1 = get_func(h)
    y2 = get_func(h * r)
    y2_reduced = y2[::r]
    min_len = min(len(y1), len(y2_reduced))
    max_error = 0
    for i in range(min_len):
        err = abs((y1[i] - y2_reduced[i]) / (r ** p - 1))
        if err > max_error:
            max_error = err
    return max_error



def lu_decomposition(A):
    n = len(A)
    LU = copy(A)
    for k in range(n):
        for i in range(k + 1, n):
            LU[i][k] /= LU[k][k]
            for j in range(k + 1, n):
                LU[i][j] -= LU[i][k] * LU[k][j]
    return LU

def forward_substitution(L, b):
    n = len(L)
    y = [0] * n

    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
    return y


def backward_substitution(U, y):
    n = len(U)
    x = [0] * n

    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]

    return x

def copy(A):
    return [[A[i][j] for j in range(len(A[i]))] for i in range(len(A))]

def solve_system(A, b):
    n = len(A)
    LU = lu_decomposition(A)

    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            L[i][j] = LU[i][j]
        for j in range(i, n):
            U[i][j] = LU[i][j]

    y = forward_substitution(L, b)

    x = backward_substitution(U, y)
    return x

def finite_difference(N):
    a, b = 0.0, 1.0
    h = (b - a) / N
    x = [a + i*h for i in range(N+1)]

    def A(x): return math.exp(x) + 1
    def B(x): return -2
    def C(x): return -math.exp(x)
    def F(x): return 0.0

    alpha = 3/4
    beta = (math.exp(2) * (math.exp(1) + 2)) / (math.exp(1) + 1)**2

    size = N + 1
    A_matrix = [[0.0 for _ in range(size)] for _ in range(size)]
    d = [0.0 for _ in range(size)]

    A_matrix[0][0] = -3/(2*h)
    A_matrix[0][1] = 4/(2*h)
    A_matrix[0][2] = -1/(2*h)
    d[0] = alpha

    for i in range(1, N):
        xi = x[i]
        Ai = A(xi)
        Bi = B(xi)
        Ci = C(xi)
        Fi = F(xi)

        A_matrix[i][i-1] = Ai / h**2 - Bi / (2*h)
        A_matrix[i][i] = -2 * Ai / h**2 + Ci
        A_matrix[i][i+1] = Ai / h**2 + Bi / (2*h)
        d[i] = Fi

    A_matrix[N][N] = 3 / (2*h)
    A_matrix[N][N-1] = -4 / (2*h)
    A_matrix[N][N-2] = 1 / (2*h)
    d[N] = beta


    y = solve_system(A_matrix,d)

    return x, y

def f_exact(x):
    return math.exp(x) - 1 + 1 / (math.exp(x) + 1)

def mae(y1, y2):
    assert len(y1) == len(y2)
    return sum(abs(y1[i] - y2[i]) for i in range(len(y1))) / len(y1)


mu, X_sh, Y_sh, all_paths = shooting_method(mu1=0.6, mu2=1.5, h=0.01)

X = [i * 0.01 for i in range(101)]
Y = [f_exact(x) for x in X]


plot_shooting_iterations(all_paths, X, Y)

print("Метод стрельбы:")
print(f"Ошибка Рунге–Ромберга: {runge_romberg_method(lambda h: shooting_method(0.6, 1.5, h)[2], 0.01, 2, 2)}")
print(f"Ошибка в сравнении с точным решением: {mae(Y_sh, Y)}")
plot_2(X_sh, Y_sh, X, Y, "Метод стрельбы")


get_fd = lambda n: finite_difference(n)[1]
X_fd, Y_fd = finite_difference(100)

print("\nМетод конечных разностей:")
print(f"Ошибка Рунге–Ромберга: {runge_romberg_method(get_fd, 100, 2, 1)}")
print(f"Ошибка в сравнении с точным решением: {mae(Y_fd, Y)}")
plot_2(X_fd, Y_fd, X, Y, "Метод конечных разностей")

