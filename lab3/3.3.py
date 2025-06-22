import matplotlib.pyplot as plt

X_i = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
Y_i = [-0.4597, 1.0, 1.5403, 1.5839, 2.010, 3.3464]


def lu_decomposition(LU):
    n = len(LU)
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


def mnk_func(x, a):
    result = 0
    for i, ai in enumerate(a):
        result += ai * (x ** i)
    return result


def sum_sq_er(Xi, Yi, a, f):
    sum = 0
    for i, y in enumerate(Yi):
        sum += (y - f(Xi[i], a)) ** 2
    return sum


def get_a(k, Xi, Yi):
    F = []
    for j, x in enumerate(Xi):
        ph = [x ** i for i in range(k + 1)]
        F.append(ph)

    G = [[0] * (k + 1) for _ in range(k + 1)]
    for i in range(k + 1):
        for j in range(k + 1):
            G[i][j] = sum(F[m][i] * F[m][j] for m in range(len(Xi)))

    z = [sum(F[m][i] * Yi[m] for m in range(len(Xi))) for i in range(k + 1)]

    LU = lu_decomposition(G)

    n = len(LU)
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            L[i][j] = LU[i][j]
        for j in range(i, n):
            U[i][j] = LU[i][j]

    y = forward_substitution(L, z)

    a = backward_substitution(U, y)

    return a


# Получаем коэффициенты для полиномов 1-й и 2-й степени
a1 = get_a(1, X_i, Y_i)
a2 = get_a(2, X_i, Y_i)

# Генерация точек для графика
X = [x * 0.01 + X_i[0] for x in range(int((X_i[-1] - X_i[0]) / 0.01) + 1)]
Y1 = [mnk_func(x, a1) for x in X]
Y2 = [mnk_func(x, a2) for x in X]


def Plot3(Xi, Yi, X, Y1, Y2, l1='', l2=''):
    fig = plt.figure(figsize=(7, 7))
    plt.scatter(Xi, Yi, c='r')
    plt.plot(X, Y1, c='g', label=l1)
    plt.plot(X, Y2, c='b', label=l2)
    plt.legend()
    plt.show()


print("Сумма квадратов ошибок для приближающего многочлена 1-й степени: ", sum_sq_er(X_i, Y_i, a1, mnk_func))
print("Сумма квадратов ошибок для приближающего многочлена 2-й степени: ", sum_sq_er(X_i, Y_i, a2, mnk_func))

Plot3(X_i, Y_i, X, Y1, Y2, "1-й степени", "2-й степени")