import matplotlib.pyplot as plt

n = 4
X_i = [0.0, 1.0, 2.0, 3.0, 4.0]
f_i = [1.0, 1.5403, 1.5839, 2.01, 3.3464]
task_function_value_point_2 = 1.5

h = [0]
for i in range(1, len(X_i)):
    h.append(X_i[i] - X_i[i - 1])

A = [[2 * (h[1] + h[2]), h[2], 0]]
_d = [3 * ((f_i[2] - f_i[1]) / h[2] - (f_i[1] - f_i[0]) / h[1])]

for i in range(3, n):
    next_row_A = [0] * len(A[0])
    next_row_A[i - 3] = h[i - 1]
    next_row_A[i - 2] = 2 * (h[i - 1] + h[i])
    next_row_A[i - 1] = h[i]
    A.append(next_row_A)

    _d.append(3 * ((f_i[i] - f_i[i - 1]) / h[i] - (f_i[i - 1] - f_i[i - 2]) / h[i - 1]))

next_row_A = [0] * len(A[0])
next_row_A[-1] = 2 * (h[n - 1] + h[n])
next_row_A[-2] = h[n - 1]
A.append(next_row_A)

_d.append(3 * ((f_i[n] - f_i[n - 1]) / h[n] - (f_i[n - 1] - f_i[n - 2]) / h[n - 1]))


def TDMAsolver(A, _d):
    n = len(_d)
    a = []
    b = []
    c = []
    d = _d
    ACopy = A

    a.append(0)
    b.append(ACopy[0][0])
    c.append(ACopy[0][1])
    for i in range(1, n - 1):
        a.append(ACopy[i][i - 1])
        b.append(ACopy[i][i])
        c.append(ACopy[i][i + 1])

    a.append(ACopy[n - 1][n - 2])
    b.append(ACopy[n - 1][n - 1])
    c.append(0)

    # creating P and Q
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    for i in range(1, n - 1):
        p.append(- c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))

    ans = [[0 for _ in range(n)] for _ in range(n)]
    ans[n - 1] = (d[n - 1] - a[n - 1] * q[n - 2]) / (b[n - 1] + a[n - 1] * p[n - 2])

    for i in range(n - 2, -1, -1):
        ans[i] = p[i] * ans[i + 1] + q[i]

    return ans


c = [0] + TDMAsolver(A, _d)
a = f_i[:-1]

b = [(f_i[i] - f_i[i - 1]) / h[i] - (1 / 3) * h[i] * (c[i] + 2 * c[i - 1]) for i in range(1, n)]
b.append((f_i[n] - f_i[n - 1]) / h[n] - (2 / 3) * h[n] * c[n - 1])

d = [(c[i] - c[i - 1]) / (3 * h[i]) for i in range(1, n)]
d.append(-c[n - 1] / (3 * h[n]))


# print("a:", a)
# print("b:", b)
# print("c:", c)
# print("d:", d)


def S(x, a, b, c, d, X):
    for i in range(1, len(X)):
        if x >= X[i - 1] and x <= X[i]:
            return a[i - 1] + b[i - 1] * (x - X[i - 1]) + c[i - 1] * (x - X[i - 1]) ** 2 + d[i - 1] * (
                    x - X[i - 1]) ** 3


print(f"Значение в точке X*: {S(task_function_value_point_2, a, b, c, d, X_i)}")

fig = plt.figure(figsize=(7, 7))
Xi = []
Yi = []
i = 0
while i <= 4:
    Xi.append(i)
    Yi.append(S(i, a, b, c, d, X_i))
    i += 0.01

plt.scatter(X_i, f_i, c='r')
plt.plot(Xi, Yi, c='g')
plt.show()
