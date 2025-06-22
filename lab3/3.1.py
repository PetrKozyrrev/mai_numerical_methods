import matplotlib.pyplot as plt
import math

def f(x):
    return math.cos(x) + x


X_a = [0, math.pi / 6, 2 * math.pi / 6, 3 * math.pi / 6]
X_b = [0, math.pi / 6, math.pi / 4, math.pi / 2]


Y_a = [f(x) for x in X_a]
Y_b = [f(x) for x in X_b]


def lagrange_interpolation(X, Y, x):
    n = len(X)
    L = 0
    for i in range(n):
        li = 1
        for j in range(n):
            if j != i:
                li *= (x - X[j]) / (X[i] - X[j])
        L += Y[i] * li
    return L



def divided_differences(X, Y):
    n = len(X)
    diff_table = [Y[:]]
    for k in range(1, n):
        diff_table.append([])
        for i in range(n - k):
            diff = (diff_table[k - 1][i + 1] - diff_table[k - 1][i]) / (X[i + k] - X[i])
            diff_table[k].append(diff)
    return diff_table


def newton_interpolation(X, Y, x):
    n = len(X)
    diff_table = divided_differences(X, Y)
    result = diff_table[0][0]
    product_term = 1
    for i in range(1, n):
        product_term *= (x - X[i - 1])
        result += diff_table[i][0] * product_term
    return result


x_star = 1.0


lagrange_result_a = lagrange_interpolation(X_a, Y_a, x_star)
lagrange_result_b = lagrange_interpolation(X_b, Y_b, x_star)

newton_result_a = newton_interpolation(X_a, Y_a, x_star)
newton_result_b = newton_interpolation(X_b, Y_b, x_star)


print(f"f(x*) = {f(x_star)}")
print("---------------------------------------------------------------------------------")
print("a)")
print(f"Результат интерполяции методом Лагранжа в точке x* = {x_star}: {lagrange_result_a}")
print(f"Погрешность интерполяции методом Лагранжа: {abs(lagrange_result_a - f(x_star))}")

print("b)")
print(f"Результат интерполяции методом Лагранжа в точке x* = {x_star}: {lagrange_result_b}")
print(f"Погрешность интерполяции методом Лагранжа: {abs(lagrange_result_b - f(x_star))}")

print("---------------------------------------------------------------------------------")
print("a)")
print(f"Результат интерполяции методом Ньютона в точке x* = {x_star}: {newton_result_a}")
print(f"Погрешность интерполяции методом Ньютона: {abs(newton_result_a - f(x_star))}")

print("b)")
print(f"Результат интерполяции методом Ньютона в точке x* = {x_star}: {newton_result_b}")
print(f"Погрешность интерполяции методом Ньютона: {abs(newton_result_b - f(x_star))}")

fig = plt.figure(figsize=(7, 7))
X_l1 = []
Y_l1 = []

X_n1 = []
Y_n1 = []
i = 0
while i <= 3:
    X_l1.append(i)
    Y_l1.append(lagrange_interpolation(X_a, Y_a, i))

    X_n1.append(i)
    Y_n1.append(newton_interpolation(X_a, Y_a, i))
    i += 0.01

plt.scatter(X_a, Y_a, c='r')
plt.scatter([x_star], [f(x_star)], c='y')
plt.plot(X_l1, Y_l1, c='g')
plt.plot(X_n1, Y_n1, c='b')
plt.title("а) Многочлены Лагранжа и Ньютона")
plt.show()

fig = plt.figure(figsize=(7, 7))
X_l2 = []
Y_l2 = []

X_n2 = []
Y_n2 = []
i = 0
while i <= 3:
    X_l2.append(i)
    Y_l2.append(lagrange_interpolation(X_b, Y_b, i))

    X_n2.append(i)
    Y_n2.append(newton_interpolation(X_b, Y_b, i))
    i += 0.01

plt.scatter(X_b, Y_b, c='r')
plt.scatter([x_star], [f(x_star)], c='y')
plt.plot(X_l2, Y_l2, c='r')
plt.plot(X_n2, Y_n2, c='g')
plt.title("б) Многочлены Лагранжа и Ньютона")
plt.show()