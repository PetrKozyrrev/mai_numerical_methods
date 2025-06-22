import math

def solve_linear_system(A, b):
    a = [A[0] + [b[0]],
         A[1] + [b[1]]]

    if a[0][0] == 0:
        if a[1][0] != 0:
            a[0], a[1] = a[1], a[0]
        else:
            raise ValueError("Система не имеет уникального решения.")

    factor = a[1][0] / a[0][0]
    for j in range(3):
        a[1][j] = a[1][j] - factor * a[0][j]

    if a[1][1] == 0:
        raise ValueError("Система не имеет уникального решения.")

    x_2 = a[1][2] / a[1][1]
    x_1 = (a[0][2] - a[0][1] * x_2) / a[0][0]

    return x_1, x_2


def jacobian(x1, x2, df1_dx1, df1_dx2, df2_dx1, df2_dx2):
    return [[df1_dx1(x1, x2), df1_dx2(x1, x2)],
            [df2_dx1(x1, x2), df2_dx2(x1, x2)]]


def max_jacobian_norm(x1, x2, dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2):
    jacobian_matrix = jacobian(x1, x2, dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2)

    jacobian_norm_1 = abs(jacobian_matrix[0][0]) + abs(jacobian_matrix[0][1])
    jacobian_norm_2 = abs(jacobian_matrix[1][0]) + abs(jacobian_matrix[1][1])

    jacobian_norm_3 = abs(jacobian_matrix[0][0]) + abs(jacobian_matrix[1][0])
    jacobian_norm_4 = abs(jacobian_matrix[0][1]) + abs(jacobian_matrix[1][1])

    return min(max(jacobian_norm_1, jacobian_norm_2),max(jacobian_norm_3, jacobian_norm_4))


def newton_method(x1_init, x2_init, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, tol=1e-4, max_iter=100):
    x1, x2 = x1_init, x2_init
    for iteration in range(max_iter):
        f_v = [f1(x1, x2), f2(x1, x2)]

        J = jacobian(x1, x2, df1_dx1, df1_dx2, df2_dx1, df2_dx2)

        try:
            d1, d2 = solve_linear_system(J, [-f_v[0], -f_v[1]])
        except ValueError as e:
            print(f"Ошибка на итерации {iteration + 1}: {e}")
            return None, None

        x1_new = x1 + d1
        x2_new = x2 + d2

        max_change = max(abs(x1_new - x1), abs(x2_new - x2))

        if max_change < tol:
            print(f"Сходимость достигнута за {iteration + 1} итераций.")
            return x1_new, x2_new

        x1, x2 = x1_new, x2_new

    print("Достигнуто максимальное число итераций.")
    return x1, x2


def simple_iteration(x1_init, x2_init, phi1, phi2, dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2, epsilon,
                     max_iter):
    x1 = x1_init
    x2 = x2_init

    try:
        q = max_jacobian_norm(x1, x2, dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2)
    except Exception as e:
        print(f"Ошибка при вычислении q: {e}")
        return (None, None), 0

    if q >= 1:
        print("Метод не сходится, т.к. q >= 1")
        return (None, None), 0

    for k in range(max_iter):
        try:
            x1_new = phi1(x1, x2)
            x2_new = phi2(x1, x2)
        except Exception as e:
            print(f"Ошибка на итерации {k + 1}: {e}")
            return (None, None), k + 1

        error = max(abs(x1_new - x1), abs(x2_new - x2))

        if (q / (1 - q)) * error < epsilon:
            return (x1_new, x2_new), k + 1

        x1, x2 = x1_new, x2_new
    print(f"Метод не сошелся за {max_iter} итераций")
    return (None, None), max_iter


def f1(x1, x2):
    return 1 / 4 * x1 ** 2 + x2 ** 2 - 1


def f2(x1, x2):
    return 2 * x2 - math.exp(x1) - x1


def df1_dx1(x1, x2):
    return 1 / 2 * x1


def df1_dx2(x1, x2):
    return 2 * x2


def df2_dx1(x1, x2):
    return -1 - math.exp(x1)


def df2_dx2(x1, x2):
    return 2


def phi1(x1, x2):
    return math.log(2*x2 - x1)


def phi2(x1, x2):
    return math.sqrt(1 - (x1**2/4))


def dphi1_dx1(x1, x2):
    return (-1)/(2*x2-x1)


def dphi1_dx2(x1, x2):
    return 2/(2*x2 - x1)


def dphi2_dx1(x1, x2):
    return (-x1)/(2 * math.sqrt(4 - x1**2))


def dphi2_dx2(x1, x2):
    return 0


x1_init = 1
x2_init = 2
epsilon = 1e-8
max_iter = 1000

print("======= Метод Ньютона =======")
solution = newton_method(x1_init, x2_init, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2)
print("Решение Ньютоном:", solution)

print("\n======= Метод простых итераций =======")
result, iterations = simple_iteration(x1_init, x2_init, phi1, phi2, dphi1_dx1, dphi1_dx2, dphi2_dx1, dphi2_dx2,
                                      epsilon, max_iter)
if result[0] is not None:
    print(f"Решение: x1 = {result[0]}, x2 = {result[1]} достигнуто за {iterations} итераций")
else:
    print("Решение не найдено")