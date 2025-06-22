# Чтение матрицы из файла
def read_input_file():
    with open('1.3.txt', 'r') as f:
        n = int(f.readline())
        A = []
        for i in range(n):
            A.append(list(map(int, f.readline().split())))
        b = list(map(int, f.readline().split()))
        epsilon = float(f.readline().strip())
    return n, A, b, epsilon


#  диагональное преобладание матрицы A (по строкам)
def diag_condition_1(A):
    n = len(A)
    for i in range(n):
        s = 0
        for j in range(len(A[i])):
            if i != j: s += abs(A[i][j])
        if abs(A[i][i]) <= s:
            return False
    return True

#  диагональное преобладание матрицы A (по столбцам)
def diag_condition_2(A):
    n = len(A)
    for i in range(n):
        s = 0
        for j in range(len(A[i])):
            if i != j: s += abs(A[j][i])
        if abs(A[i][i]) <= s:
            return False
    return True

def simple_iterations_solve(a, b, eps):
    if (diag_condition_1(a) or diag_condition_2(a)):
        print("Достаточное условие выполнено")
    else:
        print("Достаточное условие НЕ выполнено")

    eps_k = 1
    x = [b[i] / a[i][i] for i in range(len(a))]
    iter_count = 0
    print(f"Eps = {eps}")
    while (eps < eps_k and iter_count < 100):
        x_k = [0.0] * n
        for i in range(n):
            s = 0.0
            for j in range(n):
                if j != i:
                    s += (-a[i][j] / a[i][i]) * x[j]
            x_k[i] = (b[i] / a[i][i]) + s

        eps_k = max(abs(x_k[i] - x[i]) for i in range(n))

        x = [i for i in x_k]

        iter_count += 1
        print(f"iter = {iter_count}, eps_k = {eps_k}, x = {x}")

    return x, iter_count


def zeidel_solve(a, b, eps):
    n = len(a)
    if (diag_condition_1(a) or diag_condition_2(a)):
        print("Достаточное условие выполнено")
    else:
        print("Достаточное условие НЕ выполнено")

    eps_k = 1
    x = [b[i] / a[i][i] for i in range(len(a))]
    iter_count = 0
    print(f"Eps = {eps}")
    while (eps < eps_k and iter_count < 100):

        x_k = [i for i in x]

        for i in range(n):
            s = 0.0
            # Обновлённые значения для j < i
            for j in range(i):
                s += a[i][j] * x[j]
            # Старые значения для j > i
            for j in range(i + 1, n):
                s += a[i][j] * x_k[j]
            x[i] = (b[i] - s) / a[i][i]

        eps_k = max(abs(x[i] - x_k[i]) for i in range(n))

        iter_count += 1
        print(f"iter = {iter_count}, eps_k = {eps_k}, x = {x}")

    return x, iter_count



n, A, b, epsilon = read_input_file()

# Ответ: -6 2 7 -1


print("Метод простых итераций:")
x_simple, iter_simple = simple_iterations_solve(A, b, epsilon)
print(f"Решение: {x_simple}")
print(f"Количество итераций: {iter_simple}")

print("\nМетод Зейделя:")
x_seidel, iter_seidel = zeidel_solve(A, b, epsilon)
print(f"Решение: {x_seidel}")
print(f"Количество итераций: {iter_seidel}")
