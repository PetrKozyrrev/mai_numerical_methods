import math

# Произведение матрицы на вектор
def multiply_matrix_vector(A, v):
    if len(A[0]) != len(v):
        raise ValueError("Количество столбцов в матрице должно быть равно длине вектора.")

    result = [0] * len(A)

    for i in range(len(A)):
        result[i] = sum(A[i][j] * v[j] for j in range(len(A[0])))  # Скалярное произведение строки с вектором

    return result


# Перемножение матриц
def multiply_matrices(A, B):
    if len(A[0]) != len(B):
        raise ValueError("Количество столбцов в первой матрице должно быть равно количеству строк во второй матрице.")

    result = [[0] * len(B[0]) for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result


# Транспонирование матрицы
def transpose(A):
    n = len(A)
    transposed = [[A[j][i] for j in range(n)] for i in range(n)]
    return transposed


# Чтение матрицы из файла
def read_input_file():
    with open('1.4.txt', 'r') as f:
        n = int(f.readline())
        A = []
        for i in range(n):
            A.append(list(map(int, f.readline().split())))
        epsilon = float(f.readline().strip())
    return n, A, epsilon

# Норма матрицы
def mat_norm(a):
    res = 0
    n = len(a)
    for i in range(n):
        for j in range(i + 1, n):
            res += a[i][j] ** 2
    return math.sqrt(res)

# Координаты максимума
def get_coors_of_max(a):
    n = len(a)
    i_max = 0
    j_max = 1
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if abs(a[i][j]) > abs(a[i_max][j_max]):
                i_max, j_max = i, j
    return i_max, j_max

# Вычисление угла phi
def calc_phi(a, i, j):
    if a[i][i] == a[j][j]:
        return math.pi / 4
    return 0.5 * math.atan2((2 * a[i][j]), (a[i][i] - a[j][j]))

# Единичная матрица NxN
def eye(n):
    u = []
    for k in range(n):
        tmp = []
        for t in range(n):
            if (t == k):
                tmp.append(1)
            else:
                tmp.append(0)
        u.append(tmp)
    return u

# Матрица поворота
def get_rotation_matrix(n, phi, i, j):
    u = eye(n)
    u[i][i] = math.cos(phi)
    u[i][j] = -math.sin(phi)
    u[j][i] = math.sin(phi)
    u[j][j] = math.cos(phi)
    return u


def jacobi_rotations(a, eps):
    if(check_symmetry(a) == False):
        raise ValueError("Матрица НЕ симметрична")

    n = len(a)
    a_k = [row[:] for row in a]
    u = eye(n)
    iter_count = 0
    while mat_norm(a_k) > eps:
        iter_count += 1
        i, j = get_coors_of_max(a_k)
        phi = calc_phi(a_k, i, j)
        u_k = get_rotation_matrix(n, phi, i, j)
        u = multiply_matrices(u, u_k)
        tmp = multiply_matrices(transpose(u_k), a_k)
        a_k = multiply_matrices(tmp, u_k)
        print(f"iter = {iter_count} mat_norm = {mat_norm(a_k)} eps = {eps}")
    return a_k, u, iter_count

# Проверка симметричности матрицы
def check_symmetry(a):
    for i in range(n):
        for j in range(i+1,n):
            if(a[i][j] != a[j][i]):
                return False
    return True

n, a, eps = read_input_file()

a_k, u, iter_count = jacobi_rotations(a, eps)
print("Собственные значения: ") # -7.534 6.123 8.411
for i in range(n):
    print(f" l_{i + 1} = {a_k[i][i]}")
print("Собственные векторы: ")
u = transpose(u)
for i in range(n):
    print(f" x_{i + 1} = {u[i]}")

print("Количество итераций:", iter_count)


print("\nПроверка (Ax = lx):")
for i in range(n):
    print(f"{multiply_matrix_vector(a, u[i])} = {[g * a_k[i][i] for g in u[i]]}")

