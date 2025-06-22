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

def copy(A):
    return [[A[i][j] for j in range(n)] for i in range(n)]

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

    # 1. Ly = b
    y = forward_substitution(L, b)

    # 2.Ux = y
    x = backward_substitution(U, y)

    return x


def det(A):
    n = len(A)
    LU = lu_decomposition(A)
    d = 1
    for i in range(n):
        d *= LU[i][i]
    return d


def inverse(A):
    n = len(A)
    M = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        b = [1 if j == i else 0 for j in range(n)]
        x = solve_system(A, b)
        for k in range(n):
            M[k][i] = x[k]
    return M


# Чтение матрицы из файла
def read_input_file():
    with open('1.1.txt', 'r') as f:
        matrix = [list(map(int, row.split())) for row in f.readlines()]
    return matrix


matrix = read_input_file()
n = len(matrix)
A = []
b = []
for i in range(n):
    A.append(matrix[i][:len(matrix)])
    b.append((matrix[i][-1]))

# Решаем систему
x = solve_system(A, b)

# Выводим результат
print("Решение системы:") # -2 -3 -6 -7
print(x)

print("Определитель: ")  # -345
print(det(A))

print("Обратная матрица:")
M = inverse(A)
for i in M:
    print(*i)

print("Проверка A*A^-1: ")
E = multiply_matrices(A, M)
for i in range(n):
    for j in range(n):
        print(f"{round(E[i][j],2)}",end=' ')
    print('\n')
