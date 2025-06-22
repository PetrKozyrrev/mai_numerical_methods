arr = []

with open("1.2.txt", 'r') as f:
    for i in f.readlines():
        arr.append(list(map(int, i.split(' '))))


def checking(arr):
    n = len(arr)
    for i in range(1, n-1):
        if arr[i][0] == 0 or arr[i][2] == 0:
            return False

    flag = 0
    for i in range(n):
        if i == 0:
            if abs(arr[i][0]) == abs(arr[i][1]):
                flag = 1
            if abs(arr[i][0]) < abs(arr[i][1]):
                return False
        elif i == n - 1:
            if abs(arr[i][1]) == abs(arr[i][0]):
                flag = 1
            if abs(arr[i][1]) < abs(arr[i][0]):
                return False
        else:
            if abs(arr[i][1]) == (abs(arr[i][0]) + abs(arr[i][2])):
                flag = 1
            if abs(arr[i][1]) < (abs(arr[i][0]) + abs(arr[i][2])):
                return False

    if(flag):
        return True
    else:
        return False


def method_progonki(arr):
    if checking(arr):
        print("Проверка: Достаточное условие выполнено")
    else:
        print("Проверка: Достаточное условие НЕ выполнено")

    n = len(arr)

    p = [0] * n
    q = [0] * n

    p[0] = ((-1) * arr[0][1]) / arr[0][0]
    q[0] = arr[0][2] / arr[0][0]

    # прямой ход
    for i in range(1, n):
        if i != n - 1:
            p[i] = ((-1) * arr[i][2]) / (arr[i][1] + arr[i][0] * p[i - 1])
            q[i] = (arr[i][3] - arr[i][0] * q[i - 1]) / (arr[i][1] + arr[i][0] * p[i - 1])
        else:
            p[i] = 0
            q[i] = (arr[i][2] - arr[i][0] * q[i - 1]) / (arr[i][1] + arr[i][0] * p[i - 1])

    # обратный ход
    x = [0] * n
    x[n - 1] = q[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = p[i] * x[i + 1] + q[i]

    return x

print("Решение методом прогонки: ") # 7 3 -7 5 -9
print(method_progonki(arr))
