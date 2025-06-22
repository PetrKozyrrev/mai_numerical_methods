Xi = [0.2, 0.5, 0.8, 1.1, 1.4]
Yi = [12.906, 5.5273, 3.8777, 3.2692, 3.0319]
XX = 0.8


def approximation_left_right(x, X, Y):
    ind = -1
    for i in range(len(X)):
        if X[i] > x:
            ind = i - 1
            break

    p_left = (Y[ind] - Y[ind - 1]) / (X[ind] - X[ind - 1])
    p_right = (Y[ind + 1] - Y[ind]) / (X[ind + 1] - X[ind])

    return p_left, p_right, ind


def P1(x, X, Y):
    p_left, p_right, ind = approximation_left_right(x, X, Y)
    return p_left + (p_right - p_left) * (2 * x - X[ind - 1] - X[ind]) / (X[ind + 1] - X[ind - 1])


def P2(x, X, Y):
    p_left, p_right, ind = approximation_left_right(x, X, Y)
    return (p_right - p_left) * 2 / (X[ind + 1] - X[ind - 1])

print(f"Левосторонняя производная: {approximation_left_right(XX, Xi, Yi)[0]}")
print(f"Правосторонняя производная: {approximation_left_right(XX, Xi, Yi)[1]}")
print("Первая производная в точке X* : ", P1(XX, Xi, Yi))
print("Вторая производная в точке X* : ", P2(XX, Xi, Yi))
