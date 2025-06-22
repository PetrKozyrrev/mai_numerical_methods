task_function_5 = lambda x : (x**2)/(x**3 - 27)
X0 = -2
Xk = 2
h1 = 1.0
h2 = 0.5

def rectangle_method(x0, xk, h):
    intg = 0
    xi1 = x0
    xi2 = x0 + h
    while (xi2 <= xk):
        intg += task_function_5((xi1 + xi2)/2)
        xi1 += h
        xi2 += h
    return intg*h

print("Метод прямоугольника")
print(f"c шагом {h1}: ", rectangle_method(X0, Xk, h1))
print(f"c шагом {h2}: ", rectangle_method(X0, Xk, h2))


def trapezoida_method(x0, xk, h):
    n = int((xk - x0) / h)
    intg = 0
    intg += task_function_5(x0) / 2 + task_function_5(xk) / 2
    for i in range(1, n):
        xi = x0 + i * h
        intg += task_function_5(xi)

    return intg * h

print("Метод трапеций")
print(f"c шагом {h1}: ", trapezoida_method(X0, Xk, h1))
print(f"c шагом {h2}: ", trapezoida_method(X0, Xk, h2))


def simpson(x0, xk, h, f):
    n = int((xk - x0) / h)

    if n % 2 != 0:
        raise ValueError("Количество интервалов должно быть четным. Попробуйте изменить шаг h.")

    intg = f(x0) + f(xk)
    for i in range(1, n):
        x = x0 + i * h
        if i % 2 == 0:
            intg += 2 * f(x)
        else:
            intg += 4 * f(x)

    return intg * h / 3

print("Метод Симпсона")
print(f"c шагом {h1}: ", simpson(X0, Xk, h1, task_function_5))
print(f"c шагом {h2}: ", simpson(X0, Xk, h2, task_function_5))

def runge_romberg_method(x0, xk, h, r, num):
    if num == 1:
        p = 1
        return (rectangle_method(x0, xk, h) - rectangle_method(x0, xk, h * r)) / (r ** p - 1)
    elif num == 2:
        p = 2
        return (trapezoida_method(x0, xk, h) - trapezoida_method(x0, xk, h * r)) / (r ** p - 1)
    elif num == 3:
        p = 4
        return (simpson(x0, xk, h, task_function_5) - simpson(x0, xk, h * r,task_function_5)) / (r ** p - 1)


print("\n===== Погрешность Метод Рунге-Ромберга =====")
print("\nМетод прямоугольников:")
print(f"Для шага {h1}: ", runge_romberg_method(X0, Xk, h1, 2, 1))
print(f"Для шага {h2}: ", runge_romberg_method(X0, Xk, h2, 2, 1))

print("\nМетод трапеций:")
print(f"Для шага {h1}: ", runge_romberg_method(X0, Xk, h1, 2, 2))
print(f"Для шага {h2}: ", runge_romberg_method(X0, Xk, h2, 2, 2))

print("\nМетод Симпсона:")
print(f"Для шага {h1}: ", runge_romberg_method(X0, Xk, h1, 2, 3))
print(f"Для шага {h2}: ", runge_romberg_method(X0, Xk, h2, 2, 3))

print("\n===== Уточнение значений Метод Рунге-Ромберга =====")
print("\nМетод прямоугольников:")
print(f"Для шага {h1}: ", rectangle_method(X0, Xk, h1) + runge_romberg_method(X0, Xk, h1, 2, 1))
print(f"Для шага {h2}: ", rectangle_method(X0, Xk, h2) + runge_romberg_method(X0, Xk, h2, 2, 1))

print("\nМетод трапеций:")
print(f"Для шага {h1}: ", trapezoida_method(X0, Xk, h1) + runge_romberg_method(X0, Xk, h1, 2, 2))
print(f"Для шага {h2}: ", trapezoida_method(X0, Xk, h2) + runge_romberg_method(X0, Xk, h2, 2, 2))

print("\nМетод Симпсона:")
print(f"Для шага {h1}: ", simpson(X0, Xk, h1, task_function_5) + runge_romberg_method(X0, Xk, h1, 2, 3))
print(f"Для шага {h2}: ", simpson(X0, Xk, h2, task_function_5) + runge_romberg_method(X0, Xk, h2, 2, 3))
