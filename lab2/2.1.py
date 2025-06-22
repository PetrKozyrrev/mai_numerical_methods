import math

def simple_iterations(phi, phi_d, eps, start_a, start_b):
    if (f(start_a) * f(start_b) >= 0):
        raise ValueError("f(a) * f(b) должен быть < 0")
    q = -1
    if (abs(phi_d(start_a)) >= 1):
        if (abs(phi_d(start_b)) >= 1):
            raise ValueError("Достаточное условие не выполнено")
        else:
            q = phi_d(start_b)
    else:
        q = phi_d(start_a)

    iter_count = 0
    coef = q / (1 - q)
    x_k = start_b
    dx = 10e9
    while eps < coef * dx:
        x_k_next = phi(x_k)
        dx = abs(x_k_next - x_k)
        x_k = x_k_next
        iter_count += 1
    return x_k, iter_count


def newton_method(f, f_d, f_d2, eps, start_a, start_b):
    if (f(start_a) * f(start_b) >= 0):
        raise ValueError("f(a) * f(b) должен быть < 0")

    if (f_d(start_a) * f_d(start_b) < 0):
        raise ValueError("f`(x) не знакопостоянна")

    if (f_d2(start_a) * f_d2(start_b) < 0):
        raise ValueError("f``(x) не знакопостоянна")

    x0 = start_a
    if abs(f(start_a) * f_d2(start_a)) > (f_d(start_a)) ** 2:
        if abs(f(start_b) * f_d2(start_b) > (f_d(start_b)) ** 2):
            raise ValueError("abs(f(a) * f``(a)) должно быть <= f`(a)**2 или abs(f(b) * f``(b)) < (f`(b))**2 должен быть <= (f`(b))**2")
        else:
            x0 = start_b
    else:
        if abs(f(start_b) * f_d2(start_b) <= (f_d(start_b)) ** 2):
            if (abs(f(start_a)) <= abs(f(start_b))):
                x0 = start_a
            else:
                x0 = start_b
        else:
            x0 = start_a
    iter_count = 0
    x_k = x0
    dx = 10e9
    while eps < dx:
        x_k_next = x_k - f(x_k) / f_d(x_k)
        dx = abs(x_k_next - x_k)
        x_k = x_k_next
        iter_count += 1
    return x_k, iter_count


eps = 0.00000001

f = lambda x: math.log(x + 1) - 2 * x + 0.5
f_d = lambda x: 1 / (x + 1) - 2
f_d2 = lambda x: -1 / (x + 1) ** 2

phi = lambda x: (math.log(x + 1) + 0.5) / 2
phi_d = lambda x: 1 / (2 * x + 2)

start_a, start_b = 0, 2

x, iter_count = simple_iterations(phi, phi_d, eps, start_a, start_b)
print(f"Eps = {eps}")
print(f"Simple iterations | x={x} | iter_count={iter_count}")
print("--------------------------------------------------------")

x, iter_count = newton_method(f, f_d, f_d2, eps, start_a, start_b)
print(f"Newton method:    | x={x} | iter_count={iter_count}")
