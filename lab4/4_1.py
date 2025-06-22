from math import sqrt, exp, cos, tan
import matplotlib.pyplot as plt


def frange(start, stop, step):
    while start <= stop + 1e-10:
        yield round(start, 10)
        start += step


def f(x, y, z):
    return 2 * tan(x) * z + 3 * y


def g(x, y, z):
    return z


def exact_solution(x):
    return 0.25 * exp(-sqrt(2) * x) * ((2 + 3 * sqrt(2)) * exp(2 * sqrt(2) * x) + 2 - 3 * sqrt(2)) / cos(x)


def euler_method(f, g, y0, z0, interval, h):
    l, r = interval
    x = list(frange(l, r, h))
    y = [y0]
    z = z0
    for i in range(len(x) - 1):
        z += h * f(x[i], y[i], z)
        y.append(y[i] + h * g(x[i], y[i], z))
    return x, y


def runge_kutta_method(f, g, y0, z0, interval, h, return_z=False):
    l, r = interval
    x = list(frange(l, r, h))
    y = [y0]
    z = [z0]
    for i in range(len(x) - 1):
        K1 = h * g(x[i], y[i], z[i])
        L1 = h * f(x[i], y[i], z[i])
        K2 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * K1, z[i] + 0.5 * L1)
        L2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * K1, z[i] + 0.5 * L1)
        K3 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * K2, z[i] + 0.5 * L2)
        L3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * K2, z[i] + 0.5 * L2)
        K4 = h * g(x[i] + h, y[i] + K3, z[i] + L3)
        L4 = h * f(x[i] + h, y[i] + K3, z[i] + L3)
        delta_y = (K1 + 2 * K2 + 2 * K3 + K4) / 6
        delta_z = (L1 + 2 * L2 + 2 * L3 + L4) / 6
        y.append(y[i] + delta_y)
        z.append(z[i] + delta_z)

    return (x, y) if not return_z else (x, y, z)


def adams_method(f, g, y0, z0, interval, h):
    x_runge, y_runge, z_runge = runge_kutta_method(f, g, y0, z0, interval, h, return_z=True)
    x = x_runge
    y = y_runge[:4]
    z = z_runge[:4]
    for i in range(3, len(x_runge) - 1):
        z_i = z[i] + h * (55 * f(x[i], y[i], z[i]) -
                          59 * f(x[i - 1], y[i - 1], z[i - 1]) +
                          37 * f(x[i - 2], y[i - 2], z[i - 2]) -
                          9 * f(x[i - 3], y[i - 3], z[i - 3])) / 24
        z.append(z_i)
        y_i = y[i] + h * (55 * g(x[i], y[i], z[i]) -
                          59 * g(x[i - 1], y[i - 1], z[i - 1]) +
                          37 * g(x[i - 2], y[i - 2], z[i - 2]) -
                          9 * g(x[i - 3], y[i - 3], z[i - 3])) / 24
        y.append(y_i)
    return x, y


def runge_rombert_method(h1, h2, y1, y2, p):
    assert h1 == h2 * 2
    norm = 0
    for i in range(len(y1)):
        norm += (y1[i] - y2[i * 2]) ** 2
    return norm ** 0.5 / (2**p + 1)


def mae(y1, y2):
    assert len(y1) == len(y2)
    return sum(abs(a - b) for a, b in zip(y1, y2)) / len(y1)

y0 = 1
dy0 = 3
interval = (0, 1)
h = 0.1

x_euler, y_euler = euler_method(f, g, y0, dy0, interval, h)
plt.plot(x_euler, y_euler, label=f'euler, h={h}')
x_euler2, y_euler2 = euler_method(f, g, y0, dy0, interval, h / 2)
plt.plot(x_euler2, y_euler2, label=f'euler, h={h/2}')

x_runge, y_runge = runge_kutta_method(f, g, y0, dy0, interval, h)
plt.plot(x_runge, y_runge, label=f'runge-kutta, h={h}')
x_runge2, y_runge2 = runge_kutta_method(f, g, y0, dy0, interval, h / 2)
plt.plot(x_runge2, y_runge2, label=f'runge-kutta, h={h/2}')

x_adams, y_adams = adams_method(f, g, y0, dy0, interval, h)
plt.plot(x_adams, y_adams, label=f'adams, h={h}')
x_adams2, y_adams2 = adams_method(f, g, y0, dy0, interval, h / 2)
plt.plot(x_adams2, y_adams2, label=f'adams, h={h/2}')

x_exact = list(frange(interval[0], interval[1], h))
x_exact2 = list(frange(interval[0], interval[1], h / 2))
y_exact = [exact_solution(xi) for xi in x_exact]
y_exact2 = [exact_solution(xi) for xi in x_exact2]
plt.plot(x_exact, y_exact, label='real')

print("Погрешность численного решения путем сравнения с точным решением:")
print(f"Шаг = {h}")
print("Эйлера:", mae(y_euler, y_exact))
print("Рунге-Кутта:", mae(y_runge, y_exact))
print("Адамса:", mae(y_adams, y_exact))

print(f"\nШаг = {h/2}")
print("Эйлера:", mae(y_euler2, y_exact2))
print("Рунге-Кутта:", mae(y_runge2, y_exact2))
print("Адамса:", mae(y_adams2, y_exact2))

print("\nПогрешность методом Рунге – Ромберга:")
print("Эйлера:", runge_rombert_method(h, h/2, y_euler, y_euler2, 1))
print("Рунге-Кутта:", runge_rombert_method(h, h/2, y_runge, y_runge2, 4))
print("Адамса:", runge_rombert_method(h, h/2, y_adams, y_adams2, 4))

plt.legend()
plt.grid()
plt.show()
