def euler_method(f, x0, y0, x_end, h):
    """
    Метод Эйлера для решения ОДУ.

    f: Правая часть ОДУ (dy/dx = f(x, y)).
    x0: Начальная точка.
    y0: Начальное значение y.
    x_end: Конец интервала.
    h: Шаг интегрирования.
    """
    x_values = [x0]
    y_values = [y0]

    x, y = x0, y0
    while x < x_end:
        y = y + h * f(x, y)
        x = x + h
        if x<x_end:
          x_values.append(x)
          y_values.append(y)

    return x_values, y_values

# Правая часть ОДУ
def f(x, y):
    return -2 * x + y

# Начальные условия
x0 = 0
y0 = 1
x_end = 1
h = 0.1

# Решение методом Эйлера
x_vals, y_vals = euler_method(f, x0, y0, x_end, h)

# Вывод результатов
print("Решение:")
for x, y in zip(x_vals, y_vals):
    print(f"x = {x:.1f}, y = {y:.4f}")


import matplotlib.pyplot as plt
import numpy as np

# Точное решение для сравнения (если известно)
# def exact_solution(x):
#     return 3 * np.exp(x) - 2 * x - 2

# Построение графика
plt.plot(x_vals, y_vals, 'o-', label="Метод Эйлера")
# x_exact = np.linspace(x0, x_end, 100)
# y_exact = exact_solution(x_exact)
# plt.plot(x_exact, y_exact, label="Точное решение")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Сравнение численного и точного решения")
plt.legend()
plt.grid()
plt.show()


def predictor_corrector(f, x0, y0, x_end, h):
    """
    Метод предиктора-корректора для решения ОДУ.

    f: Правая часть ОДУ (dy/dx = f(x, y)).
    x0: Начальная точка.
    y0: Начальное значение y.
    x_end: Конец интервала.
    h: Шаг интегрирования.
    """
    x_values = [x0]
    y_values = [y0]

    x, y = x0, y0
    while x < x_end:
        # Шаг предиктора (метод Эйлера)
        y_pred = y + h * f(x, y)

        # Шаг корректора (метод трапеций)
        x_next = x + h
        y_corr = y + h / 2 * (f(x, y) + f(x_next, y_pred))

        # Обновление значений
        x, y = x_next, y_corr
        x_values.append(x)
        y_values.append(y)

    return x_values, y_values

# Правая часть ОДУ
def f(x, y):
    return x + y

# Начальные условия
x0 = 0
y0 = 1
x_end = 1
h = 0.1

# Решение методом предиктора-корректора
x_vals, y_vals = predictor_corrector(f, x0, y0, x_end, h)

# Вывод результатов
print("Решение:")
for x, y in zip(x_vals, y_vals):
    print(f"x = {x:.1f}, y = {y:.4f}")


import matplotlib.pyplot as plt
import numpy as np

# Точное решение (если известно)
def exact_solution(x):
    return -x - 1 + 2 * np.exp(x)

# Построение графика
plt.plot(x_vals, y_vals, 'o-', label="Метод предиктора-корректора")
x_exact = np.linspace(x0, x_end, 100)
y_exact = exact_solution(x_exact)
plt.plot(x_exact, y_exact, label="Точное решение")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Сравнение численного и точного решения")
plt.legend()
plt.grid()
plt.show()


def dydx(x, y):
  return (-2*x + y)

def rungeKutta(x0, y0, x, h):
  n = round((x-x0)/h)
  y = y0

  for i in range(1, n + 1):
    k1 = h * dydx(x0, y)
    k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)

    y = y + (1/6) * (k1 + 2 * k2)
    x0 = x0 + h
  return y

x0 = 0
y = 1
x = 1
h = 0.1
print('y(x)= ', rungeKutta(x0, y, x, h))

def runge_kutta_4(f, x0, y0, x_end, h):
    """
    Метод Рунге-Кутты 4-го порядка для решения ОДУ.

    f: Правая часть ОДУ (dy/dx = f(x, y)).
    x0: Начальная точка.
    y0: Начальное значение y.
    x_end: Конец интервала.
    h: Шаг интегрирования.
    """
    x_values = [x0]
    y_values = [y0]

    x, y = x0, y0
    while x < x_end:
        k1 = f(x, y)
        k2 = f(x + h / 2, y + h / 2 * k1)
        k3 = f(x + h / 2, y + h / 2 * k2)
        k4 = f(x + h, y + h * k3)

        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        x = x + h

        x_values.append(x)
        y_values.append(y)

    return x_values, y_values

# Правая часть ОДУ
def f(x, y):
    return -2 * x + y

# Начальные условия
x0 = 0
y0 = 1
x_end = 1
h = 0.1

# Решение методом Рунге-Кутты 4-го порядка
x_vals, y_vals = runge_kutta_4(f, x0, y0, x_end, h)

# Вывод результатов
print("Решение:")
for x, y in zip(x_vals, y_vals):
    print(f"x = {x:.1f}, y = {y:.4f}")


import numpy as np

def runge_kutta(f, x0, y0, x_end, h, order):
    """
    Универсальная реализация метода Рунге-Кутты для решения ОДУ.

    f: Функция правой части ОДУ dy/dx = f(x, y).
    x0: Начальная точка.
    y0: Начальное значение y.
    x_end: Конечная точка x.
    h: Шаг интегрирования.
    order: Порядок метода Рунге-Кутты (1, 2, 3 или 4).
    """
    if order not in [1, 2, 3, 4]:
        raise ValueError("Порядок метода должен быть 1, 2, 3 или 4.")

    # Список для хранения значений x и y
    x_values = [x0]
    y_values = [y0]

    # Текущие значения x и y
    x, y = x0, y0

    while x < x_end:
        if order == 1:  # Метод Эйлера (Рунге-Кутты 1-го порядка)
            k1 = f(x, y)
            y = y + h * k1

        elif order == 2:  # Метод Рунге-Кутты 2-го порядка
            k1 = f(x, y)
            k2 = f(x + h / 2, y + h / 2 * k1)
            y = y + h * k2

        elif order == 3:  # Метод Рунге-Кутты 3-го порядка
            k1 = f(x, y)
            k2 = f(x + h / 2, y + h / 2 * k1)
            k3 = f(x + h, y - h * k1 + 2 * h * k2)
            y = y + h / 6 * (k1 + 4 * k2 + k3)

        elif order == 4:  # Метод Рунге-Кутты 4-го порядка
            k1 = f(x, y)
            k2 = f(x + h / 2, y + h / 2 * k1)
            k3 = f(x + h / 2, y + h / 2 * k2)
            k4 = f(x + h, y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        x = x + h
        x_values.append(x)
        y_values.append(y)

    return np.array(x_values), np.array(y_values)

# Пример использования
def example_function(x, y):
    return -2 * x + y  # Пример ОДУ: dy/dx = -2x + y

# Начальные условия
x0 = 0
y0 = 1
x_end = 1
h = 0.1

# Решение методом Рунге-Кутты разных порядков
for order in range(1, 5):
    x_vals, y_vals = runge_kutta(example_function, x0, y0, x_end, h, order)
    print(f"Метод Рунге-Кутты {order}-го порядка:")
    for x, y in zip(x_vals, y_vals):
        print(f"x = {x:.2f}, y = {y:.4f}")
    print()
