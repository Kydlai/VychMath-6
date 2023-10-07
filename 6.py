# Лабораторная работа #6 (5)
# @ Кудлай Никита, PЗ2141
# Усовершенствованный метод Эйлера
# Метод Рунге-Кутта 4 порядка
# Метод Адамса

import numpy as np
import matplotlib.pyplot as plt
from math import exp

y0 = 0
def mod_euler_method(f, a, b, y0, h, E):

    dots = [(a, y0)]
    n = int((b - a) / h)
    for i in range(1, n + 1):
        dots.append((dots[-1][0] + h, dots[-1][1] + h / 2 * (f(dots[-1][0], dots[-1][1]) +
                                                                f(dots[-1][0] + h, dots[-1][1] + h * f(dots[-1][0], dots[-1][1])))))
    return dots





def runge_kutte_method(f, a, b, y0, h, E):
    dots = [(a, y0)]
    n = int((b - a) / h)
    for i in range(1, n + 1):
        k1 = h * f(dots[i-1][0], dots[i-1][1])
        k2 = h * f(dots[i-1][0] + h / 2, dots[i-1][1] + k1 / 2)
        k3 = h * f(dots[i-1][0] + h / 2, dots[i-1][1] + k2 / 2)
        k4 = h * f(dots[i-1][0] + h, dots[i-1][1] + k3)
        dots.append((dots[i-1][0] + h,
                    dots[i-1][1] + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)))
    return dots
def adams_method(f, a, b, y0, h, E):
    n = int((b - a) / h)
    b0 = min(b, a + 4 * h)
    dots = mod_euler_method(f, a, b0, y0, h, E)[:4]
    for i in range(4, n + 1):
        df = f(dots[-1][0], dots[-1][1]) - f(dots[-2][0], dots[-2][1])
        d2f = f(dots[-1][0], dots[-1][1]) - 2 * f(dots[-2][0], dots[-2][1]) + f(dots[-3][0], dots[-3][1])
        d3f = f(dots[-1][0], dots[-1][1]) - 3 * f(dots[-2][0], dots[-2][1]) + \
            3 * f(dots[-3][0], dots[-3][1]) - f(dots[-4][0], dots[-4][1])
        dots.append((dots[-1][0] + h,
                     dots[-1][1] + h * f(dots[-1][0], dots[-1][1]) +
                     (h ** 2) * df / 2 + 5 * (h ** 3) * d2f / 12 + 3 * (h ** 4) * d3f / 8))

    return dots


def plot(x, y, acc_x, acc_y, flag):

    ax = plt.gca()
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.plot(1, 0, marker=">", ms=5, color='k',
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.plot(0, 1, marker="^", ms=5, color='k',
            transform=ax.get_xaxis_transform(), clip_on=False)

    plt.plot(x, y, label="y(x)")
    if flag:
        plt.plot(acc_x, acc_y, label="Точный y(x)")

    plt.legend()
    plt.show(block=False)


def gettask(task_id):
    if task_id == '1':
        return lambda x, y: y + (1 + x) * (y ** 2), \
               lambda x: -1 / x, \
               1, \
               1.5, \
               -1
    elif task_id == '2':
        return lambda x, y: (x ** 2) - 2 * y, \
               lambda x: 0.75 * exp(-2 * x) + 0.5 * (x ** 2) - 0.5 * x + 0.25, \
               0, \
               1, \
               1
    elif task_id == '3':
        return lambda x, y: 2 * (x ** 2) + x - y, \
               lambda x: (-1) / exp(x) + 2 * (x ** 2) - 3 * x + 3, \
               0, \
               2, \
               2
    else:
        return None


def getdata_input():
    data = {}

    print("\nВыберите метод дифференцирования.")
    print(" 1 — Модифицированный метод Эйлера")
    print(" 2 - Метод Рунге-Кутта 4 порядка")
    print(" 3 — Метод Адамса")
    while True:
        try:
            method_id = input("Метод дифференцирования: ")
            if method_id != '1' and method_id != '2' and method_id != '3':
                raise AttributeError
            break
        except AttributeError:
            print("Метода нет в списке!")
    data['method_id'] = method_id

    print("\nВыберите задачу.")
    print(" 1 — y' = y + (1 + x)y²\n     на [1; 1,5] при y(1) = -1")
    print(" 2 — y' = x² - 2y\n     на [0; 1] при y(0) = 1")
    print(" 3 — y' = 2x² + x - y\n     на [0; 2] при y(0) = 2")
    while True:
        try:
            task_id = input("Задача: ")
            func, acc_func, a, b, y0 = gettask(task_id)
            if func is None:
                raise AttributeError
            break
        except AttributeError:
            print("Функции нет в списке!")
    data['f'] = func
    data['acc_f'] = acc_func
    data['flag'] = True
    data['a'] = a
    data['b'] = b
    data['y0'] = y0

    print("\nЖелаете изменить границы и начальное значение? Введите 1 если да")
    while True:
        try:
            flag = int(input())
            if flag != 1:
                flag = 0
            break
        except (ValueError, ArithmeticError):
            flag = 0
    if flag == 1:
        data['flag'] = False
        print("\nВведите границы a и b.")
        while True:
            try:
                a_b = [float(i) for i in input("Границы a и b: ").split()]
                if a_b[0] >= a_b[1]:
                    raise ArithmeticError
                break
            except (ValueError, ArithmeticError):
                print("Введите корректные значения a и b.")
        data['a'] = a_b[0]
        data['b'] = a_b[1]
        print("\nВведите значение y в точке a =", a_b[0])
        while True:
            try:
                y0 = float(input("Значение y: "))
                break
            except (ValueError, ArithmeticError):
                print("Введите корректное значение y.")
        data['y0'] = y0

    print("\nВведите шаг точек.")
    while True:
        try:
            h = float(input("Шаг точек: "))
            if h <= 0 or h > 1:
                raise ArithmeticError
            break
        except (ValueError, ArithmeticError):
            print("Шаг точек должен быть положительным числом.")
    data['h'] = h

    print("\nВведите желаемую точность (0 если не проверять на точность)")
    while True:
        try:
            E = float(input("Желаемая точность: "))
            if E < 0:
                raise ArithmeticError
            if E == 0:
                E = 100
            break
        except (ValueError, ArithmeticError):
            print("Точность должна быть положительным числом.")
    data['E'] = E

    return data


def main():
    print("\tЛабораторная работа #6 (5)")
    print("\tЧисленное дифференцирование")
    R = 0
    data = getdata_input()
    euler = mod_euler_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
    runge_kutte = runge_kutte_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
    adams = adams_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
    if data['method_id'] == '1':
        while True:
            answer = mod_euler_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            answer_1 = mod_euler_method(data['f'], data['a'], data['b'], data['y0'], data['h']/2, data['E'])
            euler = answer
            runge_kutte = runge_kutte_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            adams = adams_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            R = (answer[-1][1]-answer_1[-1][1])/(2**2-1)
            if abs(R) <= data['E']:
                break
            else:
                data['h'] /= 2
    elif data['method_id'] == '2':
        while True:
            answer = runge_kutte_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            answer_1 = runge_kutte_method(data['f'], data['a'], data['b'], data['y0'], data['h']/2, data['E'])
            euler = mod_euler_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            runge_kutte = runge_kutte_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            adams = adams_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            R = (answer[-1][1]-answer_1[-1][1])/(2**4-1)
            if abs(R) <= data['E']:
                break
            else:
                data['h'] /= 2

    elif data['method_id'] == '3':
        while True:
            answer = adams_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            euler = mod_euler_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            runge_kutte = runge_kutte_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            adams = adams_method(data['f'], data['a'], data['b'], data['y0'], data['h'], data['E'])
            R = 0
            for i in range(len(answer)):
                R = max(R, abs(answer[i][1] - data['acc_f'](answer[i][0])))
            if abs(R) <= data['E']:
                break
            else:
                data['h'] /= 2
    else:
        answer = None

    if answer is None:
        print("\n\nВо время вычисления произошла ошибка!")
    else:

        x = np.array([dot[0] for dot in answer])
        y = np.array([dot[1] for dot in answer])
        acc_x = np.linspace(np.min(x), np.max(x), 100)
        acc_y = [data['acc_f'](i) for i in acc_x]
        plot(x, y, acc_x, acc_y, data['flag'])

        print("\n\nИтоговый шаг h = ", data['h'], "\nТочность R = ", "%.6f" % R)
        print("\nРезультаты вычисления.")
        if data['flag']:
            print("%12s%12s%12s%12s%12s" % ("x", "Эйлер y", "Рунге y", "Адамс y", "Точный y"))
            for i in range(len(answer)):
                print("%12.4f%12.4f%12.4f%12.4f%12.4f" % (answer[i][0], euler[i][1], runge_kutte[i][1], adams[i][1], data['acc_f'](answer[i][0])))
        else:
            print("%12s%12s%12s%12s" % ("x", "Эйлер y", "Рунге y", "Адамс y"))
            for i in range(len(answer)):
                print("%12.4f%12.4f%12.4f%12.4f" % (answer[i][0], euler[i][1], runge_kutte[i][1], adams[i][1]))

    input("\n\nНажмите Enter, чтобы выйти.")


main()
