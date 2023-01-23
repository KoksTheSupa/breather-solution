from decimal import *
import _thread
from math import fabs, sqrt, isnan
from numba import prange

import matplotlib.pyplot as plt

N = 100
DJ = Decimal("0.16")
SQDJ = Decimal(sqrt(1 + DJ**2))
OMEGA = Decimal("-0.28")
B = Decimal("0.15")
EPS = H = Decimal("0.000001")
S = [Decimal(0.0)] * (N+1)
DELIMITER = 2_000_000_000  # double results[2000000000]
s0 = Decimal("0.509457626013908")
temp = Decimal(Decimal("0.4") / (DELIMITER))


"""Построение графиков"""


def plot(x, y):
    """x = [i for i in prange(x)]
    xmax = x[y.index(max(y))]"""
    fig, graph = plt.subplots()
    graph.scatter(x, y, 30, 'black')
    graph.plot(x, y, linewidth=1.0)
    """graph.spines['bottom'].set_position('center') """
    graph.set_ylabel('S')
    graph.set_xlabel('Omega')
    """graph.annotate(f'Local max = {max(y)}', xy=(xmax, max(y)), xytext=(xmax, max(y) + Decimal(0.10)),
                   arrowprops=dict(facecolor='black',
                                   shrink=0.05, headlength=5, width=2),
                   )
    graph.set_ylim(-max(y) - Decimal('0.15'), max(y) + Decimal('0.15'))"""
    graph.grid(True)
    plt.show()


"""Запись в файл"""


def writef(S: list):
    with open(f'e_N({N+1})_Omega{OMEGA}_B{B}_dj{DJ*100}_h{H}.nb', 'a+') as file:
        file.write(
            f'text=\"parameters: omega={OMEGA}, N(real)=({N+1}), B={B}, dj={DJ}, h={H}, S[N]={S[0]} (EDGE algorithm)\"\ndataS = {{')
        for i in range(N):

            file.write('{%s,%s},' % (i, (S[i])))

        file.write('{%s,%s}};' % (N, S[N]))
        file.write(
            "\n\nListPlot[dataS,Joined->True,Mesh->All,PlotRange->All]\n\n")


"""Расчет цепочки"""


def counter(start: Decimal, N: int, H=Decimal("0.0")):
    S[0] = Decimal(str(start)).quantize(
        Decimal("1.000000000000000"), ROUND_HALF_UP)
    a = OMEGA * S[0] + (2 * B * S[0] + H) * Decimal((sqrt(1 - S[0]**2)))
    S[1] = ((-a * SQDJ * Decimal(sqrt(1 - S[0]**2)) + S[0] * Decimal(sqrt(1 + DJ**2 * Decimal(((1 - S[0]**2))) - Decimal(a**2)))
             ) / Decimal((1 + DJ**2 * Decimal(1 - S[0]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
    for i in range(1, N):
        A = OMEGA * S[i] + (2 * B * S[i] + H) * Decimal(sqrt(1 - S[i]**2)) + S[i - 1] * SQDJ * Decimal(sqrt(1 - S[i]**2)) - S[i] * \
            Decimal(sqrt(1 - S[i-1]**2))
        S[i + 1] = ((-A * SQDJ * Decimal(sqrt(1 - S[i]**2)) + S[i] * Decimal(sqrt(Decimal(1 - A**2) + DJ**2 * Decimal((1 - S[i]**2))))
                     ) / Decimal((1 + DJ**2 * Decimal(1 - S[i]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
        if isnan(S[i + 1]):
            break


"""Проверка на сходимость"""


def check(last: Decimal, prelast: Decimal) -> bool:
    temp = fabs(Decimal(1.0) - (last * Decimal(sqrt(1 - prelast**2)) -
                                (prelast + 2 * B * last) * Decimal(sqrt(1 - last**2))) / (last * OMEGA))
    return temp < EPS


"""Получение бризерного решения"""


def get_breather(amount: int, N: int):
    for i in prange(amount):
        counter(s0 + temp * i, N+1)
        if check(S[N], S[N - 1]):
            print(s0 + temp * i)
            counter(s0 + temp * i, N+1, H)
            writef(S)
            plot(S, N+1)


"""Расчет энергии"""


def energy():
    pass


def main():
    print(f'initial: {s0}, start max: {s0 + temp * DELIMITER}')
    print(f'Omega = {OMEGA}, eps = {EPS}\n')
    get_breather(DELIMITER, N)
    plot([-0.32, -0.28, -0.22, -0.2, -0.18], [0.165036570213893, 0.509457626013908,
                                              0.736884391556128, 0.788864920607633, 0.833085737950537])


if __name__ == "__main__":
    main()
