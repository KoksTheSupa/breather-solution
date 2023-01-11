from decimal import *
from math import sqrt, isnan
from numba import prange
from random import randint

import matplotlib.pyplot as plt

N = 100
DJ = Decimal("0.16")
OMEGA = Decimal("-0.28")
B = Decimal("0.15")
EPS = H = Decimal(str(10**-6))
S = [0] * (N+1)  # np.zeros(N + 1)
DELIMITER = 2000000000  # double results[2000000000]
s0 = Decimal("0.495373257505000")
temp = Decimal(Decimal("0.4") / (DELIMITER))


"""Построение графиков"""


def plot(x, y):
    x = [i for i in prange(x)]
    xmax = x[y.index(max(y))]
    fig, graph = plt.subplots()
    graph.scatter(x, y, 30, 'black')
    graph.plot(x, y, linewidth=1.0)
    graph.spines['bottom'].set_position('center')
    graph.set_ylabel('S')
    graph.set_xlabel('N')
    graph.annotate(f'Local max = {max(y)}', xy=(xmax, max(y)), xytext=(xmax, max(y) + Decimal(0.10)),
                   arrowprops=dict(facecolor='black',
                                   shrink=0.05, headlength=5, width=2),
                   )
    graph.set_ylim(-max(y) - Decimal('0.15'), max(y) + Decimal('0.15'))
    graph.grid(True)
    plt.show()


"""Запись в файл"""


def writef(S):
    with open(f'e_N({N+1})_Omega{OMEGA}_B{B}_dj{DJ*100}_h{H}.nb', 'a+') as file:
        file.write(
            f'text=\"parameters: omega={OMEGA}, N(real)=({N+1}), B={B}, dj={DJ}, h={H}, S[N]={S[0]} (EDGE algorithm)\"\ndataS = {{')
        for i in range(N):

            file.write('{%s,%s},' % (i, (S[i])))

        file.write('{%s,%s}};' % (N, S[N]))
        file.write(
            "\n\nListPlot[dataS,Joined->True,Mesh->All,PlotRange->All]\n\n")


"""Расчет цепочки"""


def counter(start: float, N: int):
    S[0] = Decimal(str(start)).quantize(
        Decimal("1.000000000000000"), ROUND_HALF_UP)
    a = OMEGA * S[0] + 2 * B * S[0] * Decimal((sqrt(1 - S[0]**2)))
    S[1] = ((-a * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[0]**2)) + S[0] * Decimal(sqrt(1 + DJ**2 * Decimal(((1 - S[0]**2))
                                                                                                                  ) - Decimal(a**2)))) / Decimal((1 + DJ**2 * Decimal(1 - S[0]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
    for i in range(1, N-1):
        A = OMEGA * S[i] + 2 * B * S[i] * Decimal(sqrt(1 - S[i]**2)) + S[i - 1] * Decimal(
            sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) - S[i] * Decimal(sqrt(1 - S[i-1]**2))
        S[i + 1] = ((-A * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) + S[i] * Decimal(sqrt(Decimal(1 - A**2) + DJ**2 *
                                                                                                      Decimal((1 - S[i]**2))))) / Decimal((1 + DJ**2 * Decimal(1 - S[i]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
        if isnan(S[i + 1]):
            break
    return (S)


def check(last, prelast):
    temp = abs(Decimal(1.0) - (last * Decimal(sqrt(1 - last**2)) -
               (prelast + 2 * B * S[N]) * Decimal(sqrt(1 - last**2))) / (last * OMEGA))
    if temp < EPS:
        return 1
    else:
        return 0


def get_breather(amount, N):
    for i in prange(amount):
        counter(s0 + temp * i, N)
        if check(S[N], S[N - 1]):
            return (s0 + temp * i)


def main():
    print(f'initial: {(s0)}, start max: {s0 + temp * DELIMITER}')
    print(f'Omega = {OMEGA}, eps = {EPS}\n')
    plot(N + 1, counter(s0, N))
    writef(S)
    plt.plot([-0.32, -0.28, -0.22, -0.2, -0.18], [0.165036570213893, 0.509457626013908,
                                                  0.736884391556128, 0.788864920607633, 0.833085737950537])
    plt.show()


if __name__ == "__main__":
    main()
