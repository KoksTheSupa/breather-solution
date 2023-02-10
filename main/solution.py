import _thread

from decimal import *
from math import fabs, sqrt, isnan
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing as mp

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d

N = 100
DJ = Decimal("0.16")
SQDJ = Decimal(sqrt(1 + DJ**2))
OMEGA = Decimal("-0.28")
B = Decimal("0.15")
BETA = Decimal("1.5")
EPS = Decimal("0.000001")
S = [Decimal(0.0)] * (N+1)
DELIMITER = 2_000_000_000  # double results[2000000000]
s0 = Decimal("0.2")
temp = Decimal("0.4") / DELIMITER


"""Построение графиков"""


def plot(x, y):
    x = [i for i in range(x)]
    xmax = x[y.index(max(y))]
    fig, graph = plt.subplots()
    graph.scatter(x, y, 30, 'black')
    graph.plot(x, y, linewidth=1.0)
    graph.spines['bottom'].set_position('center')  # type: ignore
    graph.set_ylabel('E')
    graph.set_xlabel('N')
    graph.grid(True)
    graph.annotate(f'Local max = {max(y)}', xy=(xmax, max(y)), xytext=(xmax, max(
        y) + Decimal(0.10)), arrowprops=dict(facecolor='black', shrink=0.05, headlength=5, width=2))
    graph.set_ylim(-max(y) - Decimal('0.15'), max(y) + Decimal('0.15'))
    plt.show()


"""Построение 3D модели"""


def plot3D(x, y):
    x = [i for i in range(x)]
    z = y*2
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    for xi, yi, zi in zip(z, x, y):
        if xi < 0:
            line = art3d.Line3D(*zip((xi, yi, 0), (xi, yi, zi)),
                                marker=11, markevery=(1, 1))
        else:
            line = art3d.Line3D(*zip((xi, yi, 0), (xi, yi, zi)),
                                marker=10, markevery=(1, 1))
        ax.add_line(line)
    ax.set(xlabel="Sx", ylabel="N", zlabel="Sy")
    ax.set_xlim3d(-0.8, 0.8)  # type: ignore
    ax.set_ylim3d(0, 101)  # type: ignore
    ax.set_zlim3d(-0.8, 0.8)  # type: ignore
    ax.grid(False)
    plt.show()


"""Запись в файл"""


def writef(chain: list):
    with open(f'e_N({N+1})_Omega{OMEGA}_B{B}_dj{DJ*100}.nb', 'a+') as file:
        file.write(
            f'text=\"parameters: omega={OMEGA}, N(real)=({N+1}), B={B}, dj={DJ}, S[N]={S[0]} (EDGE algorithm)\"\ndataS = {{')
        for i in range(N):

            file.write('{%s,%s},' % (i, (S[i])))

        file.write('{%s,%s}};' % (N, S[N]))
        file.write(
            "\n\nListPlot[dataS,Joined->True,Mesh->All,PlotRange->All]\n\n")


"""Расчет цепочки"""


def counter(start: Decimal, N: int):
    S[0] = Decimal(str(start)).quantize(
        Decimal("1.000000000000000"), ROUND_HALF_UP)
    a = OMEGA * S[0] + (2 * B * S[0]) * Decimal((sqrt(1 - S[0]**2)))
    S[1] = ((-a * SQDJ * Decimal(sqrt(1 - S[0]**2)) + S[0] * Decimal(sqrt(1 + DJ**2 * Decimal(((1 - S[0]**2))) - Decimal(a**2)))
             ) / Decimal((1 + DJ**2 * Decimal(1 - S[0]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
    for i in range(1, N):
        A = OMEGA * S[i] + (2 * B * S[i]) * Decimal(sqrt(1 - S[i]**2)) + S[i - 1] * SQDJ * Decimal(sqrt(1 - S[i]**2)) - S[i] * \
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
    for i in range(amount):
        counter(s0 + temp * i, N)
        if check(S[N], S[N - 1]):
            print(s0 + temp * i)
            counter(s0 + temp * i, N)
            writef(S)
            plot(N+1, S)


"""Расчет энергии"""


def energy(chain: list):
    E = 0
    E0 = Decimal("-2.35")
    sum1, sum2, sum3, sum4 = 0, 0, 0, 0
    for i in range(N):
        sum1 = Decimal(chain[i]*chain[i+1])
        sum2 = (1-chain[i]*chain[i])
        sum3 = Decimal(sqrt((1 - chain[i]*chain[i])
                            * (1 - chain[i+1]*chain[i+1])))
        sum4 = Decimal(sqrt(1 - chain[i]*chain[i]))
        E += (-Decimal(sqrt(1 + DJ**2))*sum1 + Decimal("0.15")
              * sum2 - sum3 - Decimal("1.5")*sum4 - E0)/N


def main():
    print(f'initial: {s0}, start max: {s0 + temp * DELIMITER}')
    print(f'Omega = {OMEGA}, eps = {EPS}\n')
    p = mp.Process(get_breather(DELIMITER, N))
    mp.Pool(processes=4)
    p.start()
    energy(S)
    plot3D(N + 1, S)


if __name__ == "__main__":
    main()
