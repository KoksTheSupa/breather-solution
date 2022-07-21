from math import sqrt, pow, isnan
# import matplotlib.pyplot as plt
# from numba import prange
from functools import lru_cache
from decimal import *

N = 100
DJ = 0.16
OMEGA = -0.22
B = 0.15
EPS = H = pow(10, -6)
S = [0] * (N + 1)
DELIMITER = 2000000000 # double results[2000000000]
s0 = 0.72850
temp = (0.4) / (DELIMITER)



"""Функция, позволяющая задать нужное число знаков после запятой"""
def toFixed(numObj: float, digits: int) -> str:
    return f'{numObj:.{digits}f}'

@lru_cache(maxsize= None)
def counter(start: float) -> bool:

    S[0] = start
    a = OMEGA * S[0] + 2 * B * S[0] * sqrt(1 - pow(S[0], 2))
    S[1] = (-a * sqrt(1 + DJ**2) * sqrt(1 - pow(S[0], 2)) + S[0] * sqrt(1 + DJ**2 * (1 - pow(S[0], 2)) - pow(a, 2))) / (1 + DJ**2 * (1 - pow(S[0], 2)))  

    for i in range(1, N):
        A = OMEGA * S[i] + 2 * B * S[i] * sqrt(1 - pow(S[i], 2)) + S[i - 1] * sqrt(1 + DJ**2) * sqrt(1 - pow(S[i], 2)) - S[i] * sqrt(1 - pow(S[i - 1], 2))
        S[i + 1] = (-A * sqrt(1 + DJ**2) * sqrt(1 - pow(S[i], 2)) + S[i] * sqrt(1 - pow(A, 2) + DJ**2 * (1 - pow(S[i], 2)))) / (1 + DJ**2 * (1 - pow(S[i], 2)))

        if isnan(S[i + 1]): break

    if abs(abs(Decimal(S[0])) - abs(Decimal(S[N]))) < 0.01: return 0
    temp = abs(1.0 - (S[N] * sqrt(1 - pow(S[N - 1], 2)) - (S[N - 1] + 2 * B * S[N]) * sqrt(1 - pow(S[N], 2))) / (S[N] * OMEGA))
    if Decimal(temp) < Decimal(EPS):
        return 1
    else: return 0

def selection(amount: int) -> float:
    for i in range(amount):
        examp = toFixed(Decimal(s0 + temp * i),15)
        if counter(s0 + temp * i) > 0:
            print(s0 + temp * i)
            return(s0 + temp * i)


def main():
    print(f'initial: {(s0)}, start max: {s0 + temp * DELIMITER}')
    print(f'Omega = {OMEGA}, eps = {EPS}\n')
    S[0] = selection(DELIMITER)
    a = OMEGA * S[0] + (2 * B * S[0] + H) * sqrt(1 - pow(S[0], 2))
    S[1] = (-a * sqrt(1 + DJ**2) * sqrt(1 - pow(S[0], 2)) + S[0] * sqrt(1 + DJ**2 * (1 - pow(S[0], 2)) - pow(a, 2))) / (1 + DJ**2 * (1 - pow(S[0], 2)))
   
    for i in range(1, N):
        A = OMEGA * S[i] + (2 * B * S[i] + H) * sqrt(1 - pow(S[i], 2)) + S[i - 1] * sqrt(1 + DJ**2) * sqrt(1 - pow(S[i], 2)) - S[i] * sqrt(1 - pow(S[i - 1], 2))
        S[i + 1] = (-A * sqrt(1 + DJ**2) * sqrt(1 - pow(S[i], 2)) + S[i] * sqrt(1 - pow(A, 2) + DJ**2 * (1 - pow(S[i], 2)))) / (1 + DJ**2 * (1 - pow(S[i], 2)))
        
        if isnan(S[i + 1]): break

    print(f'Left: {S[N] * OMEGA}, right: {S[N] * sqrt(1 - pow(S[N - 1], 2)) - (S[N - 1] + 2 * B * S[N]) * sqrt(1 - pow(S[N], 2))}')
    
    if (abs(1.0 - (S[N] * sqrt(1 - pow(S[N - 1], 2)) - (S[N - 1] + 2 * B * S[N]) * sqrt(1 - pow(S[N], 2))) / (S[N] * OMEGA)) < 0.000001 * 1):
        """Запись данных в файл"""
        with open(f'e_N({N+1})_Omega{OMEGA}_B{B}_dj{DJ*100}_h{H}.nb','a+') as file:
            file.write(f'text=\"parameters: omega={OMEGA}, N(real)=({N+1}), B={B}, dj={DJ}, h={H}, S[N]={toFixed(S[0], 15)} (EDGE algorithm)\"\ndataS = {{')
            for i in range(N):

                    file.write('{%s,%s},'% (i, (toFixed(S[i], 15))))
            
            file.write('{%s,%s}};'% (N, (toFixed(S[N], 15))))
            file.write("\n\nListPlot[dataS,Joined->True,Mesh->All,PlotRange->All]\n\n")

    # fig, graph = plt.subplots()
    # graph.scatter([x for x in prange(N + 1)], S, 30, 'black' )
    # graph.plot([x for x in prange(N + 1)], S, linewidth=1.0)
    # graph.spines['bottom'].set_position('center')
    # graph.set_ylabel('S') 
    # plt.show()


if __name__ == "__main__":
    main()