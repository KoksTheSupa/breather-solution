from math import sqrt, pow, isnan
import matplotlib.pyplot as plt
from numba import prange
from decimal import *

N = 100
DJ = Decimal("0.16")
OMEGA = Decimal("-0.2")
B = Decimal("0.15")
EPS = H = Decimal(str(10**-6))
S = [0] * (N + 1)#np.zeros(N + 1)
DELIMITER = 2000000000 # double results[2000000000]
s0 = 0.8
temp = (0.4) / (DELIMITER)


"""Построение графиков"""
def plot(x, y):
    fig, graph = plt.subplots()
    graph.scatter([i for i in prange(x)], y, 30, 'black' )
    graph.plot([i for i in prange(x)], y, linewidth=1.0)
    graph.spines['bottom'].set_position('center')
    graph.set_ylabel('S') 
    plt.show()


def counter(start: float) -> bool:
    S[0] = Decimal(str(start)).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
    a = OMEGA * S[0] + 2 * B * S[0] * Decimal((sqrt(1 - S[0]**2)))
    S[1] = ((-a * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[0]**2)) + S[0] * Decimal(sqrt(1 +DJ**2 * Decimal(((1 - S[0]**2))) - Decimal(a**2)))) / Decimal((1 + DJ**2 * Decimal(1 - S[0]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)  
    for i in range(1, N):
        A = OMEGA * S[i] + 2 * B * S[i] * Decimal(sqrt(1 - S[i]**2)) + S[i - 1] * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) - S[i] * Decimal(sqrt(1 - S[i-1]**2))
        S[i + 1] = ((-A * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) + S[i] * Decimal(sqrt(Decimal(1 - A**2) + DJ**2 * Decimal((1 - S[i]**2))))) / Decimal((1 + DJ**2 * Decimal(1 - S[i]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
        if isnan(S[i + 1]): break
    return(S)
    if abs(abs(S[0]) - abs(S[N])) < 0.01: return 0
    return(S)
    temp = abs(1.0 - (S[N] * sqrt(1 - pow(S[N - 1], 2)) - (S[N - 1] + 2 * B * S[N]) * sqrt(1 - pow(S[N], 2))) / (S[N] * OMEGA))
    if temp < EPS:
        return 1
    else: return 0

def selection(amount: int) -> float:
    for i in prange(amount):
        if counter(s0 + temp * i) > 0:
            print(s0 + temp * i)
            return(s0 + temp * i)

def main():
    print(f'initial: {(s0)}, start max: {s0 + temp * DELIMITER}')
    print(f'Omega = {OMEGA}, eps = {EPS}\n')
    print(counter(0.556926333500000))
    S[0] = 0.556926333500000
    a = OMEGA * S[0] + (2 * B * S[0] + H) * sqrt(1 - pow(S[0], 2))
    S[1] = (-a * sqrt(1 + DJ**2) * sqrt(1 - pow(S[0], 2)) + S[0] * sqrt(1 + DJ**2 * (1 - pow(S[0], 2)) - pow(a, 2))) / (1 + DJ**2 * (1 - pow(S[0], 2)))
   
    for i in range(1, N):
        A = OMEGA * S[i] + (2 * B * S[i] + H) * sqrt(1 - pow(S[i], 2)) + S[i - 1] * sqrt(1 + DJ**2) * sqrt(1 - pow(S[i], 2)) - S[i] * sqrt(1 - pow(S[i - 1], 2))
        S[i + 1] = (-A * sqrt(1 + DJ**2) * sqrt(1 - pow(S[i], 2)) + S[i] * sqrt(1 - pow(A, 2) + DJ**2 * (1 - pow(S[i], 2)))) / (1 + DJ**2 * (1 - pow(S[i], 2)))
        
        if isnan(S[i + 1]): break

    print(f'Left: {S[N] * OMEGA}, right: {S[N] * sqrt(1 - pow(S[N - 1], 2)) - (S[N - 1] + 2 * B * S[N]) * sqrt(1 - pow(S[N], 2))}')
    
    if (abs(1.0 - (S[N] * sqrt(1 - pow(S[N - 1], 2)) - (S[N - 1] + 2 * B * S[N]) * sqrt(1 - pow(S[N], 2))) / (S[N] * OMEGA)) < 0.000001):
        """Запись данных в файл"""
        with open(f'e_N({N+1})_Omega{OMEGA}_B{B}_dj{DJ*100}_h{H}.nb','a+') as file:
            file.write(f'text=\"parameters: omega={OMEGA}, N(real)=({N+1}), B={B}, dj={DJ}, h={H}, S[N]={S[0]} (EDGE algorithm)\"\ndataS = {{')
            for i in range(N):

                    file.write('{%s,%s},'% (i, (S[i])))
            
            file.write('{%s,%s}};'% (N, S[N]))
            file.write("\n\nListPlot[dataS,Joined->True,Mesh->All,PlotRange->All]\n\n")

    plot(N+1, S)

if __name__ == "__main__":
    main()