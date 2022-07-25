from math import sqrt, pow, isnan
import matplotlib.pyplot as plt
from numba import prange
from decimal import *
from multiprocessing import Pool

N = 100
DJ = Decimal("0.16")
OMEGA = Decimal("-0.2")
B = Decimal("0.15")
EPS = H = Decimal(str(10**-6))
S = [0] * (N + 1)#np.zeros(N + 1)
DELIMITER = 2000000000 # double results[2000000000]
s0 = Decimal("0.0")
temp = Decimal(Decimal("0.4") / (DELIMITER))


"""Построение графиков"""
def plot(x, y):
    fig, graph = plt.subplots()
    graph.scatter([i for i in prange(x)], y, 30, 'black' )
    graph.plot([i for i in prange(x)], y, linewidth=1.0)
    graph.spines['bottom'].set_position('center')
    graph.set_ylabel('S') 
    plt.show()


def counter(start: float,N: int) -> bool:
    S[0] = Decimal(str(start)).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
    a = OMEGA * S[0] + 2 * B * S[0] * Decimal((sqrt(1 - S[0]**2)))
    S[1] = ((-a * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[0]**2)) + S[0] * Decimal(sqrt(1 +DJ**2 * Decimal(((1 - S[0]**2))) - Decimal(a**2)))) / Decimal((1 + DJ**2 * Decimal(1 - S[0]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)  
    for i in range(1, N):
        A = OMEGA * S[i] + 2 * B * S[i] * Decimal(sqrt(1 - S[i]**2)) + S[i - 1] * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) - S[i] * Decimal(sqrt(1 - S[i-1]**2))
        S[i + 1] = ((-A * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) + S[i] * Decimal(sqrt(Decimal(1 - A**2) + DJ**2 * Decimal((1 - S[i]**2))))) / Decimal((1 + DJ**2 * Decimal(1 - S[i]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
        if isnan(S[i + 1]): break
    return(S)

def check(last, prelast):
    temp = abs(Decimal(1.0) - (last * Decimal(sqrt(1 - last**2)) - (prelast + 2 * B * S[N]) * Decimal(sqrt(1 - last**2))) / (last * OMEGA))
    if temp < EPS:
        return 1
    else: return 0

def get_breather(amount,N):
    for i in range(amount):
        counter(s0 + temp * i,N)
        if check(S[N], S[N - 1]):
            return(s0 + temp * i)


def main():
    print(f'initial: {(s0)}, start max: {s0 + temp * DELIMITER}')
    print(f'Omega = {OMEGA}, eps = {EPS}\n')
    
    
    S[0] = get_breather(DELIMITER,N)
    a = OMEGA * S[0] + Decimal((2 * B * S[0] + H)) * Decimal(sqrt(1 - S[0]**2))
    S[1] = ((-a * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[0]**2)) + S[0] * Decimal(sqrt(1 + DJ**2 * (1 - S[0]**2) - a**2))) / Decimal((1 + DJ**2 * (1 - S[0]**2)))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
   
    for i in range(1, N):
        A = OMEGA * S[i] + Decimal(2 * B * S[i] + H) * Decimal(sqrt(1 - S[i]**2)) + S[i - 1] * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) - S[i] * Decimal(sqrt(1 - S[i - 1]**2))
        S[i + 1] = ((-A * Decimal(sqrt(1 + DJ**2)) * Decimal(sqrt(1 - S[i]**2)) + S[i] * Decimal(sqrt(1 - A**2 + DJ**2 * (1 - S[i]**2)))) / (1 + DJ**2 * (1 - S[i]**2))).quantize(Decimal("1.000000000000000"), ROUND_HALF_UP)
        if isnan(S[i + 1]): break
    print(f'Left: {S[N] * OMEGA}, right: {S[N] * Decimal(sqrt(1 - S[N-1]**2)) - (S[N - 1] + 2 * B * S[N]) * Decimal(sqrt(1 - S[N]**2))}')
    print(S)
    if (abs(Decimal('1.0') - (S[N] * Decimal(sqrt(1 - S[N - 1]**2)) - (S[N - 1] + 2 * B * S[N]) * Decimal(sqrt(1 - S[N]**2))) / (S[N] * OMEGA)) < EPS):
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