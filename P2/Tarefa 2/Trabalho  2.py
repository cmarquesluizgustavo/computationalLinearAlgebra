import math

def função (C, x):
    y = C[0] * math.exp(C[1] * x) + C[2] * x ** C[3]
    return y

def função_derivada(C, x):
    y = C[0] * C[1] * math.exp(C[1] * x) + C[2] * C[3] * x ** (C[3] - 1)
    return y

def bissection(C, a, b, Tolm):
    if função(C, a) * função(C, b) > 0:
        return "Não há raízes ou há múltiplas raízes no intervalo ."

    while abs(b - a) > Tolm:
        x = (a + b) / 2
        if função(C, x) == 0:
            return x
        if função(C, a) * função(C, x) < 0:
            b = x
        else:
            a = x
    return x

def newton_root (C, a, b, Tolm = 0.0001):
    maxIter = 100
    if função(C, a) * função(C, b) > 0:
        return "Não há raízes ou há múltiplas raízes no intervalo."

    x = (a + b) / 2
    for i in range(maxIter):
        xOld = x
        x = x - função(C, x) / função_derivada(C, x)
        if abs(xOld - x) < Tolm:
            return x
    return "Último resultado: " + str(x) + ", mas não houve convergência em " + str(maxIter) + " iterações."

 
def gaussian_quadrature(C, a, b, n):
    dt = (b - a) /  2
    aux = (a + b) / 2
    I = 0
    if n == 0 or n == 1: return "Erro, não há como integrar sem pontos, aumente o n!."
    if n > 10: return "Não há dados para mais que 10 pontos, diminua o n!"
    table = [ [], [[2], [0]],
            [ [1.0, 1.0], [0.57735, -0.57735] ],
            [ [0.555556, 0.88889, 0.55556], [0.77459, 0, -0.77459] ],
            [ [0.34786, 0.65215, 0.65215, 0.34786], [0.86114, 0.33998, -0.33999, -0.86114] ],
            [ [0.23693, 0.47863, 0.56889, 0.47863, 0.23693], [0.90618, 0.53847, 0, -0.53847, -0.90618] ],
            [ [0.17133, 0.36076, 0.46791, 0.46791, 0.36076, 0.17133], [0.93247, 0.66121, 0.23862, -0.23862, -0.66121, -0.93247] ],
            [ [0.12949, 0.27971, 0.38183, 0.41796, 0.38183, 0.27971, 0.12949], [0.94911, 0.74153, 0.40585, 0, -0.40585, -0.74153, -0.94911] ],
            [ [0.10123, 0.22238, 0.31371, 0.36268, 0.36268, 0.31371, 0.22238, 0.10123], [0.96029, 0.79667, 0.53553, 0.18343, -0.18343, -0.53553, -0.79667, -0.96029] ],
            [ [0.08127, 0.18065, 0.26061, 0.31235, 0.33024, 0.31235, 0.26061, 0.18065, 0.08127], [0.96816, 0.83603, 0.61337, 0.32425, 0, -0.32425, -0.61337, -0.83603, -0.96816] ],
            [ [0.06667, 0.14945, 0.21909, 0.26927, 0.29552, 0.29552, 0.26927, 0.21909, 0.14945, 0.06667], [0.97391, 0.86506, 0.67941, 0.43339, 0.14887, -0.14887, -0.43339, -0.67941, -0.86506, -0.97391] ]
            ]
    for i in range(n):
        xi = dt * table[n][1][i] + aux
        I += função(C, xi) * table[n][0][i]
    return I * dt

def polinomial_quadrature(C, a, b, n):
    L = b - a
    delta = (b - a) / (n - 1)
    x = [a + delta * i for i in range(n)]
    I = 0
    if n == 0 or n == 1: return "Erro, não há como integrar sem pontos, aumente o n!."
    if n == 2: w = [L/2,L/2]
    elif n == 3: w = [L/6,(2*L)/3,L/6]
    elif n == 4: w = [L/8,(3*L)/8, (3*L)/8, L/8]
    elif n == 5: w = [(7*L)/90, (16*L)/45, (2*L)/15, (16*L)/45, (7*L)/90]
    if n > 10: return "Não há dados para mais que 5 pontos, diminua o n!"


    for i in range(n):
        I += função(C, x[i]) * w[i]
    return I

def derivative_foward (C, x, delta):
    fx = função(C, x)
    fdx = função(C, x + delta)
    return (fdx - fx) / delta

def derivative_backwards (C, x, delta):
    f = função(C, x)
    fdx = função(C, x - delta)
    return (f - fdx) / delta

def derivative_central (C, x, delta):
    fplus = função(C, x + delta)
    fminus = função(C, x - delta)
    return (fplus - fminus) / (2 * delta)

def richard_extrapolation(C, x, delta, delta2):
    d1 = derivative_foward(C, x, delta)
    d2 = derivative_foward(C, x, delta2)
    q = delta2 / delta
    return d1 + (d1 - d2) / (q - 1)


#Analysing input.csv as requested
path = ""
with open(path + 'input.csv', 'r') as entrada:
    lineCounter = 0
    matrizA, vectorB = [], []
    for line in entrada:
        if lineCounter == 0:
            if line == "4\n" or line == "4.0\n":
                ICOD1 = 4
            else:
                ICOD1, ICOD2 = line.split(",")
                ICOD1, ICOD2 = int(ICOD1), int(ICOD2)
            lineCounter += 1
            continue
        elif lineCounter == 1:
            c1, c2, c3, c4 = line.split(",")
            c1, c2, c3, c4 = float(c1), float(c2), float(c3), float(c4)
            lineCounter += 1
            continue
        elif lineCounter == 2:
            a, b = line.split(",")
            a, b = float(a), float(b)
            if ICOD1 > 2:
                x = a
                delta = b
                if ICOD1 == 3:
                    break
            lineCounter += 1
            continue
        else:
            Tolm = float(line)
            lineCounter += 1
            if ICOD1 == 2:
                n = int(Tolm)
            elif ICOD1 == 4:
                delta2 = Tolm
            break
entrada.close()

C = [c1, c2, c3, c4]

output = open(path + "output.txt", "a")

if ICOD1 == 1:
    output.write("Constantes c1, c2, c3, c4: " + str([c1, c2, c3, c4]) + "\n" + "As margens de análise são: " + str(a) + " e " + str(b) + "\n")
    if ICOD2 == 1:
        result = bissection(C, a, b, Tolm)
        output.write("O resultado da Bisseção é " + str(result) + "\n")
    else:
        result = newton_root(C, a, b, Tolm)
        output.write("O resultado da Newton é " + str(result) + "\n")

elif ICOD1 == 2:
    output.write("Constantes c1, c2, c3, c4: " + str([c1, c2, c3, c4]) + "\n" + "As margens de análise são: " + str(a) + " e " + str(b) + ". E o número de pontos é " + str(n) + "\n")
    if ICOD2 == 1:
        result = gaussian_quadrature(C, a, b, n)
        output.write("O resultado da quadratura gaussiana é " + str(result) + "\n")
    else:
        result = polinomial_quadrature(C, a, b, n)
        output.write("O resultado da quadratura polinomial é " + str(result) + "\n")

elif ICOD1 == 3:
    output.write("Constantes c1, c2, c3, c4: " + str([c1, c2, c3, c4]) + "\n" + "O x é " + str(x) + " e o delta é " + str(delta) + "\n")
    if ICOD2 == 1:
        result = derivative_foward(C, x, delta)
        output.write("O resultado da derivada para frente é " + str(result) + "\n")
    elif ICOD2 == 2:
        result = derivative_foward(C, x, delta)
        output.write("O resultado da derivada para trás é " + str(result) + "\n")
    else:
        result = derivative_central(C, x, delta)
        output.write("O resultado da derivada central é " + str(result) + "\n")
else:
    output.write("Constantes c1, c2, c3, c4: " + str([c1, c2, c3, c4]) + "\n" + "O x é: " + str(x) + " e o delta é" + str(delta) + " o delta2 é " + str(delta2) + "\n")
    result = richard_extrapolation(C, x, delta, delta2)
    output.write("O resultado da extrapolação de Richard é " + str(result) + "\n")


output.write("\n\n")



output.close()