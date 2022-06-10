def lu(A, b, n, IDET):
    L = []
    contadorTroca = 0
    for i in range(n):
        L.append([0]*n)
        L[i][i] = 1

    for k in range(n-1):
        A, houveTroca = pivoteamento(A, k, n)
        if houveTroca: contadorTroca += 1
        for i in range(k + 1, n):
            m = - A[i][k] / A[k][k]
            L[i][k] = -m
            for j in range(k + 1, n): A[i][j] = m * A[k][j] + A[i][j]
            A[i][k] = 0

    det = 1
    for i in range(n): det = det *A[i][i]
    
    y = foward_elimination(L, b, n)
    x = backwards_elimination(A, y, n)
    
    if contadorTroca % 2 != 0: det = -det
    return (x, det)


def cholesky(A, b, n, IDET):
    g = [[0] * n for i in range(n)]
    y = [[0] * n for i in range(n)]
    gt = [[0] * n for i in range(n)]

    for i in range(n):
        if A[i][i] <= 0:
            return ("Error")
        for j in range(n):
            if A[i][j] != A[j][i]:
                return ("Error")

    for i in range(n):
        for k in range(i + 1):
            sum = 0
            for j in range(k): sum += g[i][j] * g[k][j]
            
            if (i == k): g[i][k] = (A[i][i] - sum) ** 0.5
            
            else: g[i][k] = (1.0 / g[k][k] * (A[i][k] - sum))
            gt[k][i] = g[i][k]

    y = foward_elimination(g, b, n)
    x = backwards_elimination(gt, y, n)
    det = 1
    for i in range(n): det = det * g[i][i]

    return (x, det)


def jacobi(A, b, n, IDET, tol = 0.0001, iterMax = 1000):
    error = tol + 1
    x = []
    history = []

    for i in range(n):
        for j in range(n):
            if i == j: continue
            A[i][j] = A[i][j] / A[i][i]
        b[i] = b[i] / A[i][i]
        x.append(b[i])

    for k in range(iterMax):
        previousX = x[:]
        for i in range(n):
            sum = 0
            for j in range(n):
                if i == j: continue
                sum += A[i][j] * previousX [j]
            x[i] = b[i] - sum
        differenceVector = [abs(previousXi - Xi) for (previousXi, Xi) in zip(previousX, x)]
        error = max(differenceVector) / abs(max(x, key=abs))
        history.append(error)
        if error <= tol: break

    return (x, k, history)
  

def gauss_seidel(A, b, n, IDET, tol = 0.0001, iterMax = 50):
    error = tol + 1
    history = []
    x = [0] * n
    x[0] = b[0]/A[0][0]

    for i in range(n):
        for j in range(n):
            if i == j: continue
            A[i][j] = A[i][j] / A[i][i]
        b[i] = b[i] / A[i][i]

    for k in range(iterMax):
        previousX = x[:]
        for i in range(n):
            sum = 0
            for j in range(n):
                if i == j: continue
                sum += A[i][j] * x[j]
            x[i] = b[i] - sum
        differenceVector = [abs(previousXi - Xi) for (previousXi, Xi) in zip(previousX, x)]
        error = max(differenceVector) / abs(max(x, key=abs))
        history.append(error)
        if error <= tol: break

    return (x, k, history)
    
def backwards_elimination(A, b, n):
    x = [0] * n
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i + 1, n): sum += A[i][j] * x[j]
        x[i] = (b[i] - sum) / A[i][i]
    return x

def foward_elimination(A, b, n):
    x = [0] * n
    for i in range(n):
        sum = 0
        for j in range (0, i): sum += A[i][j] * x[j]
        x[i] = (b[i] - sum) / A[i][i]
    return x

def pivoteamento(A, k, n):
    houveTroca = False
    pivo = A[k][k]
    posPivot = k
    for i in range(k+1, n):
        if abs(A[i][k]) > abs(pivo):
            pivo = A[i][k]
            posPivot = i
            houveTroca = True

    if houveTroca:
        A[k], A[posPivot] = A[posPivot], A[k]

    return A, houveTroca


path = ""
with open(path + 'input.csv', 'r') as entrada:
    lineCounter = 0
    matrizA, vectorB = [], []
    for line in entrada:
        if lineCounter == 0:
            n = int(line)
            lineCounter += 1
            continue
        elif lineCounter == 1:
            ICOD = int(line)
            lineCounter += 1
            continue
        elif lineCounter == 2:
            IDET = float(line)
            lineCounter += 1
            continue
        elif lineCounter == 3:
            Tolm = float(line)
            lineCounter += 1
            continue
        elif lineCounter > 2 and lineCounter <= n + 3:
            matrizA.append([int(x) for x in line.split(',')])
            lineCounter += 1
            continue
        else:
            vectorB = [float(x) for x in line.split(',')]
            lineCounter += 1
            continue
entrada.close()


print("Arquivos lidos, iniciando algoritmo selecionado")
if ICOD == 1:
    out = lu(matrizA, vectorB, n, IDET)
    print("Vetor resultado: " + str(out[0]) + "\n")
    output = open(path + "output.txt", "a")
    output.write("Lu:\nVetor resultado: " + str(out[0]) + "\n")
    if (IDET):
        output.write("Determinante: " + str(out[1]) + "\n\n")
        print("Determinante: " + str(out[1]) + "\n\n")
    output.close()

elif ICOD == 2:
    out = cholesky(matrizA, vectorB, n, IDET)
    if out[0] == "E":
        print("Erro: A matriz deve ser simétrica positiva\n\n")
        output = open("output.txt", "a")
        output.write("Erro: A matriz deve ser simétrica positiva.\n\n")
    else:    
        print("Vetor resultado: " + str(out[0]) + "\n")
        output = open("output.txt", "a")
        output.write("Cholesky:\nVetor resultado: " + str(out[0]) + "\n")
        if (IDET): 
            output.write("Determinante: " + str(out[1]) + "\n")
            print("Determinante: " + str(out[1]) + "\n")
    output.close()


elif ICOD == 3: 
    out = jacobi(matrizA, vectorB, n, IDET, Tolm)
    print("Vetor resultado: " + str(out[0]) + "\n")
    print("Número de iterações: " + str(out[1]) + "\n")
    output = open("output.txt", "a")
    output.write("Jacobi:\nVetor resultado: " + str(out[0]) + "\n")
    if (IDET): 
        output.write("Esse método não possui nenhuma facilitação para o cálculo da determinante.\n")
        print("Esse método não possui nenhuma facilitação para o cálculo da determinante.\n")
    output.write("Número de iterações: " +  str(out[1]) + "\n")
    output.write("Histórico do erro:\n" + str(out[2]) + "\n\n")
    output.close()

elif ICOD == 4:
    out = gauss_seidel(matrizA, vectorB, n, IDET)
    print("Gauss-seidel:\nVetor resultado: " + str(out[0]) + "\n")
    print("Número de iterações: " + str(out[1]) + "\n")
    output = open("output.txt", "a")
    output.write("Vetor resultado: " + str(out[0]) + "\n")
    if (IDET): 
        output.write("Esse método não possui nenhuma facilitação para o cálculo da determinante.\n")
        print("Esse método não possui nenhuma facilitação para o cálculo da determinante.\n")
    output.write("Número de iterações: " +  str(out[1]) + "\n")
    output.write("Histórico do erro:\n" + str(out[2]) + "\n\n")
    output.close()