def lu(A, b):
    n = len(A)
    L = []
    for i in range(n):
        L.append([0]*n)
        L[i][i] = 1 

    for k in range(n-1): 
        for i in range(k + 1, n):
            m = - A[i][k] / A[k][k]
            L[i][k] = -m
            for j in range(k + 1, n): A[i][j] = m * A[k][j] + A[i][j]
            A[i][k] = 0
    
    y = foward_elimination(L, b, n) 
    x = backwards_elimination(A, y, n) 
    
    return x

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

def jacobiana(x):
    y = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    y[0][0] = (2 * x[0])
    y[0][1] = (4 * x[1])
    y[0][2] = (12 * x[2])
	
    y[1][0] = (12 * x[1] * x[0]) + (36 * x[1] * x[2])
    y[1][1] = (24 * x[1] ** 2) + (6 * x[0] ** 2) + (36 * x[0] * x[2]) + (108 * x[2] ** 2)
    y[1][2] = (36 * x[1] * x[0]) + (216 * x[1] * x[2])
    
    y[2][0] = (120 * x[1] ** 2 * x[0]) + (576 * x[1] ** 2 * x[2]) + (504 * x[2] ** 2 * x[0]) + (1296 * x[2] ** 3) + (72 * x[0] ** 2 * x[2]) + 3
    y[2][1] = (240 * x[1] ** 3) + (120 * x[1] * x[0] ** 2) + (1152 * x[1] * x[0] * x[2]) + (4464 * x[1] * x[2] ** 2)
    y[2][2] = (576 * x[1] ** 2 * x[0]) + (4464 * x[1] ** 2 * x[2]) + (504 * x[2] * x[0] ** 2) + (3888 * x[2] ** 2 * x[0]) + (13392 * x[2] ** 3) + (24 * x[0] ** 3)

    return y

def funcao(x, teta1, teta2):
    y = [0, 0, 0]

    y[0] = 2 * x[1] ** 2 + x[0] ** 2 + 6 * x[2] ** 2 - 1
    y[1] = 8 * x[1] ** 3 + 6 * x[1] * x[0] ** 2 + 36 * x[1] * x[0] * x[2] + 108 * x[1] * x[2] ** 2 - teta1
    y[2] = 60 * x[1] ** 4 + 60 * x[1] ** 2 * x[0] ** 2 + 576 * x[1] ** 2  * x[0] * x[2] + 2232 * x[1] ** 2 * x[2] ** 2 + 252 * x[2] ** 2 * x[0] ** 2 + 1296 * x[2] ** 3 * x[0] + 3348 * x[2] ** 4 + 24 * x[0] ** 3 * x[2] + 3 * x[0] - teta2

    return y

def product2vectorsDivided(a, b, c):
    a = [a[0]/c, a[1]/c, a[2]/c]
    result = [[], [], []]
    result[0] = [a[0] * b[0], a[0] * b[1], a[0] * b[2]]
    result[1] = [a[1] * b[0], a[1] * b[1], a[1] * b[2]]
    result[2] = [a[2] * b[0], a[2] * b[1], a[2] * b[2]]

    return result

def norma(vector):
    result = 0
    for i in range(len(vector)): result += vector[i] ** 2
    return result ** 0.5

def vectorDifference(a, b):
    n = len(a)
    vector = []
    for i in range(n): vector.append(a[i] - b[i])
    return vector

def matrix_addition(A, B):
    result = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            result[i][j] = A[i][j] + B[i][j]
    return result

def matrix_vector_multiplication(M, v):
    n = len(M)
    result = [0] * n
    for i in range(n):
        for j in range(n):
            result[i] += M[i][j] * v[j]
    return result 

def newton(teta1, teta2, Tolm = 0.000000001):  
    x = [1, 0, 0]
    maxIter = 100
    for i in range(maxIter):
        A = jacobiana(x)
        fx = funcao(x, teta1, teta2)
        delta = lu(A, fx)
        x = vectorDifference(x, delta)
        if norma(delta) < Tolm : break 

    return x  

def broyden (teta1, teta2, Tolm = 1):
    x = [1, 0, 0]
    A = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    fxAnterior = 0
    maxIter = 100
    for i in range(maxIter):
        fx = funcao(x, teta1, teta2)
        delta = lu(A, fx)
        previousFx = fx
        x = vectorDifference(x, delta)
        fx = funcao(x, teta1, teta2)
        y = [fx[0] - previousFx[0], fx[1] - previousFx[1], fx[2] - previousFx[2]]
        if norma(delta) / norma(x) < Tolm : break
        Atimesdelta = matrix_vector_multiplication(A, delta)
        yminus = vectorDifference(y, Atimesdelta)
        normaDelta = norma(delta) ** 2
        product = product2vectorsDivided(yminus, delta, normaDelta)
        A = matrix_addition(A, product)

    return x



#Analysing input.csv as requested
path = ""
with open("input.csv", "r") as entrada:
    lineCounter = 0
    matrizA, vectorB = [], []
    for line in entrada:
        if lineCounter == 0:
            ICOD = int(line)
            lineCounter += 1
            continue
        elif lineCounter == 1:
            teta1, teta2 = line.split(",")
            teta1, teta2 = float(teta1), float(teta2)
            lineCounter += 1
            continue
        else:
            Tolm = float(line)
            lineCounter += 1
            break
entrada.close()




output = open(path + "output.txt", "a")

if ICOD == 1:
    result = newton(teta1, teta2, Tolm)
    output.write("Resultado do método de Newton para teta1 = " + str(teta1) + " e teta2 = " + str(teta2) + " é: " + "\n" + str(result) + "\n\n\n")
    print("Resultado do método de Newton para teta1 = " + str(teta1) + " e teta2 = " + str(teta2) + " é: " + "\n" + str(result) + "\n\n\n")
    output.close()

if ICOD == 2:
    result = broyden(teta1, teta2, Tolm)
    output.write("Resultado do método de Broyden para teta1 = " + str(teta1) + " e teta2 = " + str(teta2) + " é: " + "\n" + str(result) + "\n\n\n")
    print("Resultado do método de Broyden para teta1 = " + str(teta1) + " e teta2 = " + str(teta2) + " é: " + "\n" + str(result) + "\n\n\n")
    output.close()