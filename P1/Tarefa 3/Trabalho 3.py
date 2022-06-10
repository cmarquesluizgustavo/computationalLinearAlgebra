def lagrange(xi, yi, x, n):
    y = 0
    for i in range(n+1):
        p = 1
        for j in range(n+1):
            if i == j: continue
            p = p * (x - xi[j]) / (xi[i] - xi[j])
    
        y = y + p * yi[i]
    
    return y


def regression(xi, yi, x, n):
    meanX, meanY, sumXY, sumXX = 0, 0, 0, 0

    for i in range(n):
        meanX += xi[i] / n
        meanY += yi[i] / n
        sumXY += xi[i] * yi[i]
        sumXX += xi[i] * xi[i]

    sumXY -= n * meanX * meanY
    sumXX -= n * meanX * meanX

    b1 = sumXY / sumXX
    b0 = meanY - b1 * meanX

    y = b0 + b1 * x

    return y


path = ""
with open(path + 'input.csv', 'r') as entrada:
    lineCounter = 0
    a = []
    for line in entrada:
        if lineCounter == 0:
            ICOD = int(line)
            lineCounter += 1
            continue
        elif lineCounter == 1:
            n = int(line)
            xi, yi = [0] * n, [0] * n
            lineCounter += 1
            continue
        elif lineCounter > 1 and lineCounter <= n + 1:
            i = lineCounter - 2
            x, y = line.split(",")
            xi[i], yi[i] = float(x), float(y)
            lineCounter += 1
            continue
        elif lineCounter == n + 2:
            x = float(line)
            continue
entrada.close()
print("Arquivos lidos, iniciando algoritmo selecionado")

if ICOD == 1:
    out = lagrange(xi, yi, x, n - 1)
    print("Método da interpolação:\nValor y estimado: " + str(out))
    output = open(path + "output.txt", "a")

    if (x > max(xi) or x < min(xi)):
        print("Aviso: O ideal é escolher um valor de x para estimar o y entre valores passados no xi")
        output.write("Aviso: O ideal é escolher um valor de x para estimar o y entre valores passados no xi\n")

    output.write("Método da interpolação:\nValor y estimado: " + str(out) + "\n\n")
    output.close()

elif ICOD == 2:
    out = regression(xi, yi, x, n)
    print("Método da regressão:\nValor y estimado: " + str(out))
    output = open(path + "output.txt", "a")

    if (x > max(xi) or x < min(xi)):
        print("Aviso: O ideal é escolher um valor de x para estimar o y entre valores passados no xi")
        output.write("Aviso: O ideal é escolher um valor de x para estimar o y entre valores passados no xi\n")
    
    output.write("Método da regressão:\nValor y estimado: " + str(out) + "\n\n")
    output.close()