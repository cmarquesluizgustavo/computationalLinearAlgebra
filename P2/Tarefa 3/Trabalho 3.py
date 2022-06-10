import math

def função(t, a, w):
    return a[0] * math.sin(w[0] * t) + a[1] * math.sin(w[1] * t) + a[2] * math.sin(w[2] * t) 

def ed(t, y, yderivado, c, k, m, a, w):
    return (função(t, a, w) - c * yderivado - k * y) / m

def RKN (h = 0.1, T = 1, a = [1, 2, 1.5], w = [0.05, 1, 2], m = 1, c = 0.1, k = 2, x = 0, v = 0):
    n = T / h
    tprevious = 0
    results = [[], [], [], []]
    for i in range (1, int(n) + 1):
        t = i * h
        k1 = (h / 2) * ed(tprevious, x, v, c, k, m, a, w)
        K = (h / 2) * (v + k1 / 2)
        k2 = (h / 2) * ed(tprevious + (h / 2), x + K, v + k1, c, k, m, a, w)
        k3 = (h / 2) * ed(tprevious + (h / 2), x + K, v + k2, c, k, m, a, w)
        L = h * (v + k3)
        k4 = (h / 2) * ed(tprevious + h, x + L, v + 2 * k3, c, k, m, a, w)
        x = x + h * (v + (1 / 3) * (k1 + k2 + k3))
        v = v + (1 / 3)  * (k1 + 2 * k2 + 2 * k3 + k4)
        tprevious = t
        aceleration = ed(t, x, v, c, k, m, a, w)
        results[0].append(t)
        results[1].append(x)
        results[2].append(v)
        results[3].append(aceleration)


    return results

#Analysing input.csv as requested
path = ""
flag = False
with open(path + 'input.csv', 'r') as entrada:
    lineCounter = 0
    matrizA, vectorB = [], []
    for line in entrada:
        if lineCounter == 0:
            h = float(line)
            lineCounter += 1
            continue
        elif lineCounter == 1:
            T = float(line)
            lineCounter += 1
            continue
        if line:
            m, c, k, a1, a2, a3, w1, w2, w3 = line.split(",")
            m, c, k, a, w = float(m), float(c), float(k), [float(a1), float(a2), float(a3)], [float(w1), float(w2), float(w3)]
            flag = True
        break
entrada.close()

if not flag:
    print("hey")
    a = [1, 2, 1.5]
    w = [0.05, 1, 2]
    m = 1
    c = 0.1
    k = 2


print("Resultado do Runge-Kutta-Nystron para as entradas:")
print("m =", m, "c =", c, "k =", k)
print("a1 =", a[0], "a2 =", a[1], "a3 =", a[2])
print("w1 =", w[0], "w2 =", w[1], "w3 =", w[2])
print("Tempo = ", T, "Passo = ", h)

results = RKN(h, T, a, w, m, c, k)
output = open(path + "output.txt", "a")
output.write("Resultado do Runge-Kutta-Nystron para as entradas: \nm = " + str(m) + " c = " + str(c) + " k = " + str(k) + "\nw1 = " + str(w[0]) + " w2 = " + str(w[1]) + " w3 = " +  str(w[2]) + "\n" + "a1 = " +  str(a[0]) + " a2 = " + str(a[1]) + " a3 = " + str(a[2]) + "\n" "Tempo = " + "T" + " Passo = " + "h" + "\n\n")
output.write("Abaixo os valores de: \nTempo | Posição | Velocidade | Aceleração" + "\n")
print("Resultado do Runge-Kutta-Nystron para as entradas: \nm = " + str(m) + " c = " + str(c) + " k = " + str(k) + "\nw1 = " + str(w[0]) + " w2 =" + str(w[1]) + " w3 = " +  str(w[2]), "a")
print("Abaixo os valores de: \nTempo | Posição | Velocidade | Aceleração")

for i in range(len(results[0])):
    output.write(str(results[0][i]) + " | " + str(results[1][i]) + " | " + str(results[2][i]) + " | " + str(results[3][i]) + "\n")
    print(str(results[0][i]) + " | " + str(results[1][i]) + " | " + str(results[2][i]) + " | " + str(results[3][i]))

output.write("\n\n\n\n")
output.close()