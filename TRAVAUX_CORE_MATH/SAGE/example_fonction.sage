load('load.sage')
print("Test avec x = R107(x)\n")
L = []
L1 = []
j =0
z =0
for i in range(1000):
    x = random()
    while x> 1/2^8 or x <  0.00001:
        x = random()
    x = R107(x)
    A = R200(log((1+x).exact_rational()))
    B = R200(p(x))
    C = R200((B.exact_rational()-A.exact_rational())/A.exact_rational())
    D = R200(-log(abs(C))/log(2.0))
    if D < 104:
        L.append(x)
        j=j+1
    if D > 104:
        L1.append(x)
        z=z+1
print("le nombre de résultat qui sont superieur à 104 bits :\n",z)
print("le nombre de résultat qui sont inferieur à 104 bits :\n",j)

print("Test avec x =RR(x)\n")
L = []
L1 = []
j =0
z =0
for i in range(1000):
    x = random()
    while x> 1/2^8 or x <  0.00001:
        x = random()
    x = RR(x)
    A = R200(log((1+x).exact_rational()))
    B = R200(p(x))
    C = R200((B.exact_rational()-A.exact_rational())/A.exact_rational())
    D = R200(-log(abs(C))/log(2.0))
    if D < 104:
        L.append(x)
        j=j+1
    if D > 104:
        L1.append(x)
        z=z+1
print("le nombre de résultat qui sont superieur à 104 bits :\n",z)
print("le nombre de résultat qui sont inferieur à 104 bits :\n",j)