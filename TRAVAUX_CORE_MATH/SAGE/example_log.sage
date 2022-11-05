load('cr_log.sage')

nb = 1000
L = []
# We test 1000 values  between 1 and 300
for i in range(1000):
    x = randint(1,300)+random()
    x1 = RR(x)
    hA,lA = cr_log_fast_path(x1,f)
    hA1,lA1 = cr_log_accurate_path(x1,G)
    A = RR(hA).exact_rational()+RR(lA).exact_rational()
    A1= RR(hA1).exact_rational()+RR(lA1).exact_rational()
    B = log(x1.exact_rational())
    C = R200((A-B)/B)
    C1 = R200((A1-B)/B)
    D = -log(abs(C))/log(2.)
    D1 = -log(abs(C1))/log(2.)
    L.append((x,D,D1))
for i in range(nb): 
    print(L[i])