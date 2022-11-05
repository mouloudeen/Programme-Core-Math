load('cr_log.sage')
load('cr_log_modif.sage')

L = []
nb =1000000
n =0
for i in range(nb):
    if i==0:
       x = RR('0x1.000413d13ed1ap+0',16)
    else:
       x = randint(1,1000)+random()
       x =RR(x)
    hA,lA = cr_log_accurate_path(x,G)
    A = RR(hA).exact_rational()+RR(lA).exact_rational()
    B = log(x.exact_rational())
    C = R200((A-B)/B)
    D = -log(abs(C))/log(2.)
    if D<101:
        print("x = %la ,D = %d ,n = %d, i = %d \n" %(get_hex(x),D,n,i))
        n+=1
print("nombre d'erreur est ",n)
