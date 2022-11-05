load('cr_log.sage')

# RR(0x1.fb85251a3f26fp-1)
x = RR("0x1.fb85251a3f26fp-1",16)
hA,lA = cr_log_accurate_path(x,S)
A = RR(hA).exact_rational()+RR(lA).exact_rational()
B = log(x.exact_rational())
C = R200((A-B)/B)
D = -log(abs(C))/log(2.)
print (get_hex(x),get_hex(hA+lA),D)

