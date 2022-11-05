load('load.sage')
x = RR('0x1.93e17ed15431ap-283',16)
print("x = %la\n" %(get_hex(x)))
(s,m,e)=RR(x).sign_mantissa_exponent()
#e=e-53 
e = e+53

# If x is a subnormal
if s==1 and e<=0 and m!=0:

    v = m                             
    e = e-1

    while v< 2^52:
        v*=2

        e=e-1
    m1 =v*2.^(-52)

#If x is normal
else:
    #print(s,e,m)
    m1 = m*2.^(-52)
    e =e-1

binary = RR(m1).str(2)[2:10] # we recover the 8 bits after the initial 1.
i = int(binary,2)         # i is the 8-bit integer
print("entier du binaire = %d " %(i))
log_m = R200(log(m1)).exact_rational()
halpha_i,malpha_i,lalpha_i = RR(table_alpha_i_triple[i][0],16),RR(table_alpha_i_triple[i][1],16),RR(table_alpha_i_triple[i][2],16) #table computed for all i of alpha_i_m
                                  # such that log(alpha_i_m) is 71 bits accurate.
hlog_alpha_i,mlog_alpha_i,llog_alpha_i= RR(table_log_alpha_i_triple[i][0],16),RR(table_log_alpha_i_triple[i][1],16),RR(table_log_alpha_i_triple[i][2],16)
log_mf = u(m1*(halpha_i.exact_rational()+malpha_i.exact_rational()+lalpha_i.exact_rational())-1)+(hlog_alpha_i.exact_rational()+mlog_alpha_i.exact_rational()+llog_alpha_i.exact_rational())
# multiply (halpha_i_m,lalpha_i_m) by m1
hr,mr,lr = Mul133(m1,halpha_i,malpha_i,lalpha_i)
print("m1 = %la hr = %la, mr = %la, lr = %la , halpha_i = %la, malpha_i = %la, lalpha_i = %la" %(get_hex(m1),get_hex(hr),get_hex(mr),get_hex(lr),get_hex(halpha_i),get_hex(malpha_i),get_hex(lalpha_i)))
ti = R200(m1).exact_rational()*(R200(halpha_i).exact_rational()+R200(malpha_i).exact_rational()+R200(lalpha_i).exact_rational())
tf = R200(hr).exact_rational()+(R200(mr).exact_rational()+R200(lr).exact_rational())
print("L'erreur de Mul133(m1,halpha_i,lalpha_i) est : \n",R200((ti-tf)/ti))
# add hr,mr,lr by (-1) in triple-double
h,m,l = Add133(-1.0,hr,mr,lr)
print("h = %la, m = %la, l = %la " %(get_hex(h),get_hex(m),get_hex(l)))
ti = R200(hr).exact_rational()+(R200(mr).exact_rational()+R200(lr).exact_rational())-1
tf = R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational())                                
print("L'erreur de Add133(-1.0,hr,mr,lr) est : \n",R200((ti-tf)/ti))
# we calculated len(f.list())
lenf = len(U)

# F[i] is the coefficient of each monomial :x^(len(F)-1-i) in double,double

#We are turning F[0] in hf,lf with hf,: main value
#                                  lf: error value
hf,mf,lf = RR(U[0][0],16),RR(U[0][1],16),RR(U[0][2],16)

# multiply (hr,lr) by (hf,lf) result (h,l)
h1,m1,l1 = Mul333(hf,mf,lf,h,m,l)
print("i = 0, h1 = %la, m1 = %la, l1 = %la " %(get_hex(h1),get_hex(m1),get_hex(l1)))
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(hf).exact_rational()+(R200(mf).exact_rational()+R200(lf).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("L'erreur de Mul333(hf,mf,lf,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[1][0],16),RR(U[1][1],16),RR(U[1][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 1 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[2][0],16),RR(U[2][1],16),RR(U[2][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 2 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[3][0],16),RR(U[3][1],16),RR(U[3][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 3 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[4][0],16),RR(U[4][1],16),RR(U[4][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 4 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[5][0],16),RR(U[5][1],16),RR(U[5][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 5 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[6][0],16),RR(U[6][1],16),RR(U[6][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 6 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[7][0],16),RR(U[7][1],16),RR(U[7][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 7 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[8][0],16),RR(U[8][1],16),RR(U[8][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 8 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[9][0],16),RR(U[9][1],16),RR(U[9][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 9 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[10][0],16),RR(U[10][1],16),RR(U[10][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 10 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[11][0],16),RR(U[11][1],16),RR(U[11][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 11 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[12][0],16),RR(U[12][1],16),RR(U[12][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 12 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

#We are turning F[i] in hfi,lfi with hfi: main value
#                                    lfi: error value
hfi,mfi,lfi = RR(U[13][0],16),RR(U[13][1],16),RR(U[13][2],16)

# add (hfi,lfi) by (h1,l1) result (h2,l2)
h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)
ti1 = (R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
tf1 = R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational())
# multiply (h,l) by (h2,l2) result (h1,l1)
h1,m1,l1 = Mul333(h2,m2,l2,h,m,l)
ti = (R200(h).exact_rational()+(R200(m).exact_rational()+R200(l).exact_rational()))*(R200(h2).exact_rational()+(R200(m2).exact_rational()+R200(l2).exact_rational()))
tf = R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational())
print("i = 13 , h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
print("L'erreur de Add333(hfi,mfi,lfi,h1,m1,l1) est : \n",R200((ti1-tf1)/ti1))
print("L'erreur de Mul333(h2,m2,l2,h,m,l) est : \n",R200((ti-tf)/ti))

# Add (h1,l1) by (hlog_alpha_i,llog_alpha_i
h5,m5,l5 = Add333(hlog_alpha_i,mlog_alpha_i,llog_alpha_i,h1,m1,l1)
print("h5 = %la, m5 = %la, l5 = %la , hlog_alpha_i = %la, llog_alpha_i = %la" %(get_hex(h5),get_hex(m5),get_hex(l5),get_hex(hlog_alpha_i),get_hex(llog_alpha_i)))
ti = (R200(hlog_alpha_i).exact_rational()+R200(mlog_alpha_i).exact_rational()+R200(llog_alpha_i).exact_rational())+(R200(h1).exact_rational()+(R200(m1).exact_rational()+R200(l1).exact_rational()))
tf = R200(h5).exact_rational()+(R200(m5).exact_rational()+R200(l5).exact_rational())
print("l'erreur de Add333(hlog_alpha_i,mlog_alpha_i,llog_alpha_i,h1,m1,l1) est :" ,R200((ti-tf)/ti))
Clog = R200((log_m-tf)/tf)
Dlog = R200(-log(abs(Clog))/log(2.0))
Cm = R200((log_m-log_mf)/log_m)
Dm = R200(-log(abs(Cm))/log(2.0))
print("La précison du log(m) est ",Dlog)
print("La précison du log(m) avec u est ",Dm)
# log(2) =(h_log2,l_log2)
h_log2 =RR('0x1.62e42fefa39efp-1',16)
m_log2 =RR('0x1.abc9e3b39803fp-56',16)
l_log2 = RR('0x1.7b57a079a1934p-111',16)  

#e*log(2)
elog2_h,elog2_m,elog2_l = Mul133(RR(e),h_log2,m_log2,l_log2) 
print("elog2_h = %la, elog2_m = %la, elog2_l = %la e = %d" %(get_hex(elog2_h),get_hex(elog2_m),get_hex(elog2_l),e))
ti = R200(e).exact_rational()*(R200(h_log2).exact_rational()+R200(m_log2).exact_rational()+R200(l_log2).exact_rational())
tf = R200(elog2_h).exact_rational()+ R200(elog2_m).exact_rational()+ R200(elog2_l).exact_rational()
print("l'erreur de Mul133(RR(e),h_log2,m_log2,l_log2) est : ", R200((ti-tf)/ti))
te = R200(e).exact_rational()*R200(log(2)).exact_rational()
Ce = R200((te-tf)/te)
De = R200(-log(abs(Ce))/log(2.0))
print("la précison de e*log(2) est ", De)
hlog_fi,mlog_fi,llog_fi = Add333(elog2_h,elog2_m,elog2_l,h5,m5,l5)

print("hlog_fi = %la, mlog_fi = %la, llog_fi = %la" %(get_hex(hlog_fi),get_hex(mlog_fi),get_hex(llog_fi)))
ti = (R200(elog2_h).exact_rational()+R200(elog2_m).exact_rational()+R200(elog2_l).exact_rational())+(R200(h5).exact_rational()+(R200(m5).exact_rational()+R200(l5).exact_rational()))
tf = R200(hlog_fi).exact_rational()+(R200(mlog_fi).exact_rational()+R200(llog_fi).exact_rational())
print("L'erreur de  Add333(elog2_h,elog2_m,elog2_l,h5,m5,l5) est :",R200((ti-tf)/ti))
A = RR(hlog_fi).exact_rational()+(RR(mlog_fi).exact_rational()+RR(llog_fi).exact_rational())
B = log(x.exact_rational())
C = R200((A-B)/B)
D =R200(-log(abs(C))/log(2.))
print("Le nombre de bits de précisions : %f" %(D))