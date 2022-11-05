load('load.sage')
x = RR('0x1.fb85251a3f26fp-1',16)
print("x = %la\n" %(get_hex(x)))


# if  x<1 and close to 1 return with the function f 
if  x> RR('0x1.fba5e353f7ceep-1',16) and x<1  :
       
       
    # we calculated len(F)
    lenf = len(S)
    
    hr,lr = Add112(-1,x)
    print("hr = %la, lr = %la " %(get_hex(hr),get_hex(lr)))
    ti = R200(x).exact_rational()+R200(-1.).exact_rational()
    tf = R200(hr).exact_rational()+R200(lr).exact_rational()
    print("L'erreur de Add112(x,-1) est : \n",R200((ti-tf)/ti))
    
    # F[i] is the coefficient of each monomial :x^(len(F)-1-i)

    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,lf = RR(S[0][0],16),RR(S[0][1],16)

    # multiply (hr,lr) by (hf,lf) result (h,l)
    h,l = Mul222(hr,lr,hf,lf)
    print("i = 0, h = %la,m = %la, l = %la " %(get_hex(h),get_hex(m),get_hex(l)))
    ti = (R200(hr).exact_rational()+R200(lr).exact_rational())*(R200(hf).exact_rational()+R200(lf).exact_rational())
    tf = R200(h).exact_rational()+R200(l).exact_rational()
    print("L'erreur de Mul222(hr,lr,hf,lf) est : \n",R200((ti-tf)/ti))
    
    for i in range(1,lenf):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,lfi = RR(S[i][0],16),RR(S[i][1],16)

        # add (hfi,lfi) by (h,l) result (h1,l1)
        h1,l1   = Add222(hfi,lfi,h,l)
        ti1 = (R200(h).exact_rational()+R200(l).exact_rational())+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
        tf1 = R200(h1).exact_rational()+R200(l1).exact_rational()

        # multiply (hr,lr) by (h1,l1) result (h,l)
        h,l = Mul222(hr,lr,h1,l1)
        ti = (R200(hr).exact_rational()+R200(lr).exact_rational())*(R200(h1).exact_rational()+R200(l1).exact_rational())
        tf = R200(h1).exact_rational()+R200(l1).exact_rational()
        print("i = %d, h1 = %la,m1 =%la, l1 = %la, h = %la,m = %la,l = %la " %(i,get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h),get_hex(m),get_hex(l)))
        print("L'erreur de Add222(hfi,lfi,h,l) est : \n",R200((ti1-tf1)/ti1))
        print("L'erreur de Mul222(hr,lr,h1,l1) est : \n",R200((ti-tf)/ti))
 
# if x>1 and close to 1 return with function f
if x> 1 and x <= RR('0x1.0178e5916f543p+0',16) :
       
       
    # we calculated len(F)
    lenf = len(S)
    
    hr,lr = Add112(x,-1)
    print("hr = %la, lr = %la " %(get_hex(hr),get_hex(lr)))
    ti = R200(x).exact_rational()+R200(-1.).exact_rational()
    tf = R200(hr).exact_rational()+R200(lr).exact_rational()
    print("L'erreur de Add112(x,-1) est : \n",R200((ti-tf)/ti))
    
    # F[i] is the coefficient of each monomial :x^(len(F)-1-i)

    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,lf = RR(S[0][0],16),RR(S[0][1],16)

    # multiply (hr,lr) by (hf,lf) result (h,l)
    h,l = Mul222(hr,lr,hf,lf)
    print("i = 0, h = %la,m = %la, l = %la " %(get_hex(h),get_hex(m),get_hex(l)))
    ti = (R200(hr).exact_rational()+R200(lr).exact_rational())*(R200(hf).exact_rational()+R200(lf).exact_rational())
    tf = R200(h).exact_rational()+R200(l).exact_rational()
    print("L'erreur de Mul222(hr,lr,hf,lf) est : \n",R200((ti-tf)/ti))

    for i in range(1,lenf):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,lfi = RR(S[i][0],16),RR(S[i][1],16)

        # add (hfi,lfi) by (h,l) result (h1,l1)
        h1,l1   = Add222(hfi,lfi,h,l)
        ti1 = (R200(h).exact_rational()+R200(l).exact_rational())+((R200(hfi).exact_rational()+(R200(mfi).exact_rational()+R200(lfi).exact_rational())))
        tf1 = R200(h1).exact_rational()+R200(l1).exact_rational()

        # multiply (hr,lr) by (h1,l1) result (h,l)
        h,l = Mul222(hr,lr,h1,l1)
        ti = (R200(hr).exact_rational()+R200(lr).exact_rational())*(R200(h1).exact_rational()+R200(l1).exact_rational())
        tf = R200(h1).exact_rational()+R200(l1).exact_rational()
        print("i = %d, h1 = %la,m1 =%la, l1 = %la, h = %la,m = %la,l = %la " %(i,get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h),get_hex(m),get_hex(l)))
        print("L'erreur de Add222(hfi,lfi,h,l) est : \n",R200((ti1-tf1)/ti1))
        print("L'erreur de Mul222(hr,lr,h1,l1) est : \n",R200((ti-tf)/ti))
   
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
halpha_i,lalpha_i = RR(table_alpha_i[i][0],16),RR(table_alpha_i[i][1],16) #table computed for all i of alpha_i_m
                                      # such that log(alpha_i_m) is 71 bits accurate.
hlog_alpha_i,llog_alpha_i= RR(table_log_alpha_i[i][0],16),RR(table_log_alpha_i[i][1],16)

log_mf = sa(m1*(halpha_i.exact_rational()+lalpha_i.exact_rational())-1.)+(hlog_alpha_i.exact_rational()+llog_alpha_i.exact_rational())

# multiply (halpha_i_m,lalpha_i_m) by m1
hr,lr = Mul122(m1,halpha_i,lalpha_i)
print("m1 = %la hr = %la, lr = %la , halpha_i = %la, lalpha_i = %la" %(get_hex(m1),get_hex(hr),get_hex(lr),get_hex(halpha_i),get_hex(lalpha_i)))
ti = R200(m1).exact_rational()*(R200(halpha_i).exact_rational()+R200(lalpha_i).exact_rational())
tf = R200(hr).exact_rational()+R200(lr).exact_rational()
print("L'erreur de Mul122(m1,halpha_i,lalpha_i) est : \n",R200((ti-tf)/ti))

# add hr,lr by (-1) in double-double
h,l = Add122(-1.0,hr,lr)
print("h = %la, l = %la " %(get_hex(h),get_hex(l)))
ti = (R200(hr).exact_rational()+R200(lr).exact_rational())-1
tf = R200(h).exact_rational()+R200(l).exact_rational()                               
print("L'erreur de Add133(-1.0,hr,mr,lr) est : \n",R200((ti-tf)/ti))

# we calculated len(f.list())
lenf = len(S)

# F[i] is the coefficient of each monomial :x^(len(F)-1-i) in double,double

#We are turning F[0] in hf,lf with hf,: main value
#                                  lf: error value
hf,lf = RR(S[0][0],16),RR(S[0][1],16)

# multiply (hr,lr) by (hf,lf) result (h,l)
h1,l1 = Mul222(hf,lf,h,l)
print("i = 0, h1 = %la, l1 = %la " %(get_hex(h1),get_hex(l1)))
ti = (R200(h).exact_rational()+R200(l).exact_rational())*(R200(hf).exact_rational()++R200(lf).exact_rational())
tf = R200(h1).exact_rational()+R200(l1).exact_rational()
print("L'erreur de Mul222(hf,lf,h,l) est : \n",R200((ti-tf)/ti))

for i in range(1,lenf):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,lfi = RR(S[i][0],16),RR(S[i][1],16)

        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        h2,l2   = Add222(hfi,lfi,h1,l1)
        ti1 = (R200(h1).exact_rational()+R200(l1).exact_rational())+((R200(hfi).exact_rational()+R200(lfi).exact_rational()))
        tf1 = R200(h2).exact_rational()+R200(l2).exact_rational()

        # multiply (h,l) by (h2,l2) result (h1,l1)
        h1,l1 = Mul222(h2,l2,h,l)
        ti = (R200(h).exact_rational()+R200(l).exact_rational())*(R200(h2).exact_rational()+R200(l2).exact_rational())
        tf = R200(h1).exact_rational()+R200(l1).exact_rational()
        print("i = %d , h1 = %la, l1 = %la, h2 = %la, l2 = %la" %(i,get_hex(h1),get_hex(l1),get_hex(h2),get_hex(l2)))
        print("L'erreur de Add222(hfi,lfi,h1,l1) est : \n",R200((ti1-tf1)/ti1))
        print("L'erreur de Mul222(h2,l2,h,l) est : \n",R200((ti-tf)/ti))


# Add (h1,l1) by (hlog_alpha_i,llog_alpha_i
h5,l5 = Add222(hlog_alpha_i,llog_alpha_i,h1,l1)
print("h5 = %la, l5 = %la , hlog_alpha_i = %la, llog_alpha_i = %la" %(get_hex(h5),get_hex(l5),get_hex(hlog_alpha_i),get_hex(llog_alpha_i)))
ti = (R200(hlog_alpha_i).exact_rational()+R200(llog_alpha_i).exact_rational())+(R200(h1).exact_rational()++R200(l1).exact_rational())
tf = R200(h5).exact_rational()+R200(l5).exact_rational()
print("l'erreur de Add222(hlog_alpha_i,llog_alpha_i,h1,l1) est :" ,R200((ti-tf)/ti))

Clog = R200((log_m-tf)/tf)
Dlog = R200(-log(abs(Clog))/log(2.0))

Cm = R200((log_m-log_mf)/log_m)
Dm = R200(-log(abs(Cm))/log(2.0))

print("La précison du log(m) est ",Dlog)
print("La précison du log(m) avec s est ",Dm)

# log(2) =(h_log2,l_log2)
h_log2 =RR('0x1.62e42fefa39efp-1',16)
l_log2 = RR('0x1.abc9e3b39803fp-56',16)

#e*log(2)
elog2_h,elog2_l = Mul122(RR(e),h_log2,l_log2)
print("elog2_h = %la,  elog2_l = %la e = %d" %(get_hex(elog2_h),get_hex(elog2_l),e))
ti = R200(e).exact_rational()*(R200(h_log2).exact_rational()+R200(l_log2).exact_rational())
tf = R200(elog2_h).exact_rational()+ R200(elog2_l).exact_rational()
print("l'erreur de Mul122(RR(e),h_log2,l_log2) est : ", R200((ti-tf)/ti))

te = R200(e).exact_rational()*R200(log(2)).exact_rational()
Ce = R200((te-tf)/te)
De = R200(-log(abs(Ce))/log(2.0))
print("la précison de e*log(2) est ", De)
hlog_fi,llog_fi = Add222(elog2_h,elog2_l,h5,l5)

print("hlog_fi = %la, llog_fi = %la" %(get_hex(hlog_fi),get_hex(llog_fi)))
ti = (R200(elog2_h).exact_rational()+R200(elog2_l).exact_rational())+(R200(h5).exact_rational()+R200(l5).exact_rational())
tf = R200(hlog_fi).exact_rational()+R200(llog_fi).exact_rational()
print("L'erreur de  Add222(elog2_h,elog2_l,h5,l5) est :",R200((ti-tf)/ti))
A = RR(hlog_fi).exact_rational()+RR(llog_fi).exact_rational()
B = log(x.exact_rational())
C = R200((A-B)/B)
D =R200(-log(abs(C))/log(2.))
print("Le nombre de bits de précisions : %f" %(D))