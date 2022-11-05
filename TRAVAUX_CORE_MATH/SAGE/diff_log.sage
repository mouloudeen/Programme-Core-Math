load('load.sage')

# Input : x is a 53-bit floating-point numbers 
#         f is an approximation polynomial of log(1+x)
#         G is a list with Coefficient of f decomposed into triple-double 
def diff_cr_log(x,f,G):

    # we calculated len(F)
    F =f.list()
    F.reverse()
    F = F[:-1]
    lenf = len(F)
    # if x> 1.17549435*10^(-38) and x<1 return with the function f 
    if  x> RR('0x1.fba5e353f7ceep-1',16) and x<RR('0x1.ffffffffffffep-1',16)  :
       
       
        
        
        r,(hr,lr) = R200(-1.+x),Add112(-1.0,x)
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,mf,lf = RR(G[0][0],16),RR(G[0][1],16),RR(G[0][2],16)

         # multiply (hr,lr) by (hf,lf) result (h,l)
        hml,(h,m,l) = R200(R200(F[0]).exact_rational()*R200(r).exact_rational()), Mul233(hr,lr,hf,mf,lf)
        for i in range(1,lenf):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,mfi,lfi = RR(F[i][0],16),RR(F[i][1],16),RR(F[i][2],16)

            # add (hfi,lfi) by (h,l) result (h1,l1)
            hml1, (h1,m1,l1)   =  R200(R200(hml).exact_rational() + R200(F[i]).exact_rational()), Add333(hfi,mfi,lfi,h,m,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            hml, (h,m,l) = R200(R200(hml1).exact_rational()*R200(r).exact_rational()), Mul233(hr,lr,h1,m1,l1)
        return R200(hml),h,m,l
    
    
    
    # if x> 1.17549435*10^(-38) and x<1 return with the function f 
    if  x> RR('0x1.0000000000001p+0',16) and x <= RR('0x1.0178e5916f543p+0',16) :
       
       
       
        
        r,(hr,lr) = R200(-1.+x),Add112(x,-1.0)
       
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,mf,lf = RR(G[0][0],16),RR(G[0][1],16),RR(G[0][2],16)

        # multiply (hr,lr) by (hf,lf) result (h,l)
        hml,(h,m,l) = R200(R200(F[0]).exact_rational()*R200(r).exact_rational()), Mul233(hr,lr,hf,mf,lf)
        for i in range(1,lenf):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,mfi,lfi = RR(F[i][0],16),RR(F[i][1],16),RR(F[i][2],16)

            # add (hfi,lfi) by (h,l) result (h1,l1)
            hml1, (h1,m1,l1)   =  R200(R200(hml).exact_rational() + R200(F[i]).exact_rational()), Add333(hfi,mfi,lfi,h,m,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            hml, (h,m,l) = R200(R200(hml1).exact_rational()*R200(r).exact_rational()), Mul233(hr,lr,h1,m1,l1)
        return R200(hml),h,m,l

    #s represents the sign 
    #e  represents the exposant 
    #m represents fraction
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
   
    
    halpha_i,lalpha_i = RR(table_alpha_i[i][0],16),RR(table_alpha_i[i][1],16) #table computed for all i of alpha_i_m
                                      # such that log(alpha_i_m) is 71 bits accurate.
    alpha_i_m = RR(table_alpha_i[i][0],16).exact_rational()+ RR(table_alpha_i[i][1],16).exact_rational()#table calculé pour tous les i d'alpha_i_m
                                      # tel que log(alpha_i_m) soit de précision de 71 bits

    hlog_alpha_i,llog_alpha_i= RR(table_log_alpha_i[i][0],16),RR(table_log_alpha_i[i][1],16)
    log_alpha_i_m= RR(table_log_alpha_i[i][0],16).exact_rational() + RR(table_log_alpha_i[i][1],16).exact_rational()
    
    # multiply (halpha_i_m,lalpha_i_m) by m1
    r,(hr,mr,lr) = R200(m1*alpha_i_m), Mul123(m1,halpha_i,lalpha_i)
   
    # add hr,mr,lr by (-1) in triple-double
    hml,(h,m,l) = R200(r-1), Add133(-1.0,hr,mr,lr)
    
   

    # F[i] is the coefficient of each monomial :x^(len(F)-1-i) in double,double
    
    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,mf,lf = RR(G[0][0],16),RR(G[0][1],16),RR(G[0][2],16)

    # multiply (hr,lr) by (hf,lf) result (h,l)
    hml1,(h1,m1,l1) = R200(R200(hml).exact_rational()*R200(F[0]).exact_rational()), Mul333(h,m,l,hf,mf,lf)
    
    for i in range(1,lenf):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,mfi,lfi = RR(G[i][0],16),RR(G[i][1],16),RR(G[i][2],16)

        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        hml2, (h2,m2,l2) = R200(R200(hml1).exact_rational() + R200(F[i]).exact_rational()), Add333(hfi,mfi,lfi,h1,m1,l1)

        # multiply (h,l) by (h2,l2) result (h1,l1)
        hml1, (h1,m1,l1) = R200(R200(hml2).exact_rational()*R200(hml).exact_rational()),  Mul333(h,m,l,h2,m2,l2)
        
    
    
    # Add (h1,l1) by (hlog_alpha_i,llog_alpha_i
    hml5, (h5,m5,l5) =  R200(R200(hml1).exact_rational() + log_alpha_i_m) , Add233(hlog_alpha_i,llog_alpha_i,h1,m1,l1)
    
    
    # log(2) =(h_log2,l_log2)
    h_log2_e, h_log2 = RR('0x1.62e42fefa38p-1',16).exact_rational(), RR('0x1.62e42fefa38p-1',16)
    l_log2_e, l_log2 = RR('0x1.ef35793c7673p-45',16).exact_rational(),RR('0x1.ef35793c7673p-45',16)

    #e*log(2)
    elog2_h,elog2_l = Mul122(RR(e),h_log2,l_log2) 
   

    fi,(hlog_fi,mlog_fi,llog_fi) =  R200(R200(hml5).exact_rational()) + e*(h_log2_e+l_log2_e) ,Add233(elog2_h,elog2_l,h5,m5,l5)
    
    
    #output: hlog_fi and llog_fi are 53-bit floating-point numbers.
    # hlog_fi: main value and llog_fi: error value
    
    return R200(fi),hlog_fi, mlog_fi, llog_fi

x = RR('0x1.03e90f7831304p+1',16)
a,hA,mA,lA= diff_cr_log(x,u,U)
A = RR(hA).exact_rational()+(RR(mA).exact_rational()+RR(lA).exact_rational())
A1 = R200(a).exact_rational()
B = log(x.exact_rational())
C1 = R200((A1-B)/B)
C = R200((A-B)/B)
D1 = -log(abs(C1))/log(2.)
D = -log(abs(C))/log(2.)
print ("x = %la\n calcule en exact_rational(): (%la , %f)\n calcul en double,double : (%la,%la,%f) "  %(get_hex(x),get_hex(a),D1 , get_hex(hA),get_hex(lA),D))
