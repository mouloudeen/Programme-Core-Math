load('load.sage')

# Input : x is a 53-bit floating-point numbers 
#         f is an approximation polynomial of log(1+x)
def cr_log_fast_path(x,f):

    #If x is NAN return NAN
    if x==NaN:
        return NaN
    
    # If x=+0 or -0 return -INFINITY
    if x == 0:
        return -oo
    
    # If x = +INFINITY return +INFINITY
    if x == +oo:
        return +oo
    
    # if  x<1 and close to 1 return with the function f 
    if  (x> RR('0x1.fe7814e49392fp-1',16) and x<1) : 
       
        # we calculated len(f.list())
        F =f.list()
        F.reverse()
        lenf = len(F)
        
        hr,lr =Add112(-1,x)
        
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,lf = Split(F[0])

        # multiply (hr,lr) by (hf,lf) result (h,l)
        h,l = Mul222(hf,lf,hr,lr)
        
        for i in range(1,lenf-1):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,lfi = Split(F[i])

            # add (hfi,lfi) by (h,l) result (h1,l1)
            h1,l1   = Add222(hfi,lfi,h,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            h,l = Mul222(h1,l1,hr,lr)
            
        return h,l
    
    # if x>1 and close to 1 return with function f
    if x> 1 and x < RR('0x1.00068db8bac71p+0',16) :
         # we calculated len(f.list())
        F =f.list()
        F.reverse()
        lenf = len(F)
        
        hr,lr =Add112(x,-1)
        
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,lf = Split(F[0])

        # multiply (hr,lr) by (hf,lf) result (h,l)
        h,l = Mul222(hf,lf,hr,lr)
        
        for i in range(1,lenf-1):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,lfi = Split(F[i])

            # add (hfi,lfi) by (h,l) result (h1,l1)
            h1,l1   = Add222(hfi,lfi,h,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            h,l = Mul222(h1,l1,hr,lr)
            
        return h,l


    #s represents the sign 
    #e  represents the exposant 
    #m represents fraction
    (s,m,e)=RR(x).sign_mantissa_exponent()
    #e=e-53 
    
    
    e = e+53
    
    
    #If x is a negative number return NAN
    if s==-1 and e!=0:                   
        return NaN
    
    # If x is a subnormal
    if s==1 and e<0 and m!=0:
        
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
    
    
    halpha_i_m,lalpha_i_m = Split(RR(table_alpha_modified[i],16)) #table computed for all i of alpha_i_m
                                      # such that log(alpha_i_m) is 71 bits accurate.
    hlog_alpha_i_m,llog_alpha_i_m= Split(RR(table_log_alpha_modified[i],16)) 
    
    # multiply (halpha_i_m,lalpha_i_m) by (hm,lm)
    hr,lr = Mul122(m1,halpha_i_m,lalpha_i_m)
    
    # add hR,lR by (-1) in double,double
    h,l = Add122(-1,hr,lr)
    
    # we calculated len(f.list())
    F =f.list()
    F.reverse()
    lenf = len(F)

    # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,lf = Split(F[0])

    # multiply (hr,lr) by (hf,lf) result (h,l)
    h1,l1 =Mul222(hf,lf,h,l)

    for i in range(1,lenf-1):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,lfi = Split(F[i])

        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        h2,l2   = Add222(hfi,lfi,h1,l1)

        # multiply (h,l) by (h2,l2) result (h1,l1)
        h1,l1 = Mul222(h2,l2,h,l)
        
    
    
    # Add (h1,l1) by (hlog_alpha_i_m,llog_alpha_i_m)
    h5,l5 = Add222(hlog_alpha_i_m,llog_alpha_i_m,h1,l1)
    
    
    # log(2) =(h_log2,l_log2)
    h_log2 =RR('0x1.62e42fefa38p-1',16)
    l_log2 = RR('0x1.ef35793c7673p-45',16)
    
    #e*log(2)
    elog2_h,elog2_l = Mul122(RR(e),h_log2,l_log2)
   

    hlog_fi,llog_fi = Add222(elog2_h,elog2_l,h5,l5)
    
    
    #output: hlog_fi and llog_fi are 53-bit floating-point numbers.
    # hlog_fi: main value and llog_fi: error value
    
    return RR(hlog_fi),RR(llog_fi)



# Input : x is a 53-bit floating-point numbers 
#         F is a list with Coefficient of f decomposed into double,double 
def cr_log_accurate_path(x,F):

    # if  x<1 and close to 1 return with the function f 
    if  x> RR('0x1.fe73451b9c74fp-1',16) and x<1  :
       
       
        # we calculated len(F)
        lenf = len(F)
        
        hr,lr = Add112(-1,x)
       
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,lf = RR(F[0][0],16),RR(F[0][1],16)

        # multiply (hr,lr) by (hf,lf) result (h,l)
        h,l = Mul222(hr,lr,hf,lf)
       
        for i in range(1,lenf):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,lfi = RR(F[i][0],16),RR(F[i][1],16)

            # add (hfi,lfi) by (h,l) result (h1,l1)
            h1,l1   = Add222(hfi,lfi,h,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            h,l = Mul222(hr,lr,h1,l1)
            
        return h,l

    
    
     # if x>1 and close to 1 return with function f
    if x> 1 and x <= RR('0x1.0178e5916f543p+0',16) :
       
       
        # we calculated len(F)
        lenf = len(F)
        
        hr,lr = Add112(x,-1)
       
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,lf = RR(F[0][0],16),RR(F[0][1],16)

        # multiply (hr,lr) by (hf,lf) result (h,l)
        h,l = Mul222(hr,lr,hf,lf)
       
        for i in range(1,lenf):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,lfi = RR(F[i][0],16),RR(F[i][1],16)

            # add (hfi,lfi) by (h,l) result (h1,l1)
            h1,l1   = Add222(hfi,lfi,h,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            h,l = Mul222(hr,lr,h1,l1)
            
        return h,l

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
    hlog_alpha_i,llog_alpha_i= RR(table_log_alpha_i[i][0],16),RR(table_log_alpha_i[i][1],16)
    
    # multiply (halpha_i_m,lalpha_i_m) by m1
    hr,lr = Mul122(m1,halpha_i,lalpha_i)
    
    # add hr,lr by (-1) in double,double
    h,l = Add122(-1.0,hr,lr)
    
    # we calculated len(f.list())
    lenf = len(F)

    # F[i] is the coefficient of each monomial :x^(len(F)-1-i) in double,double
    
    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,lf = RR(F[0][0],16),RR(F[0][1],16)

    # multiply (hr,lr) by (hf,lf) result (h,l)
    h1,l1 = Mul222(hf,lf,h,l)
    
    for i in range(1,lenf):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,lfi = RR(F[i][0],16),RR(F[i][1],16)

        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        h2,l2   = Add222(hfi,lfi,h1,l1)

        # multiply (h,l) by (h2,l2) result (h1,l1)
        h1,l1 = Mul222(h2,l2,h,l)
        
    
    
    # Add (h1,l1) by (hlog_alpha_i,llog_alpha_i
    h5,l5 = Add222(hlog_alpha_i,llog_alpha_i,h1,l1)
    
    
    # log(2) =(h_log2,l_log2)
    h_log2 =RR('0x1.62e42fefa39efp-1',16)
    l_log2 = RR('0x1.abc9e3b39803fp-56',16)
    
    #e*log(2)
    elog2_h,elog2_l = Mul122(RR(e),h_log2,l_log2)
   

    hlog_fi,llog_fi = Add222(elog2_h,elog2_l,h5,l5)
    
    
    #output: hlog_fi and llog_fi are 53-bit floating-point numbers.
    # hlog_fi: main value and llog_fi: error value
    
    return hlog_fi,llog_fi



# Input : x is a 53-bit floating-point numbers 
#         F is a list with Coefficient of f decomposed into triple-double 
def cr_log_accurate_path_advanced(x,F):

    # if x> 1.17549435*10^(-38) and x<1 return with the function f 
    '''if  x> RR('0x1.fba5e353f7ceep-1',16) and x<1  :
       
       
        # we calculated len(F)
        lenf = len(F)
        
        hr,lr = Add112(-1.0,x)
        print("hr = %la, lr = %la " %(get_hex(hr),get_hex(lr)))
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,mf,lf = RR(F[0][0],16),RR(F[0][1],16),RR(F[0][2],16)

        # multiply (hr,lr) by (hf,lf) result (h,l)
        h,m,l = Mul233(hr,lr,hf,mf,lf)
        print("i = 0, h = %la,m = %la, l = %la " %(get_hex(h),get_hex(m),get_hex(l)))
        for i in range(1,lenf):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,mfi,lfi = RR(F[i][0],16),RR(F[i][1],16),RR(F[i][2],16)

            # add (hfi,lfi) by (h,l) result (h1,l1)
            h1,m1,l1   = Add333(hfi,mfi,lfi,h,m,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            h,m,l = Mul233(hr,lr,h1,m1,l1)
            ("i = %d, h1 = %la,m1 =%la, l1 = %la, h = %la,m = %la,l = %la " %(i,get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h),get_hex(m),get_hex(l)))
        return h,m,l
    
    
    
    # if x> 1.17549435*10^(-38) and x<1 return with the function f 
    if  x> 1 and x <= RR('0x1.0178e5916f543p+0',16) :
       
       
        # we calculated len(F)
        lenf = len(F)
        
        hr,lr = Add112(x,-1.0)
        print("hr = %la, lr = %la " %(get_hex(hr),get_hex(lr)))
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,mf,lf = RR(F[0][0],16),RR(F[0][1],16),RR(F[0][2],16)

        # multiply (hr,lr) by (hf,lf) result (h,l)
        h,m,l = Mul233(hr,lr,hf,mf,lf)
        print("i = 0, h = %la,m = %la, l = %la " %(get_hex(h),get_hex(m),get_hex(l)))
        for i in range(1,lenf):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,mfi,lfi = RR(F[i][0],16),RR(F[i][1],16),RR(F[i][2],16)

            # add (hfi,lfi) by (h,l) result (h1,l1)
            h1,m1,l1   = Add333(hfi,mfi,lfi,h,m,l)

            # multiply (hr,lr) by (h1,l1) result (h,l)
            h,m,l = Mul233(hr,lr,h1,m1,l1)
            print("i = %d, h1 = %la,m1 =%la, l1 = %la, h = %la,m = %la,l = %la " %(i,get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h),get_hex(m),get_hex(l)))
        return h,m,l'''

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
   
    
    halpha_i,malpha_i,lalpha_i = RR(table_alpha_i_triple[i][0],16),RR(table_alpha_i_triple[i][1],16),RR(table_alpha_i_triple[i][2],16) #table computed for all i of alpha_i_m
                                      # such that log(alpha_i_m) is 71 bits accurate.
    hlog_alpha_i,mlog_alpha_i,llog_alpha_i= RR(table_log_alpha_i_triple[i][0],16),RR(table_log_alpha_i_triple[i][1],16),RR(table_log_alpha_i_triple[i][2],16)
    
    # multiply (halpha_i_m,lalpha_i_m) by m1
    hr,mr,lr = Mul133(m1,halpha_i,malpha_i,lalpha_i)
   # print("m1 = %la  halpha_i = %la , malpha_i = %la ,lalpha_i = %la ,hr = %la,mr = %la,lr = %la" %(get_hex(m1),get_hex(halpha_i),get_hex(malpha_i), get_hex(lalpha_i),get_hex(hr),get_hex(mr),get_hex(lr)))
    # add hr,mr,lr by (-1) in triple-double
    h,m,l = Add133(-1.0,hr,mr,lr)
    #print("h = %la, m = %la, l = %la" %(get_hex(h),get_hex(m),get_hex(l)))
    # we calculated len(f.list())
    lenf = len(F)

    # F[i] is the coefficient of each monomial :x^(len(F)-1-i) in double,double
    
    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,mf,lf = RR(F[0][0],16),RR(F[0][1],16),RR(F[0][2],16)

    # multiply (hr,lr) by (hf,lf) result (h,l)
    h1,m1,l1 = Mul333(h,m,l,hf,mf,lf)
    #print("h1 = %la, m1 = %la, l1 = %la" %(get_hex(h1),get_hex(m1),get_hex(l1)))
    for i in range(1,lenf):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,mfi,lfi = RR(F[i][0],16),RR(F[i][1],16),RR(F[i][2],16)

        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        h2,m2,l2   = Add333(hfi,mfi,lfi,h1,m1,l1)

        # multiply (h,l) by (h2,l2) result (h1,l1)
        h1,m1,l1 = Mul333(h,m,l,h2,m2,l2) 
        #print("i = %d, h1 = %la, m1 = %la, l1 = %la, h2 = %la, m2 = %la, l2 = %la" %(i,get_hex(h1),get_hex(m1),get_hex(l1),get_hex(h2),get_hex(m2),get_hex(l2)))
    
    
    # Add (h1,l1) by (hlog_alpha_i,llog_alpha_i
    h5,m5,l5 = Add333(hlog_alpha_i,mlog_alpha_i,llog_alpha_i,h1,m1,l1)
    #print("h5 = %la , m5 = %la, l5 = %la , hlog_alpha_i = %la, mlog_alpha_i = %la, llog_alpha_i = %la" %(get_hex(h5),get_hex(m5),get_hex(l5),get_hex(hlog_alpha_i),get_hex(mlog_alpha_i),get_hex(llog_alpha_i)))
    
    # log(2) =(h_log2,l_log2)
    h_log2 =RR('0x1.62e42fefa39efp-1',16)
    m_log2 =RR('0x1.abc9e3b39803fp-56',16)
    l_log2 = RR('0x1.7b57a079a1934p-111',16)  
    
    #e*log(2)
    elog2_h,elog2_m,elog2_l = Mul133(RR(e),h_log2,m_log2,l_log2)
   

    hlog_fi,mlog_fi,llog_fi = Add333(elog2_h,elog2_m,elog2_l,h5,m5,l5)
    
    
    #output: hlog_fi and llog_fi are 53-bit floating-point numbers.
    # hlog_fi: main value and llog_fi: error value
    #print("h6 = %la, m6 = %la, l6 =%la" %(get_hex(hlog_fi),get_hex(mlog_fi),get_hex(llog_fi)))
    return hlog_fi, mlog_fi, llog_fi
