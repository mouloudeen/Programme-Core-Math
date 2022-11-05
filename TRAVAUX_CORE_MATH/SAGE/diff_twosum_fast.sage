load('load.sage')
 
def test_two_fast(x,f):
 
    # if x> 1.17549435*10^(-38) and x<1 return with the function f 
    if  (x> RR('0x1.fe7814e49392fp-1',16) and x<1) or (x> 1 and x < RR('0x1.00068db8bac71p+0',16)) :
       
        # we calculated len(f.list())
        F =f.list()
        F.reverse()
        lenf = len(F)
        
        (hr1,lr1),(hr2,lr2) =TwoSum(x,-1),Add12(x,-1)
        t1 = R200(((R200(hr2).exact_rational()+R200(lr2).exact_rational())-(R200(hr1).exact_rational()+R200(lr1).exact_rational()))/(R200(hr1).exact_rational()+R200(lr1).exact_rational()))
        print("hr1 = %la, hr2 = %la, lr1 = %la, lr2 = %la, diff = %f\n" %(get_hex(hr1),get_hex(hr2),get_hex(lr1),get_hex(lr2),t1))
        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
        #We are turning F[0] in hf,lf with hf,: main value
        #                                  lf: error value
        hf,lf = Split(F[0])

        # multiply (hr,lr) by (hf,lf) result (h,l)
        (ha,la),(hb,lb) =prod_dekker3(hf,lf,hr1,lr1),prod_dekker3(hf,lf,hr2,lr2)
        t2 = R200(((R200(hb).exact_rational()+R200(lb).exact_rational())-(R200(ha).exact_rational()+R200(la).exact_rational()))/(R200(ha).exact_rational()+R200(la).exact_rational()))
        print("ha = %la, hb = %la, la = %la, lb = %la, diff = %f\n" %(get_hex(ha),get_hex(hb),get_hex(la),get_hex(lb),t2))
        for i in range(1,lenf-1):
            #We are turning F[i] in hfi,lfi with hfi: main value
            #                                    lfi: error value
            hfi,lfi = Split(F[i])
            print("i = \n",i)
            # add (hfi,lfi) by (h,l) result (h1,l1)
            (h1a,l1a),(h1b,l1b)   = TwoSum_m2(hfi,lfi,ha,la),TwoSum_m2(hfi,lfi,hb,lb)
            t3 = R200(((R200(h1b).exact_rational()+R200(l1b).exact_rational())-(R200(h1a).exact_rational()+R200(l1a).exact_rational()))/(R200(h1a).exact_rational()+R200(l1a).exact_rational()))
            print("h1a = %la, h1b = %la, l1a = %la, l1b = %la, diff = %f\n" %(get_hex(h1a),get_hex(h1b),get_hex(l1a),get_hex(l1b),t3))

            # multiply (hr,lr) by (h1,l1) result (h,l)
            (ha,la),(hb,lb) = prod_dekker3(h1a,l1a,hr1,lr1),prod_dekker3(h1b,l1b,hr2,lr2)
            t4 = R200(((R200(hb).exact_rational()+R200(lb).exact_rational())-(R200(ha).exact_rational()+R200(la).exact_rational()))/(R200(ha).exact_rational()+R200(la).exact_rational()))
            print("ha = %la, hb = %la, la = %la, lb = %la, diff = %f\n" %(get_hex(ha),get_hex(hb),get_hex(la),get_hex(lb),t4))
            
        return (ha,la),(hb,lb)
    
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
    
    hm,lm = Split(m1)
    
    
    binary = RR(m1).str(2)[2:10] # we recover the 8 bits after the initial 1.
    i = int(binary,2)         # i is the 8-bit integer
    
    
    halpha_i_m,lalpha_i_m = Split(RR(table_alpha_modified[i],16)) #table computed for all i of alpha_i_m
                                      # such that log(alpha_i_m) is 71 bits accurate.
    hlog_alpha_i_m,llog_alpha_i_m= Split(RR(table_log_alpha_modified[i],16)) 
    
    # multiply (halpha_i_m,lalpha_i_m) by (hm,lm)
    hr,lr = prod_dekker3(hm,lm,halpha_i_m,lalpha_i_m)
    
    # add hR,lR by (-1) in double,double
    (ha,la),(hb,lb) = TwoSum_m(-1,hr,lr),Add122(-1,hr,lr)
    t5 = R200(((R200(hb).exact_rational()+R200(lb).exact_rational())-(R200(ha).exact_rational()+R200(la).exact_rational()))/(R200(ha).exact_rational()+R200(la).exact_rational()))
    print("ha = %la, hb = %la, la = %la, lb = %la, diff = %f\n" %(get_hex(ha),get_hex(hb),get_hex(la),get_hex(lb),t5))
    # we calculated len(f.list())
    F =f.list()
    F.reverse()
    lenf = len(F)

    # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
    #We are turning F[0] in hf,lf with hf,: main value
    #                                  lf: error value
    hf,lf = Split(F[0])

    # multiply (hr,lr) by (hf,lf) result (h,l)
    (h1a,l1a),(h1b,l1b) =prod_dekker3(hf,lf,ha,la),prod_dekker3(hf,lf,hb,lb),
    t6 = R200(((R200(h1b).exact_rational()+R200(l1b).exact_rational())-(R200(h1a).exact_rational()+R200(l1a).exact_rational()))/(R200(h1a).exact_rational()+R200(l1a).exact_rational()))
    print("h1a = %la, h1b = %la, l1a = %la, l1b = %la, diff = %f\n" %(get_hex(h1a),get_hex(h1b),get_hex(l1a),get_hex(l1b),t6))
    for i in range(1,lenf-1):
        #We are turning F[i] in hfi,lfi with hfi: main value
        #                                    lfi: error value
        hfi,lfi = Split(F[i])
        print("i =\n",i)
        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        (h2a,l2a),(h2b,l2b)   = TwoSum_m2(hfi,lfi,h1a,l1a),Add22(hfi,lfi,h1b,l1b)
        t7 = R200(((R200(h2b).exact_rational()+R200(l2b).exact_rational())-(R200(h2a).exact_rational()+R200(l2a).exact_rational()))/(R200(h2a).exact_rational()+R200(l2a).exact_rational()))
        print("h2a = %la, h2b = %la, l2a = %la, l2b = %la, diff = %f\n" %(get_hex(h2a),get_hex(h2b),get_hex(l2a),get_hex(l2b),t7))
        # multiply (h,l) by (h2,l2) result (h1,l1)
        (h1a,l1a),(h1b,l1b) = prod_dekker3(h2a,l2a,ha,la), prod_dekker3(h2b,l2b,hb,lb)
        t8 = R200(((R200(h1b).exact_rational()+R200(l1b).exact_rational())-(R200(h1a).exact_rational()+R200(l1a).exact_rational()))/(R200(h1a).exact_rational()+R200(l1a).exact_rational()))
        print("h1a = %la, h1b = %la, l1a = %la, l1b = %la, diff = %f\n" %(get_hex(h1a),get_hex(h1b),get_hex(l1a),get_hex(l1b),t8))
    
    
    # Add (h1,l1) by (hlog_alpha_i_m,llog_alpha_i_m)
    (h5a,l5a),(h5b,l5b) = TwoSum_m2(hlog_alpha_i_m,llog_alpha_i_m,h1a,l1a),Add22(hlog_alpha_i_m,llog_alpha_i_m,h1b,l1b)
    t9 = R200(((R200(h5b).exact_rational()+R200(l5b).exact_rational())-(R200(h5a).exact_rational()+R200(l5a).exact_rational()))/(R200(h5a).exact_rational()+R200(l5a).exact_rational()))
    print("h5a = %la, h5b = %la, l5a = %la, l5b = %la, diff = %f\n" %(get_hex(h5a),get_hex(h5b),get_hex(l5a),get_hex(l5b),t9))
    
    # log(2) =(h_log2,l_log2)
    h_log2 =RR('0x1.62e42fefa38p-1',16)
    l_log2 = RR('0x1.ef35793c7673p-45',16)
    
    #e*log(2)
    elog2_h,elog2_l = prod_dekker2(e,h_log2,l_log2)
   

    (hlog_fia,llog_fia),(hlog_fib,llog_fib) = TwoSum_m2(h5a,l5a,elog2_h,elog2_l),Add22(h5b,l5b,elog2_h,elog2_l)
    t10 = R200(((R200(hlog_fib).exact_rational()+R200(llog_fib).exact_rational())-(R200(hlog_fia).exact_rational()+R200(llog_fia).exact_rational()))/(R200(hlog_fia).exact_rational()+R200(llog_fia).exact_rational()))
    print("hlog_fia = %la, hlog_fib = %la, llog_fia = %la, llog_fib = %la, diff = %f" %(get_hex(hlog_fia),get_hex(hlog_fib),get_hex(llog_fia),get_hex(llog_fib),t10))
    #output: hlog_fi and llog_fi are 53-bit floating-point numbers.
    # hlog_fi: main value and llog_fi: error value
    
    return (RR(hlog_fia),RR(llog_fia)),(RR(hlog_fib),RR(llog_fib))



x =RR('0x1.b2ab940b3ed6p+8',16)
(hA1,lA1),(hA2,lA2) = test_two_fast(x,f)
A1 = RR(hA1).exact_rational()+RR(lA1).exact_rational()
A2 = RR(hA2).exact_rational()+RR(lA2).exact_rational()
B = log(x.exact_rational())
C1 = R200((A1-B)/B)
C2 = R200((A2-B)/B)
D1 = -log(abs(C1))/log(2.)
D2 = -log(abs(C2))/log(2.)
print((get_hex(x),D1,D2))
