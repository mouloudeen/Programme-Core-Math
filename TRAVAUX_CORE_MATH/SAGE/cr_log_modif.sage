load('load.sage')
 
 
def cr_log_200(x,f):
    
    if  (x> RR('0x1.fedf47fc2216dp-1',16) and x<1) or (x> 1 and x <= RR('0x1.0178e5916f543p+0',16)) :
       
        r = x-1
        # we calculated len(f.list())
        F =f.list()
        F.reverse()
        lenf = len(F)

        # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
        
        
        # multiply (hr,lr) by (hf,lf) result (h,l)
        h1 = F[0]*r

        for i in range(1,lenf-1):
            

            # add (hfi,lfi) by (h1,l1) result (h2,l2)
            h2 = h1 + F[i]

            # multiply (h,l) by (h2,l2) result (h1,l1)
            h1 = h2 *r
        return R200(h1)
        
    #s represents the sign 
    #e  represents the exposant 
    #m represents fraction
    (s,m,e)=RR(x).sign_mantissa_exponent()
    #e=e-53 
    
    
    e = e+53
 
    # If x is a subnormal
    if s==1 and e==0 and m!=0:
        v = m
        e = e-1023
        v=v*2
        while v< 2^52:
            v*=2
            e=e-1
        m1 =v*2.^(-52)
        
    #If x is normal
    else:
        #print(s,e,m)
        m1 = m*2.^(-52)
        e =e-1
    #print("m1 = ",m1)
   
    M = R200(m1).exact_rational()
    binary = RR(m1).str(2)[2:10] # on récupère les 8 bits apres le 1 initial.
    i = int(binary,2)         # i est l'entier des 8 bits
    #print("i = ",i)
    alpha_i_m = RR(table_alpha_i[i][0],16).exact_rational()+ RR(table_alpha_i[i][1],16).exact_rational()#table calculé pour tous les i d'alpha_i_m
                                      # tel que log(alpha_i_m) soit de précision de 71 bits
    log_alpha_i_m= RR(table_log_alpha_i[i][0],16).exact_rational() + RR(table_log_alpha_i[i][1],16).exact_rational()
    r = M*alpha_i_m-1

    # we calculated len(f.list())
    F =f.list()
    F.reverse()
    lenf = len(F)

    # F[i] is the coefficient of each monomial :x^(len(F)-1-i)
    
    
    # multiply (hr,lr) by (hf,lf) result (h,l)
    h1 = F[0]*r

    for i in range(1,lenf-1):
        

        # add (hfi,lfi) by (h1,l1) result (h2,l2)
        h2 = h1 + F[i]

        # multiply (h,l) by (h2,l2) result (h1,l1)
        h1 = h2 *r
    A=R200(h1+log_alpha_i_m).exact_rational() 
    #print("A = ",A)
    h_log2 =RR('0x1.62e42fefa38p-1',16).exact_rational()
    l_log2 = RR('0x1.ef35793c7673p-45',16).exact_rational()
    B = A+ e*(h_log2+l_log2)
    #print("B = ",B)
    return R200(B)