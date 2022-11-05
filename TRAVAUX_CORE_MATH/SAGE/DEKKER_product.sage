load('fma.sage')

#x is a 53-bit floating-point numbers and n is the power of 2
def condition(x,n):
    #Output: True or False
    return x>=2^n and x<2^(n+1)


#x is a 53-bit floating-point numbers
def Split(x):
    # We use algorithm of Veltkampf split
    C = 1<<27
    z = x*C
    ah = z -(z-x)
    al = x - ah
    
    #output: ah and al are 53-bit floating-point numbers.
    # ah: main value and al: error value
    return RR(ah),RR(al)

#x is a 107-bit floating-point numbers
def Split2(x):
    # We use algorithm of Veltkampf split
    C = 1<<54
    z = x*C
    ah = z -(z-x)
    al = x - ah
    
    #output: ah and al are 107-bit floating-point numbers.
    # ah: main value and al: error value
    return R107(ah),R107(al)


################### Double Multiplication result double-double##########
#Input : a and b are 53-bit floating-point numbers
def prod_dekker(a,b):
    r1 = a*b 
    r2 = fma(a,b,-r1)
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)


#Input : a and b are 53-bit floating-point numbers
def Mul112(a,b):
    r1 = a*b 
    r2 = fma(a,b,-r1)
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)


def Mul112Cond(a,b):
    if abs(b)>abs(a):
        return Mul112(b,a)
    else:
        return Mul112(a,b)
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)

################### multiplication Double with double-double result double-double#######
#Input : a is 53-bit floating-point numbers
# bh: main value and  bl: error value
def prod_dekker2(a,bh,bl):
    ah,al =Split(a) # decomposition into 53-bits precision
    ''' We are turning b in bh,bl with bh: main value
                                    bl: error value'''
    
    # we use algorithm of DEKKER product
    r1 = a*(bh+bl)
    t1 = -r1 +(ah*bh)
    t2 = t1 +(ah*bl)
    t3 = t2 +(al*bh)
    r2 = t3 +(al*bl)
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)

#Input : a is 53-bit floating-point numbers
# bh: main value and  bl: error value
def Mul122(a,bh,bl):
    t1,t2 = Mul112(a,bh)
    t3 = a*bl
    t4 = t2+t3
    r1,r2 = Add112(t1,t4)
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)

################### Double-double multiplication ###################
#Input : ah and bh are main values
#  al and bl are error values
def prod_dekker3(ah,al,bh,bl):
    r1,r2a = prod_dekker(ah,bh)
    r2b = ah*bl
    r2c = bh*al
    r2d = al*bl
    r2 = r2a+r2b+r2c+r2d
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)

# Input: ah,al a double-double
# Output: ah.exact_rational()+al.exact_rational()
def exact2(ah,al):
   return ah.exact_rational() + al.exact_rational()

#Input : ah and bh are main values
#  al and bl are error values
def Mul222(ah,al,bh,bl):
    t1,t2 = Mul112(ah,bh)
    t3 = ah*bl
    t4 = bh*al
    t5,t6 = Add112Cond(t3,t4)
    r1,r2 = Add222Cond(t1,t2,t5,t6)
    # output: r1 and r2 are 53-bit floating-point numbers.
    # r1: main value and r2: error value
    return (r1,r2)




################### Double Addition  result double-double##########
#Input : a and b are 53-bit floating-point numbers
def FastTwoSum(a,b):
    s = a+b
    t = b -(s-a)
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return RR(s),RR(t)

#Input : a and b are 53-bit floating-point numbers
def TwoSum(a,b):
    s = a+b
    aprime = s-b
    bprime = s-aprime
    gammaa = a-aprime
    gammab = b-bprime
    t = gammaa+gammab
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return RR(s),RR(t)


#Input : a and b are 53-bit floating-point numbers
def Add112(a,b):
    s = a+b 
    z = s-a 
    t = b-z 
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t

#Input : a and b are 53-bit floating-point numbers
def Add112Cond(a,b):
    if abs(b)>abs(a):
        return Add112(b,a)
    else:
        return Add112(a,b)
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value

################### Addit Double with double-double result double-double#######
#Input : a is 53-bit floating-point numbers
# bh: main value and  bl: error value
def FastTwoSum_m(a,bh,bl):
    s,l1 =FastTwoSum(a,bh)
    t =l1+bl
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t


#Input : a is 53-bit floating-point numbers
# bh: main value and  bl: error value
def TwoSum_m(a,bh,bl):
    s,l1 =TwoSum(a,bh)
    t =l1+bl
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t



#Input : a is 53-bit floating-point numbers
# bh: main value and  bl: error value
def Add122(a,bh,bl):
    s,l = Add112(a,bh)
    t = l+bl
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t

def Add122Cond(a,bh,bl):
    s,l1 = Add12Cond(a,bh)
    t = l1+bl
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t


################### Double-double addition ##############################
#Input : ah and bh are main values
#  al and bl are error values
def FastTwoSum_m2(ah,al,bh,bl):
    s,l1 =FastTwoSum(ah,bh)
    t =l1+bl+al
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t

#Input : ah and bh are main values
#  al and bl are error values
def TwoSum_m2(ah,al,bh,bl):
    s,l1 =TwoSum(ah,bh)
    t =l1+bl+al
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t


#Input : ah and bh are main values
#  al and bl are error values
def Add222(ah,al,bh,bl):
    s,l = Add112(ah,bh)
    m = l+al
    t = m+bl
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t


#Input : ah and bh are main values
#  al and bl are error values
def Add222Cond(ah,al,bh,bl):
    s,l = Add112Cond(ah,bh)
    m = l+al
    t = m+bl
    # output: s and t are 53-bit floating-point numbers.
    # s: main value and t: error value
    return s,t

############################# Addition and Multiplication in one function##########
#Input : g are 107-bit floating-point numbers
#        h and h1 are main values
#        l and l1 are error values
def TwoSum_DEKKER(g,h1,l1,h,l):
    
    gh,gl = Split(g)
    # add gh,gl by (h1,l1) result (h2,l2)
    h2,l2 = TwoSum_m2(gh,gl,h1,l1)

    # multiply (h,l) by (h2,l2) result (h3,l3)
    h3,l3 = prod_dekker3(h2,l2,h,l)    

    # output: h3 and l3 are 53-bit floating-point numbers.
    # h3: main value and l3: error value
    return h3,l3


#Input : g are 107-bit floating-point numbers
#        h and h1 are main values
#        l and l1 are error values
def DEKKER_TwoSum(g,h1,l1,h,l):
    
    gh,gl = Split(g)
    # multiply gh,gl by (h1,l1) result (h2,l2)
    h2,l2 = prod_dekker3(h1,l1, gh,gl)

    # add (h,l) by (h2,l2) result (h3,l3)
    h3,l3 = TwoSum_m2(h2,l2,h,l)    

    # output: h3 and l3 are 53-bit floating-point numbers.
    # h3: main value and l3: error value
    return h3,l3


#Input : gh , h and g1h are main values
#        gl , l and g1l are error values
def MulAdd22(gh,gl,h,l,g1h,g1l):
    t1,t2 = Mul112(gh,h)
    t3,t4 = Add112(g1h,t1)
    t5 = gh*l
    t6 = gl*h 
    t7 = t2 + g1l
    t8 = t4 + t7
    t9 = t5 + t6 
    t10 = t8 + t9
    h3,l3 = Add112(t3,t10)
    # output: h3 and l3 are 53-bit floating-point numbers.
    # h3: main value and l3: error value
    return h3,l3


############################# triple-double numbers #####################
#Input : ah,am,al and bh,bm,bl are triple-double numbers
def Add333(ah,am,al,bh,bm,bl):
    rh,t1 = Add112Cond(ah,bh)
    t2,t3 =Add112Cond(am,bm)
    t7,t4 = Add112Cond(t1,t2)
    t6 = al+bl
    t5 = t3+t4
    t8 = t5+t6
    rm,rl = Add112Cond(t7,t8)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl

#Input : ah,am,al and bh,bm,bl are triple-double numbers
def Add33Cond(ah,am,al,bh,bm,bl):
    if abs(bh)>abs(ah):
        return Add333(bh,bm,bl,ah,am,al)
    else :
        return Add333(ah,am,al,bh,bm,bl)
 #output :rh,rm,rl is a triple-double numbers
     

#Input : ah,al  is a double-double numbers
#        bh,bm,bl is a triple-double numbers
def Add233(ah,al,bh,bm,bl):
    rh,t1 = Add112Cond(ah,bh)
    t2,t3 = Add112Cond(bm,al)
    t4,t5 = Add112Cond(t1,t2)
    t6 = t3+bl
    t7 = t6+t5
    rm,rl = Add112Cond(t4,t7)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl

#Input : ah,al  is a double-double numbers
#        bh,bm,bl is a triple-double numbers
def Add233change(ah,al,bh,bm,bl):
    rh,t1 = Add112(bh,ah)
    t2,t3 = Add112(bm,al)
    t4,t5 = Add112(t1,t2)
    t6 = t3+bl
    t7 = t6+t5
    rm,rl = Add112(t4,t7)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl

#Input : ah,al  is a double-double numbers
#        bh,bm,bl is a triple-double numbers
def Add233Cond(ah,al,bh,bm,bl):
    rh,t1 = Add12Cond(ah,bh)
    t2,t3 = Add12Cond(al,bm)
    t4,t5 = Add12Cond(t1,t2)
    t6 = t3+bl
    t7 = t6+t5
    rm,rl = Add12Cond(t4,t7)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl


#Input : a is a double numbers
#        bh,bm,bl is a triple-double numbers
def Add133(a,bh,bm,bl):
    rh,t1 = Add112Cond(a,bh)
    t2,t3 = Add112Cond(t1,bm)
    t4 = t3+bl
    rm,rl = Add112Cond(t2,t4)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl


#Input : a is a double numbers
#        bh,bm,bl is a triple-double numbers
def Add133Cond(a,bh,bm,bl):
    rh,t1 = Add12Cond(a,bh)
    t2,t3 = Add12Cond(t1,bm)
    t4 = t3+bl
    rm,rl = Add12Cond(t2,t4)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl


#Input : ah,am,al and bh,bm,bl are triple-double numbers
def Mul333(ah,am,al,bh,bm,bl):
    rh,t1 = Mul112Cond(ah,bh)
    t2,t3 = Mul112Cond(ah,bm)
    t4,t5 = Mul112Cond(bh,am)
    t6,t7 = Mul112Cond(am,bm)
    t8 = ah*bl
    t9 = al*bh
    t10 = am*bl
    t11 = al*bm 
    t12 = t8+t9
    t13 = t10+t11
    t14,t15 = Add112Cond(t1,t6)
    t16 = t7+t15
    t17 = t12+t13
    t18 = t16+t17
    t19,t20 = Add112Cond(t14,t18)
    t21,t22 = Add222Cond(t2,t3,t4,t5)
    rm,rl = Add222Cond(t21,t22,t19,t20)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl


#Input : ah,al and bh,bl are  double-double numbers
#        
def Mul223(ah,al,bh,bl):
    rh,t1 = Mul112(ah,bh)
    t2,t3 = Mul112(ah,bl)
    t4,t5 = Mul112(bh,al)
    t6 = al*bl
    t7,t8 = Add222Cond(t2,t3,t4,t5)
    t9,t10 = Add112Cond(t1,t6)
    rm,rl = Add222Cond(t7,t8,t9,t10)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl


#Input : ah,al  is a double-double numbers
#        bh,bm,bl is a triple-double numbers
def Mul233(ah,al,bh,bm,bl):
    rh,t1 = Mul112(ah,bh)
    t2,t3 = Mul112(ah,bm)
    t4,t5 = Mul112(ah,bl)
    t6,t7 = Mul112(bh,al)
    t8,t9 = Mul112(al,bm)
    t10 = al*bl
    t11,t12 = Add222Cond(t2,t3,t4,t5)
    t13,t14 = Add222Cond(t6,t7,t8,t9)
    t15,t16 = Add222Cond(t11,t12,t13,t14)
    t17,t18 = Add112Cond(t1,t10)
    rm,rl = Add222Cond(t17,t18,t15,t16)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl

#Input : a is a double numbers
#        bh,bl is a double-double numbers
def Mul123(a,bh,bl):
    rh,t1 = Mul112(a,bh)
    t2,t3 = Mul112(a,bl)
    t5,t4 = Add112Cond(t1,t2)
    t6 = t3+t4
    rm,rl = Add112Cond(t5,t6)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl  

#Input : a is a double numbers
#        bh,bm,bl is a triple-double numbers
def Mul133(a,bh,bm,bl):
    rh,t2 = Mul112(a,bh)
    t3,t4 = Mul112(a,bm)
    t5 = a*bl
    t9,t7 = Add112Cond(t2,t3)
    t8 = t4+t5
    t10 = t7+t8
    rm,rl = Add112Cond(t9,t10)
    #output :rh,rm,rl is a triple-double numbers
    return rh,rm,rl



