load('load.sage')

#Test for Split
print("test for Split\n")
x1 =RR('0x1.906f0ca2f82fbp+2',16)
h,l = Split(x1)
#print(" x = %la , h = %la , l = %la\n" %(get_hex(x1),get_hex(h),get_hex(l)))
assert x1.exact_rational() == h.exact_rational()+l.exact_rational() ,"Split"



#test for prod_dekker 
print("Test for prod_dekker\n")
a =RR('0x1.421f273fc712bp+9',16)
b =RR('0x1.e707af8e689a7p+9',16)
h,l = prod_dekker(a,b)
#print("a = %la, b = %la, h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(h),get_hex(l))) 
assert a.exact_rational()*b.exact_rational() == h.exact_rational()+l.exact_rational() ,"prod_dekker"


#test for Mul112
print("Test for Mul112\n")
a =RR('0x1.421f273fc712bp+9',16)
b =RR('0x1.e707af8e689a7p+9',16)
h,l = Mul112(a,b)
#print("a = %la, b = %la, h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(h),get_hex(l))) 
assert a.exact_rational()*b.exact_rational() == h.exact_rational()+l.exact_rational() ,"Mul112"


#test for prod_dekker2
print("Test for prod_dekker2\n")
a =RR('0x1.c2148ef302f3fp+9',16)
b =RR('0x1.b9a87a94da0c1p+9',16)
bh,bl = Split(b)
h,l = prod_dekker2(a,bh,bl)
#print("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert a.exact_rational()*(bh.exact_rational()+bl.exact_rational()) == h.exact_rational()+l.exact_rational(),"prod_dekker2"

# test for Mul122
print("Test for Mul122\n")
a =RR('0x1.c2148ef302f3fp+9',16)
b =RR('0x1.b9a87a94da0c1p+9',16)
bh,bl = Split(b)
h,l = Mul122(a,bh,bl)
#print("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert a.exact_rational()*(bh.exact_rational()+bl.exact_rational()) -(h.exact_rational()+l.exact_rational()) < l.ulp() ," Mul122"

#test for prod_dekker3 
print(" Test for prod_dekker3\n");
a =RR('0x1.39ad06459895p+9',16)
ah,al = Split(a)
b =RR('0x1.5bc29786125eap+7',16)
bh,bl = Split(b)
h,l = prod_dekker3(ah,al,bh,bl)
#print("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert (ah.exact_rational()+ al.exact_rational())*(bh.exact_rational()+bl.exact_rational()) - (h.exact_rational()+l.exact_rational()) < l.ulp() ,"prod_dekker3"

#test for Mul222
print("Test for Mul222\n")
a =RR('0x1.39ad06459895p+9',16)
ah,al = Split(a)
b =RR('0x1.5bc29786125eap+7',16)
bh,bl = Split(b)
h,l = Mul222(ah,al,bh,bl)
#print("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert (ah.exact_rational()+ al.exact_rational())*(bh.exact_rational()+bl.exact_rational()) - (h.exact_rational()+l.exact_rational()) < l.ulp() ,"Mul222"


#test for Mul133
print("Test for Mul133\n")
a = RR('0x1.5b2057ad879acp+8',16)
bh = RR('0x1.c0b5b34724d3p+17',16)
bm = RR('-0x1.18440552f325bp-12',16)
bl = RR('0x1.f109p-66',16)
h,m,l = Mul133(a,bh,bm,bl)
#print("a = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n " %(get_hex(a),get_hex(bh),get_hex(bm),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert a.exact_rational()* (bh.exact_rational()+bm.exact_rational()+bl.exact_rational()) - (h.exact_rational()+m.exact_rational() +l.exact_rational()) < l.ulp() ,"Mul133"

#test for Mul123
print("Test for Mul123\n")
a = RR('0x1.158922b181abfp+11',16)
bh = RR('0x1.65580ca8ecfe5p+8',16)
bl = RR('0x1.65580ca8ecfe5p-46',16)
h,m,l = Mul123(a,bh,bl)
#print("a = %la, bh = %la, bl = %la , h = %la, m = %la, l = %la\n " %(get_hex(a),get_hex(bh),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert abs(a.exact_rational()*(bh.exact_rational()+bl.exact_rational()) - (h.exact_rational()+m.exact_rational()+l.exact_rational()))< l.ulp() ,"Mul123"

#test for Mul223
print("Test for Mul223\n")
a = RR('0x1.d9bb8bb2434a7p+8',16)
ah,al = Split(a)
b = RR('0x1.969146a9b9d69p+8',16)
bh,bl = Split(b)
h,m,l = Mul223(ah,al,bh,bl)
#print("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, m = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert (ah.exact_rational()+ al.exact_rational())*(bh.exact_rational()+bl.exact_rational()) - (h.exact_rational()+m.exact_rational()+l.exact_rational()) < l.ulp() ,"Mul223"

#test for Mul233
print("Test for Mul233\n")
a = RR('0x1.e716676d322c6p+9',16)
ah,al = Split(a)
bh = RR('0x1.81260f7a7c404p+18',16)
bm = RR('0x1.db484eca53ea2p-9',16)
bl = RR('-0x1.73c54ep-63',16)
h,m,l = Mul233(ah,al,bh,bm,bl)
#print("a = %la, ah = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n " %(get_hex(a),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bm),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert (ah.exact_rational()+al.exact_rational())* (bh.exact_rational()+bm.exact_rational()+bl.exact_rational()) - (h.exact_rational()+m.exact_rational() +l.exact_rational()) < l.ulp() ,"Mul233"

#test for Mul333
print("Test for Mul333\n")
ah = RR('0x1.4b7b1763d5c36p+19',16)
am = RR('-0x1.2eac725f98608p-7',16)
al = RR('0x1.d58524p-61',16)
bh = RR('0x1.6906123c043bp+19',16)
bm = RR('-0x1.1b638b9cda375p-8',16)
bl = RR('-0x1.eae1fp-62',16)
h,m,l = Mul333(ah,am,al,bh,bm,bl)
#print("ah = %la, am = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n " %(get_hex(ah),get_hex(am),get_hex(al),get_hex(bh),get_hex(bm),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert (ah.exact_rational()+am.exact_rational()+al.exact_rational()) * (bh.exact_rational()+bm.exact_rational()+bl.exact_rational()) -(h.exact_rational()+m.exact_rational() +l.exact_rational()) < l.ulp() ,"Mul333"


#test for FastTwoSum
print("Test for FastTwoSum\n")
a =RR('0x1.09e5334316b0ap+9',16)
b =RR('0x1.6afe66d00be08p+8',16)
h,l = FastTwoSum(a,b)
#print("a = %la, b = %la, h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(h),get_hex(l))) 
assert a.exact_rational()+b.exact_rational() == h.exact_rational()+l.exact_rational() ,"FastTwoSum"

#test for TwoSum
print("Test for TwoSum\n")
a =RR('0x1.09e5334316b0ap+9',16)
b =RR('0x1.6afe66d00be08p+8',16)
h,l = TwoSum(a,b)
#print("a = %la, b = %la, h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(h),get_hex(l))) 
assert a.exact_rational()+b.exact_rational() == h.exact_rational()+l.exact_rational() ,"TwoSum"


#test for Add112
print("Test for Add112\n")
a =RR('0x1.09e5334316b0ap+9',16)
b =RR('0x1.6afe66d00be08p+8',16)
h,l = Add112(a,b)
#print("a = %la, b = %la, h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(h),get_hex(l))) 
assert a.exact_rational()+b.exact_rational() == h.exact_rational()+l.exact_rational() ,"Add112"



#test for FastTwoSum_m
print("Test for FastTwoSum_m\n")
a =RR('0x1.01f03e5b06beep+9',16)
b =RR('0x1.3f20ef1cf5991p+7',16)
bh,bl = Split(b)
h,l = FastTwoSum_m(a,bh,bl)
#print("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert a.exact_rational()+bh.exact_rational()+bl.exact_rational() == h.exact_rational()+l.exact_rational() ,"FastTwoSum_m"

#test for TwoSum_m
print("Test for TwoSum_m\n")
a =RR('0x1.01f03e5b06beep+9',16)
b =RR('0x1.3f20ef1cf5991p+7',16)
bh,bl = Split(b)
h,l = TwoSum_m(a,bh,bl)
#print("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert a.exact_rational()+bh.exact_rational()+bl.exact_rational() == h.exact_rational()+l.exact_rational() ,"TwoSum_m"


#test for Add122
print("Test for Add122\n")
a =RR('0x1.01f03e5b06beep+9',16)
b =RR('0x1.3f20ef1cf5991p+7',16)
bh,bl = Split(b)
h,l = Add122(a,bh,bl)
#print("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert a.exact_rational()+bh.exact_rational()+bl.exact_rational() == h.exact_rational()+l.exact_rational() ,"Add122"


#test for FastTwoSum_m2
print("Test for FastTwoSum_m2\n")
a =RR('0x1.1b33e8b39d9f1p+4',16)
ah,al = Split(a)
b =RR('0x1.e40fc62966a64p+9',16)
bh,bl = Split(b)
h,l = FastTwoSum_m2(ah,al,bh,bl)
#print("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert ah.exact_rational()+al.exact_rational()+bh.exact_rational()+bl.exact_rational() == h.exact_rational()+l.exact_rational(),"FastTwoSum_m2"

#test for TwoSum_m2
print("Test for TwoSum_m2\n")
a =RR('0x1.1b33e8b39d9f1p+4',16)
ah,al = Split(a)
b =RR('0x1.e40fc62966a64p+9',16)
bh,bl = Split(b)
h,l = TwoSum_m2(ah,al,bh,bl)
#print("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert ah.exact_rational()+al.exact_rational()+bh.exact_rational()+bl.exact_rational() == h.exact_rational()+l.exact_rational(),"TwoSum_m2"


#test for Ad222
print("Test for Add222\n")
a =RR('0x1.1b33e8b39d9f1p+4',16)
ah,al = Split(a)
b =RR('0x1.e40fc62966a64p+9',16)
bh,bl = Split(b)
h,l = Add222(ah,al,bh,bl)
#print("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " %(get_hex(a),get_hex(b),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bl),get_hex(h),get_hex(l)))
assert ah.exact_rational()+al.exact_rational()+bh.exact_rational()+bl.exact_rational() == h.exact_rational()+l.exact_rational(),"Add222"



#test Add133 
print("Test for Add133\n")
a = RR('0x1.ce7079d3ea68p+7',16)
bh = RR('0x1.23b1a4cf5f82p+18',16)
bm = RR('0x1.80d37956f706fp-9',16)
bl = RR('-0x1.2ccap-67',16)
h,m,l = Add133(a,bh,bm,bl)
#print("a = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n " %(get_hex(a),get_hex(bh),get_hex(bm),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert a.exact_rational()+ (bh.exact_rational()+bm.exact_rational()+bl.exact_rational()) -(h.exact_rational()+m.exact_rational() +l.exact_rational()) < l.ulp() ,"Add133"

#test for Add233
print("Test for Add233\n")
a = RR('0x1.20a247f9999edp+9',16)
ah,al = Split(a)
bh = RR('0x1.0dd04f65b5954p+18',16)
bm = RR('-0x1.9b0c2cc926733p-9',16)
bl = RR('-0x1.46ccd4p-65',16)
h,m,l = Add233(ah,al,bh,bm,bl)
#print("a = %la, ah = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n " %(get_hex(a),get_hex(ah),get_hex(al),get_hex(bh),get_hex(bm),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert (ah.exact_rational()+al.exact_rational())+(bh.exact_rational()+bm.exact_rational()+bl.exact_rational()) - (h.exact_rational()+m.exact_rational() +l.exact_rational()) < l.ulp() ,"Add233"


#test for Add333
print("Test for Add333\n")
ah = RR('0x1.2937ce70794fp+17',16)
am = RR('0x1.2ffbef6672a88p-11',16)
al = RR('0x1.6bd29cp-65',16)
bh = RR('0x1.be1e53af7546p+17',16)
bm = RR('0x1.031eb811ecc19p-9',16)
bl = RR('-0x1.16769p-64',16)
h,l,m = Add333(ah,am,al,bh,bm,bl)
#print("ah = %la, am = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n " %(get_hex(ah),get_hex(am),get_hex(al),get_hex(bh),get_hex(bm),get_hex(bl),get_hex(h),get_hex(m),get_hex(l)))
assert (ah.exact_rational()+am.exact_rational()+al.exact_rational()) + (bh.exact_rational()+bm.exact_rational()+bl.exact_rational()) -(h.exact_rational()+m.exact_rational() +l.exact_rational()) < l.ulp() ,"Add333"


# test for TwoSum_DEKKER
print("Test for TwoSum_DEKKER\n")
a =RR('0x1.1b33e8b39d9f1p+4',16)
h1,l1 = Split(a)
b = RR('0x1.acfd611270449p+6',16)
h,l = Split(b)
g10h,g10l = Split(g10)
x,y = TwoSum_DEKKER(g10,h1,l1,h,l)
#print("a = %la, b = %la, h1= %la, l1 = %la, h = %la, l = %la , x = %la, y = %la\n " %(get_hex(a),get_hex(b),get_hex(h1),get_hex(l1),get_hex(h),get_hex(l),get_hex(x),get_hex(y)))
assert ((h1.exact_rational()+l1.exact_rational()+g10h.exact_rational()+g10l.exact_rational())*(h.exact_rational()+l.exact_rational())) -(x.exact_rational()+y.exact_rational()) < y.ulp(), "TwoSum_DEKKER"

# test for DEKKER_TwoSum
print("test for DEKKER_TwoSum\n")
a = RR('0x1.cd5216a848233p+9',16)
h1,l1 = Split(a)
b = RR('0x1.7289960c1bc21p+8',16)
h,l = Split(b)
g10h,g10l = Split(g10)
x,y = DEKKER_TwoSum(g10,h1,l1,h,l)
#print("a = %la, b = %la, h1= %la, l1 = %la, h = %la, l = %la , x = %la, y = %la\n " %(get_hex(a),get_hex(b),get_hex(h1),get_hex(l1),get_hex(h),get_hex(l),get_hex(x),get_hex(y)))
assert ((h1.exact_rational()+l1.exact_rational())*(g10h.exact_rational()+g10l.exact_rational())+(h.exact_rational()+l.exact_rational())) -(x.exact_rational()+y.exact_rational()) < y.ulp(), "DEKKER_TwoSum"

# test for MulAdd22
print("test for MulAdd22\n")
a = RR('0x1.cd5216a848233p+9',16)
h1,l1 = Split(a)
g10h,g10l = Split(g10)
g9h,g9l = Split(g9)
x,y = MulAdd22(g10h,g10l,h1,l1,g9h,g9l)
#print("a = %la,  h1= %la, l1 = %la, x = %la, y = %la\n " %(get_hex(a),get_hex(h1),get_hex(l1),get_hex(x),get_hex(y)))
assert ((h1.exact_rational()+l1.exact_rational())*(g10h.exact_rational()+g10l.exact_rational())+(g9h.exact_rational()+g9l.exact_rational())) -(x.exact_rational()+y.exact_rational()) < y.ulp(), "MulAdd22"