load('load.sage')
load('cr_log.sage')

o1 = open('logDN.wc','r')
r = o1.read()
WC = r.split('\n')
WC = WC[:-1]
o1.close

hfi=[]
mfi = []
lfi = []
mlfi=[]
for x in WC:
    h,m,l=cr_log_accurate_path_advanced(RR(x,16),U)
    ml = m+l
    hfi.append(get_hex(h))
    mfi.append(get_hex(m))
    lfi.append(get_hex(l))
    mlfi.append(get_hex(ml))

o1 = open('hlogDN.wc','r')
r = o1.read()
H6 = r.split('\n')
H6 = H6[:-2]
o1.close

n=0
for i in range(len(H6)):
    if RR(hfi[i],16)!=RR(H6[i],16):
        n+=1
        print(i)
print("%d different h entre C et Sage" %(n))

o1 = open('mlogDN.wc','r')
r = o1.read()
M6 = r.split('\n')
M6 =M6[:-2]
o1.close

n=0
for i in range(len(M6)):
    if RR(mfi[i],16)!=RR(M6[i],16):
        n+=1
        print(WC[i])
print("%d different m entre C et Sage" %(n))

o1 = open('llogDN.wc','r')
r = o1.read()
L6 = r.split('\n')
L6 =L6[:-2]
o1.close

n=0
for i in range(len(L6)):
    if RR(lfi[i],16)!=RR(L6[i],16):
        n+=1
        print(WC[i])
print("%d different l entre C et Sage" %(n))

o1 = open('mllogDN.wc','r')
r = o1.read()
ML6 = r.split('\n')
ML6 = ML6[:-1]
o1.close

n=0
for i in range(len(M6)):
    if RR(mlfi[i],16)!=RR(ML6[i],16):
        n+=1
        print(WC[i])

print("%d different ml entre C et Sage" %(n))

