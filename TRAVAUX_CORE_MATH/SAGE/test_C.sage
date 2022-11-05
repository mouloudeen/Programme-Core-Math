load('cr_log.sage')

oh = open('hlogDN.wc','r')
r = oh.read()
h = r.split('\n')
h = h[:-2]
om = open('mlogDN.wc','r')
r = om.read()
m = r.split('\n')
m = m[:-2]
ol = open('llogDN.wc','r')
r = ol.read()
l = r.split('\n')
l = l[:-2]
oh.close
om.close
ol.close


o1h = open('hslogDN.wc','r')
r = o1h.read()
hs = r.split('\n')
hs = hs[:-1]
o1m = open('mslogDN.wc','r')
r = o1m.read()
ms = r.split('\n')
ms = ms[:-1]
o1l = open('lslogDN.wc','r')
r = o1l.read()
ls = r.split('\n')
ls = ls[:-1]
o1h.close
o1m.close
o1l.close

print("Test rounding mode RNDN")
j=0
for i in range(35072):
    if RR(h[i],16) == RR(hs[i],16) and RR(m[i],16) == RR(ms[i],16) and RR(l[i],16) == RR(ls[i],16):
        j+=1
        print("j = %d" %(j))
    else :
        print("i = %d\n" %(i))
    
print("%d on 35073 cases" %(j))
