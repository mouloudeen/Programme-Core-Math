o1 = open('logDN.wc','r')
r = o1.read()
WC = r.split('\n')
WC = WC[:-1]
o1.close

o1 = open('hlogDN.wc','r')
r = o1.read()
H6 = r.split('\n')
H6 = H6[:-2]
o1.close

o1 = open('mlogDN.wc','r')
r = o1.read()
M6 = r.split('\n')
M6 =M6[:-2]
o1.close

o1 = open('llogDN.wc','r')
r = o1.read()
L6 = r.split('\n')
L6 =L6[:-2]
o1.close

err = RR('0x1p-151',16)
print("left = h + (m + (l - err)) et right = h + (m + (l + err)) avec err = 0x1p-151\n")
print("Les cas où right!=left :\n")

for i in range(len(M6)):
    left = RR(H6[i],16) + (RR(M6[i],16) + (RR(L6[i],16) - err))
    right = RR(H6[i],16) + (RR(M6[i],16) + (RR(L6[i],16) + err))
    if left!=right:
        print(WC[i])
        print("\n")



