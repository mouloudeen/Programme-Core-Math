load('cr_log.sage')
 
o1 = open('../src/log.wc','r')
r = o1.read()
WC = r.split('\n')
WC = WC[:-1]
n=0

L1 = []
L2 = []
L3=[]
L4 = []
L6 = []
o1.close
print("Reading stage finished")
i=0
for x in WC:
    
    x =RR(x,16)
    hA,lA = cr_log_accurate_path(x,G)
    A = RR(hA).exact_rational()+RR(lA).exact_rational()
    B = log(x.exact_rational())
    C = R200((B-A)/B)
    D = R200((-log(abs(C)))/log(2.0))
    L1.append(D)
    L2.append(get_hex(x))
    if D<101:
        #print("x = %la ,D = %e ,n = %d, i = %d \n" %(get_hex(x),D,n,i))
        n+=1
        L3.append(x)
        L6.append(D)
    i=i+1 
dico ={}
for i in range(len(WC)):
    dico[L1[i]]=L2[i]
L1.sort()  
n1 = -1073
while n1<1015:
    j = sum(condition(x,n1) for x in L3)
    L4.append((j,n1,n1+1))
    n1=n1+1
L5 = [L4[i] for i in range(len(L4)) if L4[i][0]!=0]  

print("Comparing stage finished")     
o2 = open('resultat.txt','w')
for i in range(len(L5)):
    o2.write("(%d, entre 2^(%d) et 2^(%d))\n" %(L5[i][0],L5[i][1],L5[i][2]))
o2.write("\n")
o2.write("\n")
o2.write("nombre d'erreur en dessous de 101 bits de précisions est %d\n" %(n))
min_L = L1[0]
max_L =L1[-1]
o2.write("le minimum de précisions est %d pour %la\n" %(min_L,dico[min_L]))
o2.write("le maximum de précision est %d pour %la\n" %(max_L,dico[max_L]))
o2.write("Les valeurs avec leurs précisons < 101:\n")
L6.sort()
if len(L6)<100:
    for i in range(len(L6)):
        o2.write("%d pour %la\n" %(L6[i],dico[L6[i]]))
o2.close
