load('load.sage')
m= RR(1.01709940124125)
binary = RR(m).str(2)[2:10] # on récupère les 8 bits apres le 1 initial.
i = int(binary,2) 
#print("i = ",i)
    
alpha_i_m = RR(table_alpha_modified[i],16).exact_rational() #table calculé pour tous les i d'alpha_i_m
                                      # tel que log(alpha_i_m) soit de précision de 71 bits
log_alpha_i_m= RR(table_log_alpha_modified[i],16).exact_rational() 
    
r = m * alpha_i_m - 1

F =f.list()
lenf = len(F)
    
h,l =prod_dekker(F[-1],r,53)
for i in range(1,lenf-1):
    h1,l1   = TwoSum_m(F[-1-i],h,l)
    h,l = prod_dekker2(r,h1,l1,53)
hfi,lfi = TwoSum_m(log_alpha_i_m,h,l)
A = RR(hfi).exact_rational()+RR(lfi).exact_rational()
b = log(RR(m))
B = R200(b).exact_rational()
C = R200((A-B)/B)
D=-log(abs(C))/log(2.)
print(D)