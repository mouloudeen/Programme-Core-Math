#Notre méthode avec l'algorithme de Tang et la méthode Gal avec 
# la fonction d'approximation
load('e.sage')
load('table_alpha_mo.sage') # table des valeurs alpha_i et log(alpha_i)
R200 = RealField(200)

f7 = (0.141270897304995757259860056365141645073890686035156).exact_rational()
f6 = (-0.16665893238270268472689394911867566406726837158203).exact_rational()
f5 = (0.19999998214690262177128943221759982407093048095703).exact_rational()
f4 = (-0.24999999998080402185962611838476732373237609863281).exact_rational()
f3 = (0.33333333333332582082420003644074313342571258544922).exact_rational()
f2 = (-0.5).exact_rational()
f1 = (1.0).exact_rational()

f(x)=f7 * x^7 + f6 * x^6 + f5 * x^5 + f4 * x^4 +f3 * x^3 + f2 * x^2 + f1 * x

k=8 # on prend k = 8, on va récupérer les 8 bits apres le 1 initial.
m=1+random() # On prend 1<m<2
M = R200(m).exact_rational()
binary = RR(m).str(2)[2:10] # on récupère les 8 bits apres le 1 initial.
i = int(binary,2)         # i est l'entier des 8 bits
alpha_i_m = R200(table_alpha_modified[i],16).exact_rational() #table calculé pour tous les i d'alpha_i_m
                                      # tel que log(alpha_i_m) soit de précision de 71 bits
log_alpha_i_m= R200(table_log_alpha_modified[i],16).exact_rational() 
A=R200(f(M*alpha_i_m-1)+log_alpha_i_m) #
B = R200(log(R200(M)))
print("m =",m)
print("i =",i)
print ("alpha_iprime = ",alpha_i_m )
print("log_alpha_i_prime :",log_alpha_i_m)
print("sol   =" ,A)                                                     
print("log(m)=",B)
C =R200((A.exact_rational()-B.exact_rational())/B.exact_rational())
print("erreur relative = ",C)
print("précision de %f bits" %(-log(C)/log(2.)))