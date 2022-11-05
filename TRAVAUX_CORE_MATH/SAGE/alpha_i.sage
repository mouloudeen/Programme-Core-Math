load('decomposition.sage')
load('e.sage')

def calculation_of_alpha(k):
    L = [decomposition((2**k)/(2**k +i),53) for i in range(2**8)]
    L1 = [(get_hex(L[i][0]), get_hex(L[i][1])) for i in range(2**8)]
    return L1


def calculation_of_log_alpha_i(k):
    R200 = RealField(200)
    L = [decomposition(log(R200((2**k)/(2**k +i))),53) for i in range(256)]
    L1 = [(get_hex(L[i][0]), get_hex(L[i][1])) for i in range(2**8)]
    return L1