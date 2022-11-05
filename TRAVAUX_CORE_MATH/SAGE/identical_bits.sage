def identical_bits_aux(s,prec):
    s = s.replace('.','')
    # skip minus sign if any
    j = 0
    if s[0]=='-':
        j += 1
    # skip leading zeros
    while s[j]=='0':
        j += 1
    assert s[j] == '1', "s[j] == '1'" # leading bit
    # skip next prec+1 bits
    for i in range(prec+1):
        assert s[j] in ['0','1'], "s[j] in ['0','1']"
        j += 1
    k = j + 1
    while k < len(s) and s[j] == s[k]:
        k += 1
    assert k < len(s), "k < len(s)"
    # bits s[j]...s[k-1] are identical
    k -= j
    return k



# return the number of identical bits after the rounding bit
# for function f and argument x (of type RealField)
def identical_bits(f,x,rnd='all'):
    assert rnd in ['all','directed']
    prec = x.prec()
    X = x.exact_rational()
    large_prec = 4*prec
    s = n(f(X),large_prec).str(2)
    if rnd=='all':
        return identical_bits_aux(s,prec)
    else:
        return identical_bits_aux(s,prec-1)