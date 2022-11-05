def nbitid(x): 
    X = x.exact_rational()
    i = 1 
    a =RealField(53)(log(X))
    b = RealField(53+i)(log(X))
    while a.exact_rational() == b.exact_rational(): 
        i=i+1 
        b = RealField(54+i)(log(X))
    return i