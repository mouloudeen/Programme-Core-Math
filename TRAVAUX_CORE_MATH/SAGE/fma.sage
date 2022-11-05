# return x*y+z with only one rounding
# assumes x, y and z have the same parent (same precision and same rounding)
# x=RR(pi); y=RR(log(2)); z=RR(exp(1))
# t=fma(x,y,z) # '0x1.3955e66527c1p+2'
def fma(x,y,z):
   R = x.parent()
   assert y.parent() == R
   assert z.parent() == R
   X = x.exact_rational()
   Y = y.exact_rational()
   Z = z.exact_rational()
   return R(X*Y+Z)
