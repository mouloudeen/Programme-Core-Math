# print like in C: 0x1.xxxp+e
def get_hex(x):
   if x < 0:
      return '-' + get_hex(-x)
   s = x.hex()
   if s[2] in ['0','1']:
      e = 0
   elif s[2] in ['2','3']:
      s = (x/2).hex()
      e = 1
   elif s[2] in ['4','5','6','7']:
      s = (x/4).hex()
      e = 2
   else:
      s = (x/8).hex()
      e = 3
   # add e to exponent
   s = s.split('p')
   e = ZZ(s[1])+e
   if e >= 0:
      s[1] = '+' + e.str()
   else:
      s[1] = e.str()
   return s[0] + 'p' + s[1]
