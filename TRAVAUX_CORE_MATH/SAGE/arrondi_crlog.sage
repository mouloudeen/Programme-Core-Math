load('cr_log.sage')
load('load.sage')

o1 = open('logDN.wc','r')
r = o1.read()
WCDN = r.split('\n')
o1.close
o1h = open('hslogDN.wc','w')
o1m = open('mslogDN.wc','w')
o1l = open('lslogDN.wc','w')
for x in WCDN:
    h,m,l = cr_log_accurate_path_advanced(RR(x,16),U)
    o1h.write(get_hex(h))
    o1h.write('\n')
    o1m.write(get_hex(m))
    o1m.write('\n')
    o1l.write(get_hex(l))
    o1l.write('\n')
o1h.close
o1m.close
o1l.close
