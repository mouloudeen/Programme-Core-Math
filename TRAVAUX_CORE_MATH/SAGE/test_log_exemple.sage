load('load.sage')
load('cr_log.sage')

x=RR('0x1.ffffffffffffep-1',16)
print("x = %la" %(get_hex(x)))
h,m,l =cr_log_accurate_path_advanced(x,U)
ml = m+l
print("h = %la, m = %la, l = %la, ml = %la\n" %(get_hex(h),get_hex(m),get_hex(l),get_hex(ml)))
a = log(x)
crlog_x= h+ml
print("cr_log(x) = %la, log(x) = %la\n" %(get_hex(crlog_x),get_hex(RR(a))))