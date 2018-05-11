#1 parameter bifurcation - cusps
ineq = load(e = 'ineqtechtotal', c = 'ineqtechtotal.1')
rf = run(ineq)
print rf
rb = run(ineq, DS = '-')
print rb
t_1db = rf + rb
print t_1db
save(t_1db, 't_1db')

#2 parameter bifurcation - cusps
lp1 = load(t_1db, c = 'ineqtechtotal_cusp', IRS = 2) #first cusp
cuspf = run(lp1)
print cuspf
cuspb = run(lp1, DS='-')
print cuspb
cusp = cuspf + cuspb
print cusp
lp2 = load(t_1db, c = 'ineqtechtotal_cusp', IRS = 8) #continuing thrid limit point
lpf = run(lp2)
print lpf
lpb = run(lp2, DS='-')
print lpb
lpt = lpf + lpb
t_2db = cusp + lpt
save(t_2db, 't_2db')
