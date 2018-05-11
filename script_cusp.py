#Python script
demo('ineqtech')

#1 parameter bifurcation - folds
ineq = load(e = 'ineqtech', c = 'ineqtech.1')
rf = run(ineq)
print rf
ineq = load(e = 'ineqtech', c = 'ineqtech.2')
rb = run(ineq)
print rb
onedbif = rf + rb
print onedbif
save(onedbif, 'onedbif')

#2 parameter bifurcation - cusps
lp1 = load(onedbif, c = 'ineqtech_cusp1', IRS = 2) #first cusp
cusp1f = run(lp1)
print cusp1f
cusp1b = run(lp1, DS='-')
print cusp1b
cusp1 = cusp1f + cusp1b
print cusp1
save(cusp1, 'cusp1')
lp2 = load(onedbif, c = 'ineqtech_cusp1', IRS = 16) #second cusp
cusp2f = run(lp2)
print cusp2f
cusp2b = run(lp2, DS='-')
print cusp2b
cusp2 = cusp2f + cusp2b
save(cusp2, 'cusp2')
twodbif = cusp1 + cusp2
lp3 = load(onedbif, c = 'ineqtech_cusp1', IRS = 10) #continuing thrid limit point
lp3f = run(lp3)
print lp3f
lp3b = run(lp3, DS='-')
print lp3b
lp3 = lp3f + lp3b
save(lp3, 'lp3')
twodbif = twodbif + lp3
po1 = load(onedbif, c = 'limitcycle.1', IRS = 9) #continuing thrid limit point
po = run(po1)
print po
save(po, 'po')
twodbif = twodbif + po
merge(twodbif)
save(twodbif,'twodbif')
