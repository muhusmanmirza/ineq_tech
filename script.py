#Python script

demo('ineqtech')
ineq = load(e = 'ineqtech', c = 'ineqtech.1')
rf = run(ineq)
print rf
ineq = load(e = 'ineqtech', c = 'ineqtech.2')
rb = run(ineq)
print rb
both = rf + rb
r = merge(both)
save(r, 'bif')
