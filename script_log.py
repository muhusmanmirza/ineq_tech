#Python script

demo('logineqtech')
ineq = load(e = 'logineqtech', c = 'logineqtech.1')
rf = run(ineq)
print rf
#ineq = load(e = 'ineqtech', c = 'ineqtech.2')
rb = run(ineq, DS = '-')
print rb
both = rf + rb
r = merge(both)
save(r, 'biflog')
