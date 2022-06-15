#!/usr/bin/env python3.7

import gurobipy as gp

threshold = 187.1315205471496
# bin = [7002.400071806625, 2692.665193952964]
bin = [6786.228072149728, 2700.314222982642]
strain = [5637.38, 213.27, 2692.67, 643.94]
# bin = [643.9357959024691, 9051.12946985712]
# strain = [347.15, 639.18, 207.64, 5598.88, 2776.47]
# strain = [5637.38, 2692.67, 213.27, 643.94, 208.52]
strain_sum = sum(strain)
alpha = threshold
print("strain sum: ", strain_sum, "bin sum: ", sum(bin))
print("alpha: ", alpha)
m = gp.Model("BubbleSwap")
xs = []
obj = gp.LinExpr()
for b in bin:
    le = gp.LinExpr()
    x = m.addVars([i for i in range(len(strain))] ,vtype=gp.GRB.BINARY)
    xs.append(x)
    for i, v in x.items():
        le += v*strain[i]
    m.addConstr(le <= (b+alpha))
    obj += (b - le)**2

for i in range(len(strain)):
    le = gp.LinExpr()
    for item in xs:
        # print(type(item), item, item[1])
        le += item[i]
    m.addConstr(le == 1)

m.setObjective(obj, gp.GRB.MINIMIZE)
m.optimize()
m.printAttr('x')
print([int(var.x) for var in m.getVars()])
print(m.objVal)
print(m.display())