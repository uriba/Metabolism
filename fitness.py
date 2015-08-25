from sympy import *
from sympy.interactive import printing
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

g,x,c,E_i,k_i,E_a,k_a,E_b,k_b = symbols('g,x c E_i k_i E_a k_a E_b k_b')

growth_rate = Eq(g,(c*E_b/(c+k_b)))

fixed_resources = Eq(E_b+E_a+E_i,1)

steady_state = Eq(c*E_a/(c+k_a)+x*E_i/(x+k_i)-c*E_b/(c+k_b))

x = solve([growth_rate,steady_state],[x,c])[0][0] # we only need x, but with substitution for c
x = simplify(x.subs(E_i,1-E_a-E_b))
c_max = 10.0 #maximal allowed concentration of c
g_val = 0.25
k_a_val = 20

kis = np.linspace(0.01,2)
kbs = np.linspace(0.01,2)
res = []

for k_b_val in kbs:
    # check if can be grown autotrophically:
    E_b_min = solve(growth_rate.subs([(c,c_max),(g,g_val),(k_b,k_b_val),(k_a,k_a_val)]),E_b)[0]
    if E_b_min >=1:
        continue
    for k_i_val in kis:
        if E_b_min < 1:
            if c_max*(1-E_b_min)/(c_max+k_a_val) >= g_val:
                res.append((k_i_val,k_b_val,0))
            E_i_max = 1-E_b_min
            if E_i_max > g_val:
                X = g_val*k_i_val/(E_i_max-g_val)
                res.append((k_i_val,k_b_val,X))

        # find the minimal 
fig = plt.figure()
ax = fig.add_subplot(111,projection = '3d')
xs,ys,zs = zip(*res)
ax.scatter(xs,ys,zs)
plt.draw()
#ax.set_zlim3d(top=5)           
