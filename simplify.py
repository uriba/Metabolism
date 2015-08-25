from sympy import *
c1,c2,fi,ka,va,kb,vb,kc,vc,kb2,vb2 = symbols('c1 c2 fi ka va kb vb kc vc kb2 vb2 ')
expr1 = Eq(0,fi+c2*vc/(c2+kc)-c1*va/(c1+ka)-c1*vb/(c1+kb))
expr2 = Eq(0,-c2*vc/(c2+kc)+2*c1*va/(c1+ka)-c2*vb2/(c2+kb2))
solve([expr1,expr2],[c1,c2])


