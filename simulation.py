from matplotlib.pyplot import *
km1 = 1.0
km2 =10.0 
vmax1 = 1.0
vmax2 = 2.0
stoch = 6.0/5
c0 = 1.0
influx = 0.1
c = [c0]

for i in range(1,1000):
    f1 = c[-1]*vmax1/(km1+c[-1])
    f2 = c[-1]*vmax2/(km2+c[-1])
    c.append(c[-1]+stoch*f1-f2+influx)

plot(c)
show()

