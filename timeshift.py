from pylab import *
from scipy.interpolate import interp1d
from scipy.signal import correlate, argrelmax

N0 = 11
N1 = 1000
x0 = rand(N0)
t0 = linspace(0, 1, N0)
t1 = linspace(0, 1, N1)
x_interp = interp1d(t0, x0, 'linear')

x1 = x_interp(t1)
t2 = t1[int(round(N1*2/5)):int(round(N1*4/5))]
x2 = x1[int(round(N1*2/5)):int(round(N1*4/5))]
figure()
plot(t1, x1)
plot(t2, x2)
y = correlate(x1, x2, mode='valid')
figure()
plot(y)
figure()
plot(x1)
for n in argrelmax(y)[0]:
    plot(arange(0, len(x2))+n, x2, 'r')
