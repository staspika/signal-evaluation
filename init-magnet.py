#/usr/bin/env python2.7
# -*- encoding: utf-8-*-

from __future__ import division, unicode_literals, print_function
from pylab import *
from numpy import round
from scipy.io import loadmat
from fractions import gcd
import h5py

exec(open('./scripts/lib-tiara.py').read())
exec(open('./scripts/lib-dq0.py').read())
exec(open('./scripts/lib.py').read())

mu0 = 4e-7*pi
mu = mu0
D = 6
h = 6
r = 11
I = 83 # Ampere
I1 = I
I2 = I

phi1 = arctan((D-h)/r)
r1 = sqrt((D-h)**2 + r**2)
B1=mu*I1/(2*pi*r1)

phi2 = arctan(h/r)
r2 = sqrt(h**2 + r**2)
B2=mu*I2/(2*pi*r2)

Bx = B1*cos(pi/2-phi1) + B2*cos(pi/2-phi2)
By = B1*sin(pi/2-phi1) - B2*sin(pi/2-phi2)

Bt=sqrt(Bx**2 + By**2)
print(Bt)
