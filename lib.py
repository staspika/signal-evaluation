#/usr/bin/env python2.7
# -*- encoding: utf-8-*-

from __future__ import division, unicode_literals, print_function
from pylab import *
from numpy import round
from scipy.io import loadmat
from fractions import gcd
import h5py

def plotBode(f, X):
    subplot(2,1,1)
    plot(f, abs(X), '_')
    xlim([0, 50])
    vlines(f, 0, abs(X))
    grid(True)
    subplot(2,1,2)
    plot(f, angle(X), '-')
    xlim([0, 50])
    #~ vlines(f, 0, abs(X))
    grid(True)

def lcm(a, b):
    return a*b/gcd(a, b)

def is_almost_int(x, d=1e-6):
    if abs(x - int(round(x))) < d:
        return True
    else:
        return False

def guess_rational(x, N=1000, d=1e-6):
    n = 0
    m_max = 128
    while True:
        n += 1
        m = n*x
        if True==is_almost_int(m, d):
            break
    return(n, int(round(m)))

def UT_ABB(fN=50/3, fT=1, dt=1e-6, t1=1, Ud0=0, Uq0=0, Uda=0, Udb=0, Uqa=0, Uqb=0):
    """
    Generate a signal according to ABB specification
    """
    omegaN = 2*pi*fN
    omegaT = 2*pi*fT
    NN = (1/fN)/dt
    NT = (1/fT)/dt
    t = arange(0, t1, dt) # Note that TT and TN must fit into signal a natural number of times for succesful dq0 transform.
    UT = (  Ud0*cos(omegaN*t) - Uq0*sin(omegaN*t)
          + Uda*cos(omegaN*t)*cos(omegaT*t)
          - Udb*cos(omegaN*t)*sin(omegaT*t)
          - Uqa*sin(omegaN*t)*cos(omegaT*t)
          + Uqb*sin(omegaN*t)*sin(omegaT*t))
    Ud = Uda + 1j*Udb
    Uq = Uqa + 1j*Uqb
    #Verify signal
    print('Settings:')
    print('|Ud| = {0}, ∠Ud = {1}°'.format(abs(Ud), angle(Ud)*180/pi))
    print('|Uq| = {0}, ∠Uq = {1}°'.format(abs(Uq), angle(Uq)*180/pi))
    return UT
