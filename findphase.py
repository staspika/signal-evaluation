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

# Read data record
hdf5_file = h5py.File('./data/test04_20140627_131439.tdms.hdf5', 'r')
x = hdf5_file['data']['u15']
dt = x.attrs['wf_increment']
fN = 50/3
fT = 10*fN/128
s = Signal(x.value, dt)
#~ # Mock data record
#~ fN = 50/3
#~ fT = 10*fN/128
#~ dt = 1e-5
#~ t1 = 56.22
#~ Ud0 = 15750
#~ Uq0 = 0
#~ Uda = 0
#~ Udb = 150
#~ Uqa = 0
#~ Uqb = 0
#~ s = Signal(UT_ABB(fN=fN, fT=fT, dt=dt, t1=t1, Ud0=Ud0, Uq0=Uq0, Uda=Uda, Udb=Udb, Uqa=Uqa, Uqb=Uqb), dt)
#~ x = loadmat('./data/chunks/chunk16.mat')
#~ dt = mean(diff(x['t']))
#~ fN = float(x['fN'])
#~ fT = float(x['fT'])
#~ s = Signal(x['u15'].flatten(), dt)

# Trim signal
(a, b) = guess_rational(fN/fT, d=1e-2)
TN = 1/fN
K = int(floor(len(s.x)*dt/(b*TN))) # Times minimal sample in record
M = int(round(K*b*TN/dt))
s = Signal(s.x[:M], dt)
N = int(rand()*len(s.x))
#~ s.x = append(s.x[:N], s.x[N:])

#~ # Reduce phase shift
#~ x1_test = cos(2*pi*fT*arange(0, 1/fN, dt) + 0)
#~ x2_test = cos(2*pi*fT*arange(0, 1/fT, dt) + 0)
#~ y1 = correlate(s.x, x1_test)
#~ a1 = argrelmax(y1)[0]
#~ y2 = correlate(s.x, x2_test)
#~ a2 = argrelmax(y2)[0]
#~ for i in range(0, len(a1)):
    #~ for j in range(0, len(a2)):
        #~ if abs(a1[i]-a2[j])*dt < 0.10/fN:
            #~ print(a1[i], a2[j])
#~ shift1=a1[0][0]
#~ shift2=a2[0][0]
#~ phase0 = 2*pi*shift*dt*fT
#~ print('Phase0 of disturbance signal is {0}°'.format(phase0*180/pi))
#~ shift = 10358
#~ s.x = append(s.x[shift:], s.x[:shift])

# Find Ud and Uq
u_d, u_q, u_0 = park(s.x, s.time_track(), fN)
Ud = temp_fft(u_d, dt, fT)
Uq = temp_fft(u_q, dt, fT)
print('Calculations:')
print('|Ud| = {0}, ∠Ud = {1}°'.format(abs(Ud), angle(Ud)*180/pi))
print('|Uq| = {0}, ∠Uq = {1}°'.format(abs(Uq), angle(Uq)*180/pi))

# Plot Fourier transform of u_d
f = rfftfreq(len(u_d), dt)
X = (2/len(u_d))*rfft(u_d)
fig1 = figure()
plotBode(f, X)
n = argmax(abs(X[1:]))
xlim([0, 50])
# Plot Fourier transform of u_q
f = rfftfreq(len(u_q), dt)
X = (2/len(u_q))*rfft(u_q)
fig1 = figure()
plotBode(f, X)
n = argmax(abs(X[1:]))
xlim([0, 50])
