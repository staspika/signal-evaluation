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


#~ # Read data record
#~ x = loadmat('./data/chunks/chunk16.mat')
#~ dt = mean(diff(x['t']))
#~ fN = float(x['fN'])
#~ fT = float(x['fT'])
#~ s = Signal(x['u15'].flatten(), dt)
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
#~ dt = 1e-4
#~ t1 = 56.22
#~ Ud0 = 15750
#~ Uq0 = 0
#~ Uda = 0
#~ Udb = 150
#~ Uqa = 0
#~ Uqb = 0
#~ s = Signal(UT_ABB(fN=fN, fT=fT, dt=dt, t1=t1, Ud0=Ud0, Udb=Udb), dt)
# Plot data record before trimming
fig1 = figure()
plot(s.time_track(), s.x)
f = rfftfreq(len(s.x), dt)
X = (2/len(s.x))*rfft(s.x)
fig2 = figure()
plotBode(f, X)
xlim([0, 50])
suptitle('Uncut signal record')
# Trim data record
(a, b) = guess_rational(fN/fT, d=1e-2)
TN = 1/fN
K = int(floor(len(s.x)*dt/(b*TN))) # Times minimal sample in record
M = int(round(K*b*TN/dt))
s = Signal(s.x[:M], dt)
# Test of lib-dq0
u_d, u_q, u_0 = park(s.x, s.time_track(), fN)
# Plot Fourier transform of the signal
f = rfftfreq(len(s.x), dt)
X = (2/len(s.x))*rfft(s.x)
n = [argmin(abs(f-(fN-fT))), argmin(abs(f-(fN+fT)))]
figure()
vlines(f[n], 0, max(abs(X)), linestyles='dotted')
hlines(abs(X[n]), 0, 50, linestyles='dotted')
plotBode(f, X)
xlim([0, 50])
suptitle('Trimmed signal record')
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

# Find Ud and Uq
Ud = temp_fft(u_d, dt, fT)
Uq = temp_fft(u_q, dt, fT)
print('Calculations:')
print('|Ud| = {0}, ∠Ud = {1}°'.format(abs(Ud), angle(Ud)*180/pi))
print('|Uq| = {0}, ∠Uq = {1}°'.format(abs(Uq), angle(Uq)*180/pi))
