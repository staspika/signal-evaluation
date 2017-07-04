from __future__ import division, unicode_literals, print_function
from pylab import *
from scipy.signal import butter

# Functions concerning dq0 transformation

def get_quadrature_signal(a):
    """
    Generate a signal b in quadrature with signal a.
    @TODO: This is a very time-consuming procedure for long records. Can
    it be done easily?
    """
    # Approach 1: Rotating the signal a by pi/2 in the frequency domain
    A = rfft(a)
    B = A*exp(-1j*pi/2)
    n = len(a)
    b = irfft(B, n)
    # Approach 2: Quarter-a-period delay (Danielsen 2010, p. 76, 5.4.3.)
    #b = a[:-N/4]
    #a = a[N/4:] 
    #t = t[N/4:]
    #t = t - t[0]
    return b

def clarke2park(alpha, beta, gamma, t, fN):
    """
    Clarke-to-Park transformation, as described by Concordia (1951).
    """
    alpha = array(alpha)
    beta = array(beta)
    gamma = array(gamma)
    wN = 2*pi*fN
    theta0 = 0
    theta = array(wN*t) + theta0
    d = cos(theta)*alpha + sin(theta)*beta + 0*gamma
    q = -sin(theta)*alpha + cos(theta)*beta + 0*gamma
    o = 0*alpha + 0*beta + 1*gamma
    return d, q, o

def clarke(x):
    alpha = array(x)
    beta = get_quadrature_signal(alpha) # This leads to a discontinuity at the boundaries
    gamma = 0
    return alpha, beta, gamma

def park(x, t, fN):
    alpha, beta, gamma = clarke(x)
    d, q, o = clarke2park(alpha, beta, gamma, t, fN)
    return d, q, o

def temp_fft(x, dt, fT):
    """
    Get the amplitude of signal x at frequency fT.
    """
    f = rfftfreq(len(x), dt)
    X = (2/len(x))*rfft(x)
    n = argmin(abs(f-fT))
    return X[n]
    
def temp_fft1(x, t, fT):
    dt = average(diff(t))
    f_nyquist = 1/dt/2
    Wn = fT/f_nyquist
    b, a = cheby2(8, 40, Wn, 'low', analog=False)
    x_filtered = filtfilt(b, a, x)
    t_new = arange(min(t), max(t), dt)
    x_new = interp(t_new, t, x_filtered)

def t2dq50388(u_t, i_t, t, fN, fT):
    """
    Find D and Q components as specified in EN 50388.
    """
    dt = average(diff(t))
    u_d, u_q, u_0 = park(u_t, t, fN)
    Ud = temp_fft(u_d, dt, fT)
    Uq = temp_fft(u_q, dt, fT)
    i_d, i_q, i_0 = park(i_t, t, fN)
    Id = temp_fft(i_d, dt, fT)
    Iq = temp_fft(i_q, dt, fT)
    return [[Ud, Uq], [Id, Iq]]
