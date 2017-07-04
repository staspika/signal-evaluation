from pylab import *
import scipy.stats
from scipy.signal import *

def time_track(n, offset=0, dt=1):
    """
    Computes the time track of a signal, provided the time is 
    equally-discrete.
    
    Parameters
    ----------
    n : integer
        Length of time vector.
    offset : float
        Start time (useful when the start time is not zero).
    step : float
        Time quantum (step of discretization).
    
    Returns
    -------
    t : array_like
        Time vector.
    
    Notes
    -----
    """
    return linspace(offset, offset + dt*n, n)

def guess_res(a):
    """
    Guess the resolution of signal a, provided the signal is 
    evenly-quantized and the deviance in quanta size is 
    normally-distributed.
    
    Parameters
    ----------
    a : array_like
        Signal, whose resolution is to be guessed.
    
    Returns
    -------
    (m, v) : tuple, where
        m : float
            Arithmetic mean of quanta size.
        v : float
            Variance in quanta size.
        
    Notes
    -----
    `guess_res` tends to guess the resolution of the signal quite good 
    when signal vector is large, and there are few cases when the 
    signal's value jumps over several quanta at once. The author knows 
    there are several things that could be done better.
    """
    a = list(set(a)) # Find unique values in signal a
    a.sort() # Sort the array of unique values
    x = diff(a) # Find signal's quanta
    x = list(set(x))
    x.sort()
    x = array(x)

    # Discard "multiple quanta"
    m = scipy.stats.tmean(x) # Find a mean of quanta size
    #v = scipy.stats.tstd(x) # Find a standard deviation in quanta size
    #mask = numpy.array([(n>m-2*v and n<m+2*v) for n in x])
    mask = array([(n>0.9*m and n<1.1*m) for n in x]) # Allow for 10% variation in quanta size.
    outliers = x[invert(mask)]
    x = x[mask]
    #print (outliers/m)
    m = median(x) # Find median of the quanta population
    return m

def guess_period(x):
    """
    Determine the fundamental frequency of a semiperiodic signal. If 
    more than one period of a signal is provided, returns the average 
    frequency.
    
    Works best with long records of single-frequency signals.
    """
    y = correlate(x, x, mode='full')
    y = y[int(len(y)/2):]
    #~ y = cumsum(y) # Smooth the correlation results
    a = argrelmax(y)
    b = diff(a[0])
    return average(b)

def rms(x, t, f, n=1):
    dt = average(diff(t))
    T = 1./f
    N = int(round(T/dt))
    y = []
    for i in range(N, len(t), n):
        y.append(sqrt(sum(x[i-N:i]**2)/N))
    y = array(y)
    return t[N:], y

class Signal:
    def __init__(self, x, dt):
        self.x = array(x) # TODO: Find a better variable name. "Values"?
        self.dt = dt
    def intmean(self):
        y = sqrt(cumsum(x**2)/arange(1, len(x)+1))
        return y
    def time_track(self):
        return time_track(len(self.x), offset=0, dt=self.dt)
    def guess_res(self):
        return guess_res(self.x)
    def guess_period(self):
        return guess_period(self.x)*dt
    def rms(self):
        return rms(self.x, self.time_track(), 1/self.guess_period())

def getfreqsaw(x, t, f0):
    """
    Determine the frequency of a sawtooth signal, given that it 
    contains relatevely little noise and assuming the frequency in the 
    neighbourhood of f0.
    """
    #~ x = x - mean(x)
    y = cumsum(x)
    T0 = 1./f0
    dt = average(diff(t))
    N0 = T0/dt
    n = int(floor(N0/2))
    # Search local maxima in the neighbourhood of half of a period at
    # the assumed frequency f0
    z = argrelmax(y, order = n)
    T = average(diff(t[z]))
    f = 1./T
    return f 
