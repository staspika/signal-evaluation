import numpy
from pylab import *

class Signal(object):

    def __init__(self, time=[], data=[]):
        self.time = time
        self.values = data

    def reject_outliers(self, m = 2.):
        d = numpy.abs(self.values - numpy.median(self.values))
        mdev = numpy.median(d)
        s = d/mdev if mdev else 0.
        return Signal(self.time[s<m], self.values[s<m])

    def reject_nan(self):
        mask = [x != nan for x in self.values]
        return Signal(self.time[mask], self.values[mask])

    def cleanup(self):
        return self.reject_nan().reject_outliers(1000000)
