#!/usr/bin/env python

import numpy as np
from math import exp, sqrt, pi
from params import *

Count = 1000

#end of data

t0 = sqrt(6.0)/(pi*freq)
t0 = t0 * 1.3
dt = (3.0*t0)/((float)(Count))

f = open( "SOURCES/impulse.txt", 'wt')

for i in range (Count+1):
	t = i*dt - t0
	func = pi*freq*t
	func = func*func
	func = (1.0 - 2.0*func)*exp(0.0 - func)
	f.write("%s %s \n" %(t + t0, func))

f.close()

