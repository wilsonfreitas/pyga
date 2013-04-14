#!/usr/bin/env python

from numpy import *

def schafferF6(x, y):
   t1 = math.sin(math.sqrt(x**2 + y**2));
   t2 = 1.0 + 0.001*(x**2 + y**2);
   score = 0.5 + (t1*t1 - 0.5)/(t2*t2)
   return score

X = linspace(-100, 100, 5000)
Y = linspace(-100, 100, 100)

# for x in X:
# 	for y in Y:
# 		print x, y, schafferF6(x,y);

for x in X:
	print x, schafferF6(x,0);

