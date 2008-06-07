#!/usr/bin/env python
#
# random - random number generator
#
# Copyright 2005 Wilson Freitas
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, version 2.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You have received a copy of the GNU General Public License along
#   with this program, on the COPYING file.
#
# Multiplicative linear congruencial generator:
#  I_t = I_{t+1}a + c (mod M)
#
# M = 2^31 - 1
# a = 16807
# c = 0
# Period = 2^32 - 2
#

from math import log, ceil, sqrt, sin, cos, pi, exp
import time

TWOPI = 2*pi
LOG2  = log(2)

def random_seed():
    return 2 * long(time.time() * 123456789) - 1

class Random(object):
    __MASK = (1<<31) - 1

    def __init__(self, seed=None, a=16807):
        self.seed = (seed or random_seed())
        self.__a = a
        self.__R0 = self.seed  & self.__MASK

    def rand(self):
        self.__R0 = (self.__R0 * self.__a) & self.__MASK
        return self.__R0

    def random(self):
        return self.rand()/float(self.__MASK)

    def randint(self, a, b):
        n = abs(b-a)
        x = int( ceil(log(n)/LOG2) )
        i = self.rand()>>(31-x)
        while i > n:
            i = self.rand()>>(31-x);
        return i + a

    def biased_coin_toss(self, p):
        return (self.random() < p)

    def coin_toss(self):
        return bool(self.rand() >> 30)

Uniform = Random

class Gaussian(object):
    def __init__(self, seed=None):
        self.seed = (seed or random_seed())
        self.__rng1 = Random(seed=self.seed)
        self.__rng2 = Random(seed=self.seed*987654321)
        self.__next_val = None

    def random(self):
        r = self.__next_val
        self.__next_val = None
        if r is None:
            z = self.__rng1.random()
            w = self.__rng2.random()
            slnz = sqrt(-2*log(z))
            r = slnz*sin(TWOPI*w)
            self.__next_val = slnz*cos(TWOPI*w)
        return r

class LogNormal(object):
    def __init__(self, seed=None):
        self.__rng = Gaussian(seed)

    def random(self):
        return exp(self.__rng.random())

class Exponencial(object):
    def __init__(self, seed=None):
        self.__rng = Random(seed=seed)

    def random(self):
        return -log(self.__rng.random())

