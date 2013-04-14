#!/usr/bin/env python
# encoding: utf-8
"""
chromossome.py

Created by Wilson on 2010-11-20.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

from numpy import random as npr
from numpy import zeros
from pyga import Operator, Individual
import math

def count_bits(n):
    c = 0
    while n != 0:
        n = n & (n-1)
        c += 1
    return c

def binary_mask(l):
    return long(2**l) - 1

def reset_bit(n, i):
    return n&~(1L<<i)

def set_bit(n, i):
    return n|(1L<<i)

def flip_bit(n, i):
    return n^(1L<<i)

def str_bit(v, l=0):
    if not v:
        s = "0"
    else:
        s = ""
    while v:
        s = str(v&1)+s
        v >>= 1
    if l:
        s = "0"*(l-len(s))+s
    else:
        s = s or "0"
    return s



class BinarySegment(object):
	# begin, end reffers to range value
	def __init__(self, length, rng):
		self.begin = rng[0]
		self.end = rng[1]
		self.length = length
		self.mask = binary_mask(length)
	
	def __len__(self):
		return self.length
	


class BinaryChromosome(object):
	def __init__(self, segments):
		self.segments = segments
		self.length = sum(len(n) for n in segments)
		self.mask = binary_mask(self.length)
		self.genes = 0
	
	def encode(self, data):
		s = g = 0L
		seg = self.segments
		for i, value in enumerate(data):
			e = seg[i]
			i = long( math.floor( (value - e.begin)*(2.0**e.length-1.0)/(e.end - e.begin) ) )
			g |= (i & e.mask) << s
			s += e.length
		self.genes = g
	
	def decode(self):
		data = zeros(len(self.segments), int)
		s = 0
		for i, e in enumerate(self.segments):
			n = (self.genes>>s) & e.mask
			data[i] = e.begin + (e.end-e.begin)*n/(2.0**e.length-1.0)
			s += e.length
		return data
	
	def __str__(self):
		return str_bit( self.genes, self.length )
	
	def __len__(self):
		return self.length
	


class BinaryIndividual(Individual):	
	def initialize(self):
		self.chromo = BinaryChromosome(self._segments)
		code = [ npr.randint(e.begin, e.end+1) for e in self._segments ]
		self.chromo.encode( code )
		self.evaluate()
	


class OverwriteBinaryMutationOperator(Operator):
	def __call__(self, indie, indiu):
		if npr.random() < self.rate:
			mutation = [ set_bit, reset_bit ]
			pos = npr.randint(0, indie.chromo.length)
			idx = int(npr.random() < 0.5)
			indie.chromo.genes = mutation[idx](indie.chromo.genes, pos)
		if npr.random() < self.rate:
			mutation = [ set_bit, reset_bit ]
			pos = npr.randint(0, indiu.chromo.length)
			idx = int(npr.random() < 0.5)
			indiu.chromo.genes = mutation[idx](indiu.chromo.genes, pos)
	


class FlipBinaryMutationOperator(Operator):
	def __call__(self, indie, indiu):
		if npr.random() < self.rate:
			pos = npr.randint(indie.chromo.length-1)
			indie.chromo.genes ^= (1L<<pos)
		if npr.random() < self.rate:
			pos = npr.randint(indiu.chromo.length-1)
			indiu.chromo.genes ^= (1L<<pos)
	


class BinaryCrossoverOperator(Operator):
	def __call__(self, indie_1, indie_2):
		if npr.random() < self.rate:
			chromo_1 = indie_1.chromo
			chromo_2 = indie_2.chromo
			pc = npr.randint(chromo_2.length-1)
			mask_1 = long(2**pc)-1
			mask_2 = long(2**(chromo_2.length-pc))-1
			g_1 = chromo_1.genes
			g_2 = chromo_2.genes
			chromo_1.genes = g_1 & mask_1 | ((g_2 >> pc) & mask_2) << pc
			chromo_2.genes = g_2 & mask_1 | ((g_1 >> pc) & mask_2) << pc
			indie_1.chromo = chromo_1
			indie_2.chromo = chromo_2
	

	

