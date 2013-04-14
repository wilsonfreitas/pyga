#!/usr/bin/env python
# encoding: utf-8
"""
individual.py

Created by Wilson on 2010-11-20.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

from numpy import random as npr
import numpy as np

__all__ = ['Operator', 'Individual', 'Population']

class Operator(object):
	def __init__(self, rate, fitness):
		self.rate = rate
		self.fitness = fitness
	

class Individual(object):
	def __init__(self, lazyeval=False, chromo=None):
		self.chromo = chromo
		self.fitness = None
		if not lazyeval:
			self.initialize()
	
	def initialize(self):
		raise NotImplementedError("This method must be implemented!")
	
	def evaluate(self):
		raise NotImplementedError("This method must be implemented!")
	
	def __cmp__(self, other):
		return cmp( self.evaluate(), other.evaluate() )
	
	def __str__(self):
		return str( self.chromo )
	
	def __repr__(self):
		return str( self.chromo )
	


class Population(object):
	def __init__(self, size, clazz):
		self.size = size
		self.clazz = clazz
		self.individuals = []
		while len(self.individuals) < size:
			indie = clazz()
			if not indie in self.individuals:
				self.individuals.append(indie)
		# self.individuals = [clazz() for i in range(size)]
		self.roulette = np.zeros((size, 2))
		self.sumaptitude = 0.0
	
	def select(self):
		low = 0
		high = len(self.roulette) - 1
		apt = npr.random() * self.sumaptitude
		while low <= high:
			mid = (low + high)/2
			if apt < self.roulette[mid,0]:
				high = mid - 1
			elif apt > self.roulette[mid,1]:
				low = mid + 1
			else:
				return self.individuals[mid]
		# if none individual is selected raise exception
		raise Exception('None individual selected.')
	
	def rank(self): # build ruler
		temp = 0.0
		for i, indie in enumerate(self.individuals):
			self.roulette[i,:] = (temp, temp+indie.fitness)
			temp += indie.fitness
		self.sumaptitude = temp
	
	def performance(self):
		online = 0.0
		offline = self.individuals[0].fitness
		for indie in self.individuals:
			online += indie.fitness/float(self.size)
			if indie.fitness > offline:
				offline = indie.fitness
		return online, offline
	
	def __repr__(self):
		s = ''
		for i in self.individuals:
			s += str(i) + '\n'
		return s
	
	def __len__(self):
		return self.size
	

class ClassicGA(object):
	def __init__(self, *args, **kwargs):
		for k, v in kwargs.iteritems():
			setattr(self, k, v)
		if not hasattr(self, 'popsize'):
			self.popsize = 100
		if not hasattr(self, 'generations'):
			self.generations = 50
		if not hasattr(self, 'rounds'):
			self.rounds = 20
		self.statistics = {}
		self.op_roulette = np.zeros( (len(self.operators), 2) )
		self._op_rank()
	
	def _op_rank(self):
		temp = 0.0
		for i, op in enumerate(self.operators):
			self.op_roulette[i,:] = (temp, temp+op.fitness)
			temp += op.fitness
		self.op_aptitude = temp
	
	def _op_select(self):
		low = 0
		high = len(self.op_roulette) - 1
		apt = npr.random() * self.op_aptitude
		while low <= high:
			mid = (low + high)/2
			if apt < self.op_roulette[mid,0]:
				high = mid - 1
			elif apt > self.op_roulette[mid,1]:
				low = mid + 1
			else:
				return self.operators[mid]
		# if none individual is selected raise exception
		raise Exception('None operator selected.')
	
	def run(self):
		_round = 0
		while _round < self.rounds:
			self.statistics[_round] = {}
			self.statistics[_round]['online'] = []
			self.statistics[_round]['offline'] = []
			# creates population
			pop = Population(self.popsize, self.clazz)
			_gen = 0
			while _gen < self.generations:
				# create roulette
				pop.rank()
				# compute performances
				online, offline = pop.performance()
				# log performance
				self.statistics[_round]['online'].append(online)
				self.statistics[_round]['offline'].append(offline)
				# reproduce (crossover and mutation) generate next population
				individuals = []
				while len(individuals) < len(pop): # select parents for the next generation
					i1 = pop.select()
					i2 = pop.select()
					self._op_select()(i1, i2)
					individuals.append(i1)
					individuals.append(i2)
				# reset pop
				pop.individuals = individuals
				# next generation
				_gen += 1
			# next round
			_round += 1
	

