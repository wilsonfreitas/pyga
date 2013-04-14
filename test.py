#!/usr/local/bin/python
# encoding: utf-8
"""
test.py

Created by Wilson on 2010-11-20.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import unittest
import numpy.random as npr
import numpy as np
from scipy.stats import chisquare
from pyga.binary import BinarySegment, BinaryChromosome, BinaryIndividual, str_bit
from pyga.binary import FlipBinaryMutationOperator, BinaryCrossoverOperator, OverwriteBinaryMutationOperator
from pyga import Individual, Population, ClassicGA
from math import sqrt


class TestBinaryChromosome(unittest.TestCase):
	def setUp(self):
		pass
	
	def test_BinaryChromosomeInitialization(self):
		"""Testing BinaryChromosome initialization"""
		chrome = BinaryChromosome( (BinarySegment(10, (-10,10)),) )
		self.assertEqual(chrome.length, 10)
		self.assertEqual(str_bit(chrome.mask), '1111111111')
		chrome = BinaryChromosome( (BinarySegment(4, (1,4)), BinarySegment(4, (1,4))) )
		self.assertEqual(chrome.length, 8)
		self.assertEqual(str_bit(chrome.mask), '11111111')
	
	def test_BinaryChromosomeEncode(self):
		"""Testing BinaryChromosome encode method"""
		chrome = BinaryChromosome( (BinarySegment(4, (0,15)),) )
		chrome.encode((0,))
		self.assertEqual(chrome.genes, 0)
		chrome.encode((1,))
		self.assertEqual(str_bit(chrome.genes), '1')
		chrome.encode((2,))
		self.assertEqual(str_bit(chrome.genes), '10')
		chrome.encode((4,))
		self.assertEqual(str_bit(chrome.genes), '100')
		chrome.encode((15,))
		self.assertEqual(str_bit(chrome.genes), '1111')
	
	def test_BinaryChromosomeDecode(self):
		"""Testing BinaryChromosome decode method"""
		chrome = BinaryChromosome( (BinarySegment(4, (0,15)),) )
		chrome.encode((0,))
		self.assertEqual(0, chrome.decode()[0])
		chrome.encode((1,))
		self.assertEqual(1, chrome.decode()[0])
		chrome.encode((2,))
		self.assertEqual(2, chrome.decode()[0])
		chrome.encode((4,))
		self.assertEqual(4, chrome.decode()[0])
		chrome.encode((15,))
		self.assertEqual(15, chrome.decode()[0])
	


class DummyIndividual(BinaryIndividual):
    _segments = [ BinarySegment(4, (0,16)), BinarySegment(4, (0,16)) ]
    def evaluate(self):
        u, v = self.chromo.decode()
        self.fitness = u + v
        return self.fitness
	


class TestBinaryIndividual(unittest.TestCase):
	def setUp(self):
		pass
	
	def test_IndividualInitialization(self):
		"""Testing Individual initialization"""
		indi = DummyIndividual()
		self.assertTrue( 0 <= indi.evaluate() <= 30 )
	


class TestPopulation(unittest.TestCase):
	def setUp(self):
		pass
	
	def test_PopulationInitialization(self):
		"""Testing Population initialization"""
		pop = Population(10, DummyIndividual)
		self.assertTrue( len(pop.individuals) == 10 )
	
	def test_PopulationSelection(self):
		"""Testing Population selection"""
		pop = Population(32, DummyIndividual)
		# print sorted(pop.individuals)
		pop.rank()
		hist = {}
		N = len(pop)
		for i in range(N):
			indie = pop.select()
			try:
				hist[indie.fitness] += 1
			except:
				hist[indie.fitness] = 1
		p_obs = []
		p_exp = []
		for k,v in hist.iteritems():
			p_exp.append(k/pop.sumaptitude)
			p_obs.append(v/float(N))
		chisq, p_value = chisquare(np.array(p_obs), np.array(p_exp))
		self.assertTrue( p_value > 0.99 )
	

class TestClassicGA(unittest.TestCase):
	def setUp(self):
		"""docstring for setUp"""
		pass
	
	def test_ClassicGA_OpSelect(self):
		"""Testing operator selection"""
		crossover = BinaryCrossoverOperator(0.6, 100)
		mutation = FlipBinaryMutationOperator(0.08, 10)
		mutation_2 = OverwriteBinaryMutationOperator(0.05, 5)
		ga = ClassicGA(clazz=DummyIndividual, popsize=20, generations=10, rounds=1,
			operators=(crossover, mutation, mutation_2) )
		print ga.op_roulette
	

if __name__ == '__main__':
	unittest.main()

