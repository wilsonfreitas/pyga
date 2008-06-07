#!/usr/bin/env python
#
# welga - genetic algorithms
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

import math
from rng import Random # use default random
from copy import copy

__all__ = ["Individual", "Segment", "BinaryChromosome", "ClassicGA",
           "SteadyStateGA", "DummyIndividual"]

################################################################################
# Functions
################################################################################

# XXX to implement Hamming distance
# TODO create a class to handle BitString issues
def mask_factory(l):
    return long(2**l) - 1

def reset_bit(n, i):
    return n&~(1L<<i)

def set_bit(n, i):
    return n|(1L<<i)

def flip_bit(n, i):
    return n^(1L<<i)

################################################################################
# Exceptions
################################################################################

class WrongDataLength(Exception):
    def __init__(self, length, sum):
        self.length = length
        self.sum = sum

    def __str__(self):
        return 'Soma dos tamanhos dos elementos do cromossomo eh diferente do \
                tamanho do cromossomo: not %d == %d' % (self.length, self.sum)

class InvalidValueError(TypeError):
    def __init__(self, element, value):
        self.element = element
        self.value = value

    def __str__(self):
        return 'Valor for do intervalo do elemento do cromossomo: not %g <= %g <= %g' %\
                (self.element.begin, self.value, self.element.end)

class WrongType(Exception):
    pass

################################################################################
# BinaryChromosome
################################################################################

class Segment(object):
    def __init__(self, begin, end, length):
        self.begin = begin
        self.end = end
        self.length = length
        self.mask = mask_factory(self.length)

class BinaryChromosome(object):
    __rng = Random()
    def __init__(self, elements):
        self.elements = elements
        t = 0
        for n in elements:
            t += n.length
        self.length = t
        self.genes = 0L
        self.mask = mask_factory(self.length)

    def encode(self, data):
        if type(data) in ( list, tuple ):
            self.encode_sequence(data)
        elif type(data) in ( int, float ):
            self.encode_number(data)
        else:
            raise WrongType()

    def encode_sequence(self, data):
        if not len(data) == len(self.elements):
            raise WrongDataLength(len(self.elements), len(data))
        s = g = 0L
        for i, item in enumerate(data):
            e = self.elements[i]
            if not e.begin <= item <= e.end:
                raise InvalidValueError(e, item)
            i = long( math.floor( (item - e.begin)*(2.0**e.length-1.0)/(e.end - e.begin) ) )
            g |= (i & e.mask) << s
            s += e.length
        self.genes = g

    def encode_number(self, data):
        if not len(self.elements) == 1:
            raise WrongDataLength(len(self.elements), 1)
        e = self.elements[0]
        if not e.begin <= data <= e.end:
            raise InvalidValueError(e, data)
        i = long( math.floor( (data - e.begin)*(2.0**e.length-1.0)/(e.end - e.begin) ) )
        self.genes = i & e.mask

    def decode(self):
        if len(self.elements) == 1:
            e = self.elements[0]
            n = self.genes & e.mask
            return e.begin + (e.end-e.begin)*n/(2.0**e.length-1.0)
        else:
            data = []
            s = 0
            for e in self.elements:
                n = (self.genes>>s) & e.mask
                v = e.begin + (e.end-e.begin)*n/(2.0**e.length-1.0)
                data.append(v)
                s += e.length
            return data

    def crossover(self, other, rate):
        if self.__rng.biased_coin_toss(rate):
            pc = self.__rng.randint(self.length)
            mask_1 = mask_factory(pc)
            mask_2 = mask_factory(self.length - pc)
            c_1 = BinaryChromosome(self.elements)
            c_1.genes = self.genes & mask_1 | ((other.genes >> pc) & mask_2) << pc
            c_2 = BinaryChromosome(self.elements)
            c_2.genes = other.genes & mask_1 | ((self.genes >> pc) & mask_2) << pc
            return c_1, c_2
        else:
            return self, other

    def mutate_overwrite_one_point(self, rate):
        if self.__rng.coin_toss():
            mutation = set_bit
        else:
            mutation = reset_bit
        p = self.__rng.randint(self.length)
        if self.__rng.biased_coin_toss(rate):
            self.genes = mutation(self.genes, p)

    def mutate_flip_one_point(self, rate):
        p = self.__rng.randint(self.length)
        if self.__rng.biased_coin_toss(rate):
            self.genes = flip_bit(self.genes, p)

    def mutate_flip(self, rate):
        for i in xrange(self.length):
            if self.__rng.biased_coin_toss(rate):
                self.genes = flip_bit(self.genes, i)

    def mutate_overwrite(self, rate):
        if self.__rng.coin_toss():
            mutation = set_bit
        else:
            mutation = reset_bit
        for i in xrange(self.length):
            if self.__rng.biased_coin_toss(rate):
                self.genes = mutation(self.genes, i)
    mutate = mutate_flip_one_point

    def __str__(self):
        s = ''
        for i in xrange(self.length):
            s += str((self.genes>>(self.length-i-1)) & 1)
        return s

################################################################################
# Individual
################################################################################

class Individual(object):
    def __init__(self):
        self.chromo = None
        self.fitness = 0.0
        self.aptitude = 0.0

    def mutate(self, rate):
        self.chromo.mutate(rate)

    def crossover(self, other, rate):
        cc_1, cc_2 = self.chromo.crossover(other.chromo, rate)
        child_1 = copy(self)
        child_1.chromo = cc_1
        child_1.fitness = 0.0
        child_1.aptitude = 0.0
        child_2 = copy(self)
        child_2.chromo = cc_2
        child_2.fitness = 0.0
        child_2.aptitude = 0.0
        return child_1, child_2

    def reproduce(self, other, crossover_rate, mutation_rate):
        child_1, child_2 = self.crossover(other, crossover_rate)
        child_1.mutate(mutation_rate)
        child_2.mutate(mutation_rate)
        return child_1, child_2

    def load_at_random(self):
        raise NotImplementedError("This method must be implemented!")

    def evaluate(self):
        raise NotImplementedError("This method must be implemented!")

    def __cmp__(self, other):
        if self.evaluate() > other.evaluate():
            return 1
        else:
            return -1
        #return int(self.fitness > other.fitness) or -1

    def __str__(self):
        return '%s --> %20.5f' % (str(self.chromo), self.evaluate())

class DummyIndividual(Individual):
    __rng = Random()
    _range = (-100.0, 100.0)
    _elements = [Segment(-100.0, 100.0, 22), Segment(-100.0, 100.0, 22)]
    def load_at_random(self):
        self.chromo = BinaryChromosome(self._elements)
        self.chromo.encode([self.__rng.irand(self._range), self.__rng.irand(self._range)])

    def f(self, x, y):
        return 0.5 - (math.sin(math.sqrt(x*x+y*y))**2.0 - 0.5)/(1.0 + 0.001*(x*x + y*y))**2.0

    def evaluate(self):
        u, v = self.chromo.decode()
        fitness = self.fitness = self.f(u, v)
        return fitness


################################################################################
# ClassicGA
################################################################################

class ClassicGA(object):
    # ordem: crossover_rate, mutation_rate, IndividualClazz, popsize, generations, rounds
    # (self, crossover_rate, mutation_rate, IndividualClazz, popsize=100, generations=50, rounds=50):
    __slots__ = ['crossover_rate', 'mutation_rate', 'clazz', 'popsize',\
            'generations', 'rounds', 'statistics', 'elitist', 'linfit', 'best']
    def __init__(self, *args, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
        if not hasattr(self, 'popsize'):
            self.popsize = 100
        if not hasattr(self, 'generations'):
            self.generations = 50
        if not hasattr(self, 'rounds'):
            self.rounds = 20
        if not hasattr(self, 'elitist'):
            self.elitist = False
        if not hasattr(self, 'linfit'):
            self.linfit = None
        self.statistics = None

    def run(self):
        _round = 0
        statistics = {}
        best = None
        while _round < self.rounds:
            statistics[_round] = {}
            statistics[_round]['online'] = []
            statistics[_round]['offline'] = []
            population = ClassicPopulation(self.popsize, self.clazz,\
                    self.crossover_rate, self.mutation_rate,\
                    elitist=self.elitist, linfit=self.linfit)
            # inicializa populacao
            # avalia populacao
            population.initialize()
            population.rank()
            # LOOP
            generation = 0
            # optimization tricks
            select = population.select
            rank = population.rank
            online_append = statistics[_round]['online'].append
            offline_append = statistics[_round]['offline'].append
            online = population.get_online_performance
            offline = population.get_offline_performance
            genbest = population.individuals[0]
            # optimization tricks
            while generation < self.generations:
                # seleciona populacao
                # reproduz (crossover e mutacao)
                select()
                # avalia populacao
                rank()
                online_append(online())
                offline_append(offline())
                generation += 1
                if population.best.fitness > genbest.fitness:
                    genbest = population.best
            _round += 1
            if best is None or genbest.fitness > best.fitness:
                best = genbest
            del population
        self.statistics = statistics
        self.best = best
        return best

class ClassicPopulation(object):
    __rng = Random()
    def __init__(self, size, clazz, crossover_rate, mutation_rate, elitist=False, linfit=None):
        self.size = size
        self.clazz = clazz
        self.individuals = []
        self.sum_fitness = 0.0
        self.sum_aptitude = 0.0
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.roulette = []
        self.__online_performance = []
        self.__offline_performance = []
        self.best = None
        self.elitist = elitist
        self.linfit = linfit

    def initialize(self):
        sum = 0.0
        clazz = self.clazz
        append = self.individuals.append
        # criar um individuo fora do loop para poder encontrar o melhor da populacao inicial
        indie = clazz()
        indie.load_at_random()
        fitness = indie.evaluate()
        bestfit = fitness
        best = indie
        sum += fitness
        append(indie)
        for i in xrange(self.size-1):
            indie = clazz()
            indie.load_at_random()
            fitness = indie.evaluate()
            if fitness > bestfit:
                bestfit = fitness
                best = indie
            sum += fitness
            append(indie)
        self.sum_fitness = sum
        self.best = best

    def rank(self):
        if self.linfit is not None:
            self.individuals.sort()
            compute = lambda: linfit[0] + (linfit[1]-linfit[0])*i/(size-1.0)
        else:
            compute = lambda: indie.fitness
        sum = 0.0
        del self.roulette[:]
        append = self.roulette.append
        linfit = self.linfit
        size = self.size
        # optimization tricks
        for i, indie in enumerate(self.individuals):
            aptitude = compute()
            append((sum, sum+aptitude))
            sum += aptitude
            indie.aptitude = aptitude
        self.sum_aptitude = sum

    def set_online_performance(self, value):
        self.__online_performance.append(value)
    def get_online_performance(self):
        l = len(self.__online_performance)
        return self.__online_performance[l-1]
    online_performance = property(get_online_performance, set_online_performance)

    def set_offline_performance(self, value):
        self.__offline_performance.append(value)
    def get_offline_performance(self):
        l = len(self.__offline_performance)
        return self.__offline_performance[l-1]
    offline_performance = property(get_offline_performance, set_offline_performance)

    def select(self):
        i = 0
        sum = 0.0
        newpop = []
        append = newpop.append
        selection = self.selection
        crossover = self.crossover_rate
        mutation = self.mutation_rate
        worstfit = bestfit = None
        worst = best = None
        # optimization tricks
        while i < self.size:
            p_1, p_2 = selection(), selection()
            c_1, c_2 = p_1.reproduce(p_2, crossover, mutation)
            sum += c_1.evaluate()
            sum += c_2.evaluate()
            append(c_1)
            append(c_2)
            if best is None:
                if c_1.fitness > c_2.fitness:
                    best, bestfit = c_1, c_1.fitness
                    worst, worstfit = c_2, c_2.fitness
                else:
                    best, bestfit = c_2, c_2.fitness
                    worst, worstfit = c_1, c_1.fitness
            else:
                if bestfit < c_1.fitness:
                    best, bestfit = c_1, c_1.fitness
                elif worstfit > c_1.fitness:
                    worst, worstfit = c_1, c_1.fitness
                if bestfit < c_2.fitness:
                    best, bestfit = c_2, c_2.fitness
                elif worstfit > c_2.fitness:
                    worst, worstfit = c_2, c_2.fitness
            i += 2
        if self.elitist:
            if bestfit > self.best.fitness:
                self.best = best
            else:
                sum += self.best.fitness - worst.fitness
                n = newpop.index(worst)
                newpop[n] = self.best
        else:
            self.best = best
        self.sum_fitness = sum
        self.online_performance = sum/self.size
        self.offline_performance = self.best.fitness
        del self.individuals[:]
        self.individuals.extend(newpop)

    def selection(self):
        low = 0
        high = len(self.roulette) - 1
        apt = self.__rng.rand() * self.sum_aptitude
        while low <= high:
            mid = (low + high)/2
            if apt < self.roulette[mid][0]:
                high = mid - 1
            elif apt > self.roulette[mid][1]:
                low = mid + 1
            else:
                return self.individuals[mid]
        raise Exception()

################################################################################
# SteadyStateGA
################################################################################

class SteadyStateGA(object):
    # ordem: crossover_rate, mutation_rate, IndividualClazz, popsize, generations, rounds
    # (self, crossover_rate, mutation_rate, IndividualClazz, popsize=100, generations=50, rounds=50):
    __slots__ = ['crossover_rate', 'mutation_rate', 'clazz', 'popsize',\
            'generations', 'rounds', 'statistics', 'gap', 'linfit', 'best']
    def __init__(self, *args, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
        if not hasattr(self, 'popsize'):
            self.popsize = 100
        if not hasattr(self, 'generations'):
            self.generations = 50
        if not hasattr(self, 'rounds'):
            self.rounds = 20
        if not hasattr(self, 'gap'):
            self.gap = 2
        if not hasattr(self, 'linfit'):
            self.linfit = None
        self.statistics = None

    def run(self):
        _round = 0
        statistics = {}
        best = None
        while _round < self.rounds:
            statistics[_round] = {}
            statistics[_round]['online'] = []
            statistics[_round]['offline'] = []
            population = SteadyStatePopulation(self.popsize, self.clazz,\
                    self.crossover_rate, self.mutation_rate,\
                    self.gap, linfit=self.linfit)
            # inicializa e avalia populacao
            population.initialize()
            # ordena e normaliza (se quiser) populacao
            population.rank()
            # LOOP
            generation = 0
            gengap = 0
            # optimization tricks
            popsize = self.popsize
            select = population.select
            merge = population.merge
            rank = population.rank
            append_online = statistics[_round]['online'].append
            append_offline = statistics[_round]['offline'].append
            online = population.get_online_performance
            offline = population.get_offline_performance
            # optimization tricks
            while generation < self.generations:
                # seleciona populacao
                # reproduz (crossover e mutacao)
                gengap = gengap + select()
                # avalia populacao
                merge()
                rank()
                generation += gengap%popsize == 0
                if gengap % popsize == 0:
                    append_online(online())
                    append_offline(offline())
            _round += 1
            if best is None or population.best.fitness > best.fitness:
                best = population.best
            del population
        self.statistics = statistics
        self.best = best
        return best

class SteadyStatePopulation(object):
    __rng = Random()
    def __init__(self, size, clazz, crossover_rate, mutation_rate, gap, linfit=None):
        self.size = size
        self.IndividualClazz = clazz
        self.individuals = []
        self.kids = []
        self.sum_fitness = 0.0
        self.sum_aptitude = 0.0
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.roulette = []
        self.__online_performance = []
        self.__offline_performance = []
        self.best = None
        self.gap = gap
        self.linfit = linfit

    def initialize(self):
        sum = 0.0
        clazz = self.IndividualClazz
        append = self.individuals.append
        size = self.size
        # criar um individuo fora do loop para poder encontrar o melhor da populacao inicial
        indie = clazz()
        indie.load_at_random()
        fitness = indie.evaluate()
        bestfit = fitness
        best = indie
        sum += fitness
        append(indie)
        i = 1
        while i < size:
            indie = clazz()
            indie.load_at_random()
            fitness = indie.evaluate()
            if fitness > bestfit:
                bestfit = fitness
                best = indie
            sum += fitness
            append(indie)
            i += 1
        self.best = best
        self.sum_fitness = sum

    def rank(self):
        self.individuals.sort()
        self.best = self.individuals[self.size-1]
        self.offline_performance = self.best.fitness
        if self.linfit is not None:
            compute = lambda: linfit[0] + (linfit[1]-linfit[0])*i/(size-1.0)
        else:
            compute = lambda: indie.fitness
        del self.roulette[:]
        sum = 0.0
        append = self.roulette.append
        linfit = self.linfit
        size = self.size
        for i, indie in enumerate(self.individuals):
            aptitude = compute()
            append((sum, sum+aptitude))
            sum += aptitude
            indie.aptitude = aptitude
        self.sum_aptitude = sum

    def set_online_performance(self, value):
        self.__online_performance.append(value)
    def get_online_performance(self):
        l = len(self.__online_performance)
        return self.__online_performance[l-1]
    online_performance = property(get_online_performance, set_online_performance)

    def set_offline_performance(self, value):
        self.__offline_performance.append(value)
    def get_offline_performance(self):
        l = len(self.__offline_performance)
        return self.__offline_performance[l-1]
    offline_performance = property(get_offline_performance, set_offline_performance)

    def select(self):
        newpop = []
        append = newpop.append
        selection = self.selection
        gap = self.gap
        crossover = self.crossover_rate
        mutation = self.mutation_rate
        # optimization tricks
        i = 0
        while i < gap:
            p_1, p_2 = selection(), selection()
            c_1, c_2 = p_1.reproduce(p_2, crossover, mutation)
            c_1.evaluate()
            c_2.evaluate()
            append(c_1)
            append(c_2)
            i += 2
        del self.kids[:]
        self.kids.extend(newpop)
        return i

    def evaluate(self):
        sum = 0.0
        for indie in self.individuals:
            sum += indie.evaluate()
        self.sum_fitness = sum
        return sum

    def selection(self):
        low = 0
        high = len(self.roulette) - 1
        apt = self.__rng.rand() * self.sum_aptitude
        while low <= high:
            mid = (low + high)/2
            if apt < self.roulette[mid][0]:
                high = mid - 1
            elif apt > self.roulette[mid][1]:
                low = mid + 1
            else:
                return self.individuals[mid]
        raise Exception()

    def merge(self):
        individuals = self.individuals
        kids = self.kids
        delta = 0.0
        # optimization tricks
        for i in xrange(self.gap):
            delta += kids[i].fitness - individuals[i].fitness
            individuals[i] = kids[i]
        self.sum_fitness += delta
        self.online_performance = self.sum_fitness/self.size

    def __len__(self):
        return self.size

################################################################################
# Tests
################################################################################

import struct

def count_bits(n):
    c = 0
    while n != 0:
        n = n & (n-1)
        c += 1
    return c

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

def str_bit1(v):
		s = ""
    while v:
        s = str(v&1)+s
        v >>= 1
    return s or "0"

def test(obj, test, msg=''):
    s = '%3s : %-30s : %s'
    if test:
        print s % ('OK', obj.__name__, msg)
    else:
        print s % ('BAD', obj.__name__, msg)

if __name__ == '__main__':
    mask = mask_factory(10);
    test(mask_factory, len(str_bit(mask)) == 10, 'comprimento da mascara em 32 bits')
    mask = mask_factory(100);
    test(mask_factory, len(str_bit(mask)) == 100, 'comprimento da mascara com + de 32 bits')

    e = Segment(-10, 10, 32)
    test(Segment, len(str_bit(e.mask)) == 32, 'construtor')

    c = BinaryChromosome([e])
    test(BinaryChromosome, c.length == 32, 'comprimento do cromossomo no construtor')

    c.encode(10)
    test(BinaryChromosome, len(str_bit(c.genes)) == 32, 'encode da maior solucao' )

    c.encode(-10)
    test(BinaryChromosome, c.genes == 0, 'encode da menor solucao' )

    c.encode(0)
    test(BinaryChromosome, len(str_bit(c.genes)) == 32-1, 'encode da solucao do meio' )

    try:
        c.encode(20)
        test(BinaryChromosome, False, 'encode de numero fora do range' )
    except InvalidValueError:
        test(BinaryChromosome, True, 'decode de numero fora do range' )

    # TODO: implementar cromossomo por presisao para melhorar isso
    c.encode(5.5)
    d = c.decode()
    test(BinaryChromosome, d-5.5 < 0.1, 'decode de numero de ponto flutuante' )

    c_1 = BinaryChromosome([e])
    c_1.encode(5.5)
    c_2 = BinaryChromosome([e])
    c_2.encode(5.5)
    test(BinaryChromosome, c_1.genes == c_2.genes, 'cromossomos identicos: %d == %d' % (c_1.genes, c_2.genes))
    c_11, c_22 = c_1.crossover(c_2, 1)
    test(BinaryChromosome, c_11.length == 32, 'comprimento do cromossomo apos crossover')
    test(BinaryChromosome, c_22.length == 32, 'comprimento do cromossomo apos crossover')
    test(BinaryChromosome, c_11.genes == c_1.genes, 'crossover de cromossomos identicos: %d == %d' % (c_11.genes, c_1.genes))
    test(BinaryChromosome, c_22.genes == c_2.genes, 'crossover de cromossomos identicos: %d == %d' % (c_22.genes, c_2.genes))
    c_1.encode(-3)
    print c_1, c_1.genes
    c_2.encode(5.5)
    print c_2, c_2.genes
    c_11, c_22 = c_1.crossover(c_2, 1)
    print c_11, c_11.genes
    print c_22, c_22.genes
    del c_11, c_22, c_1, c_2
    c.mutate(1)
    test(BinaryChromosome, c.genes == 0 or len(str_bit(c.genes)) == 32, 'mutacao com  taxa 100%')

    class TestIndividual(Individual):
        rng = Random()
        _range = (-10.0, 10.0)
        _elements = [Segment(-10.0, 10.0, 30)]
        def load_at_random(self):
            self.chromo = BinaryChromosome(self._elements)
            self.chromo.encode(self.rng.irand(self._range))

        def f(self, x):
            return (10.0-x)*(10.0+x)

        def evaluate(self):
            v = self.chromo.decode()
            fitness = self.fitness = self.f(v)
            return fitness

    i_1 = TestIndividual()
    i_1.load_at_random()
    i_1.evaluate()
    print i_1
    i_2 = TestIndividual()
    i_2.load_at_random()
    i_2.evaluate()
    print i_2
    i_11, i_22 = i_1.crossover(i_2, 1)
    print i_11
    print i_22

    del i_1, i_2, i_11, i_22

