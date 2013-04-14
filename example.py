

import pyga


class F6Individual(Individual):
    __rng = Random()
    _range = (-100.0, 100.0)
    _elements = [Segment(-100.0, 100.0, 22), Segment(-100.0, 100.0, 22)]
    def init(self):
        self.chromo = BinaryChromosome(self._elements)
        self.chromo.encode([self.__rng.randint(self._range), self.__rng.randint(self._range)])
	
    def f(self, x, y):
        return 0.5 - (math.sin(math.sqrt(x*x+y*y))**2.0 - 0.5)/(1.0 + 0.001*(x*x + y*y))**2.0
	
    def evaluate(self):
        u, v = self.chromo.decode()
        self.fitness = self.f(u, v)
        return self.fitness
	


class DummyIndividual(Individual):
    __rng = Random()
    _range = (-10.0, 10.0)
    _elements = [ Segment(-10.0, 10.0, 22), Segment(-10.0, 10.0, 22) ]
    def init(self):
        self.chromo = BinaryChromosome(self._elements)
        self.chromo.encode([self.__rng.randint(self._range), self.__rng.randint(self._range)])
	
    def f(self, x, y):
        return (200 - y*y -x*x)/200.0
	
    def evaluate(self):
        u, v = self.chromo.decode()
        self.fitness = self.f(u, v)
        return self.fitness
	



crossover = pyga.BinaryCrossoverOperator(0.6)
mutation  = pyga.FlipBinaryMutationOperator(0.1)

ga = pyga.ClassicGA(crossover=crossover, mutation=mutation, clazz=pyga.DummyIndividual,
                    rounds=20, popsize=200, generations=500)
ga.run()

from pylab import *

aux = zeros(500, typecode='d')
for i in ga.statistics.keys():
    for j, k in enumerate(ga.statistics[i]['online']):
        aux[j] += k/float(20)
  
figure()
plot(aux)

aux = zeros(500, typecode='d')
for i in ga.statistics.keys():
    for j, k in enumerate(ga.statistics[i]['offline']):
        aux[j] += k/float(20)

figure()
plot(aux)

show()

# vim: expandtab sw=4 ts=4
