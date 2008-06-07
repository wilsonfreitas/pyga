
import pyga

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
