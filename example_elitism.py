
import pyga

crossover = pyga.BinaryCrossoverOperator(0.7)
mutation  = pyga.FlipBinaryMutationOperator(0.1)

ga = pyga.ElitistGA(crossover=crossover, mutation=mutation, clazz=pyga.F6Individual,
                    rounds=20, popsize=100, generations=200)
ga.run()

from pylab import *

aux = zeros(ga.generations, typecode='d')
for i in ga.statistics.keys():
    for j, k in enumerate(ga.statistics[i]['online']):
        aux[j] += k/float(ga.rounds)
  
figure()
plot(aux[1:])

aux = zeros(ga.generations, typecode='d')
for i in ga.statistics.keys():
    for j, k in enumerate(ga.statistics[i]['offline']):
        aux[j] += k/float(ga.rounds)

figure()
plot(aux[1:])

show()

# vim: expandtab sw=4 ts=4
