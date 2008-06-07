
import pyga

pop = pyga.ClassicPopulation(20, pyga.DummyIndividual)
pop.initialize()
pop.rank()
for i in range(10):
	print pop.selection()

