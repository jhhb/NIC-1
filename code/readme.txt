GA v PBIL on MAXSAT 
James Boyle, Justin Wallace, Martin Beranrd

To run GA:

GA usage: ./main fileName numberOfIndividuals selectionType crossoverType crossoverProbability mutationProbability numberOfGenerations ga 

To run PBIL:

PBIL usage: ./main fileName numberOfIndividuals positiveLearningRate negativeLearningRate mutationProbability mutationAmount numberOfGenerations p

Settings that we used to test (and our result):

PBIL: NumInd=100 plr=0.05 nlr=0.075 mProb=0.02 mAmount=0.1 NumIt=500 average fitness=585.9
GA: NumInd=100 selection=bs crossover=1c cProb=0.7 mProb=0.01 NumGen=500 average fitness=597.2

We follow the command line argument requests as laid out in the handout.

We do not have much error checking on inputs.

Everything works fairly accurately, but our algorithms run a little slow because of a bottleneck in one of our shared functions.
