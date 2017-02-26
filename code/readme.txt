To run GA:

GA usage: ./main fileName numberOfIndividuals selectionType crossoverType crossoverProbability mutationProbability numberOfGenerations ga 

To run PBIL:

PBIL usage: ./main fileName numberOfIndividuals positiveLearningRate negativeLearningRate mutationProbability mutationAmount numberOfGenerations p

We follow the command line argument requests as laid out in the handout.

We do not have much error checking on inputs.

Everything works fairly accurately, but our algorithms run a little slow because of a bottleneck in one of our shared functions.
