
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "utility.h"
#include <algorithm>
#include "ga.h"
#include <cmath>


using namespace std;

const double EulerConstant = std::exp(1.0);

vector<vector<int> > rankSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses);
vector<vector<int> > tournamentSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses);
vector<vector<int> > boltzmannSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses, int numberOfClauses);
vector<int> gaMutate(vector<int> child, double probabilityOfMutation);
void runGA(string fileName, int numberOfIndividuals, string selection, string crossover, double crossoverProbability, double mutationProbability, int numberOfGenerations);

int getFitness(vector< vector<int> > cnf, vector<bool> cs);
void printVofV(vector< vector<int> > cnf);
void pbil(string fileName, int numIndividuals, double plr, double nlr, double mutProb, double mutAmount, int numIterations);
vector< vector<bool> > genCS(int numIndividuals, vector<double> probVector);
vector<double> mutate(vector<double> pv, double mutProb, double mutAmount);

bool descending (candidateFitnessAndPosition i, candidateFitnessAndPosition j) { return (i.fitnessScore < j.fitnessScore); }

int main( int argc, const char* argv[] )
{
    /*check these*/
    if(argc != 9){
        cout<<"GA usage: fileName numberOfIndividuals selectionType crossoverType crossoverProbability mutationProbability numberOfGenerations ga"<<endl;
        cout<<"PBIL usage: fileName numberOfIndividuals positiveLearningRate negativeLearningRate mutationProbability mutationAmount numberOfGenerations p"<<endl;

        cout<<"Program exiting"<<endl;
        exit(-1);
    }
    if(strncmp(argv[8],"ga",1) == 0){

        string fileName = argv[1];
        int numberOfIndividuals = atoi(argv[2]);
        string selectionType = argv[3];
        string crossoverType = argv[4];
        double crossoverProbability = stod(argv[5]);
        double mutationProbability = stod(argv[6]);
        int numberOfGenerations = atoi(argv[7]);
        //argv 8 == g
        runGA(fileName, numberOfIndividuals, selectionType, crossoverType, crossoverProbability, mutationProbability, numberOfGenerations);
    }
    else if(strncmp(argv[8],"p",1) == 0) {

        string fileName = argv[1];
        int numberOfIndividuals = atoi(argv[2]);
        double plr = stod(argv[3]);
        double nlr = stod(argv[4]);
        double mutationProbability = stod(argv[5]);
        double mutationAmount = stod(argv[6]);
        int numberOfGenerations = atoi(argv[7]);
        //argv is numberOfGenerations

        pbil(fileName, numberOfIndividuals, plr, nlr, mutationProbability, mutationAmount, numberOfGenerations);
    }
    else{
        exit(-1);
    }
    /*Check these*/

}



/*While we dont have all the individuals we want, we select two random parents from the breeding pool, cross them using a random position in the array, and add to
new array.

We then return the newCandidates
*/
vector<vector<int> > crossoverWrapper(string crossover, int numberOfIndividuals, vector<vector<int> > breedingPool, double crossoverProbability, double mutationProbability){
	
	vector<vector<int> > newCandidates;
	vector<int> child1;
	vector<int> child2;

	//1point crossover
	if(crossover == "1c"){

			int i = 0; //offspring counter and breeding pool iterator

			//generate numberOfIndividuals offspring from the breeding pool
			while(i < numberOfIndividuals - 1){
	
				double randomProbability = (double)rand()/(double)RAND_MAX;
    			//used to be randomProb <= prob
				if(randomProbability <= crossoverProbability){

					int crossPoint = rand() % breedingPool[i].size();
					child1 = onePointCrossover(breedingPool[i], breedingPool[i+1], crossPoint);
					child2 = onePointCrossover(breedingPool[i+1], breedingPool[i], crossPoint);
				}
				else {
					child1 = breedingPool[i];
					child2 = breedingPool[i+1];
				}

				child1 = gaMutate(child1, mutationProbability);
				child2 = gaMutate(child2, mutationProbability);

				newCandidates.push_back(child1);
				newCandidates.push_back(child2);

				i += 2;

				if (newCandidates.size() >= breedingPool.size()) {
					if (breedingPool.size() % 2 == 0) {
						break;
					} else {
						newCandidates.push_back(breedingPool[i]);
						i += 1;
					}
				}
			}
		}

		//Uniform Crossover
		else {

			int i = 0;

			//generate numberOfIndividuals offspring from the breeding pool
			while(i < numberOfIndividuals - 1){
				
				double randomProbability = (double)rand()/(double)RAND_MAX;
    			
				if(randomProbability <= crossoverProbability){

					child1 = uniformCrossover(breedingPool[i], breedingPool[i+1]);
					child2 = uniformCrossover(breedingPool[i], breedingPool[i+1]);
				}
				else {

					child1 = breedingPool[i];
					child2 = breedingPool[i+1];
				}

				child1 = gaMutate(child1, mutationProbability);
				child2 = gaMutate(child2, mutationProbability);

				newCandidates.push_back(child1);
				newCandidates.push_back(child2);

				i += 2;

				if (newCandidates.size() >= breedingPool.size()) {
					if (breedingPool.size() % 2 == 0) {
						break;
					} else {
						newCandidates.push_back(breedingPool[i]);
						i += 1;
					}
				}
			}
		}
		return newCandidates;
}

void runGA(string fileName, int numberOfIndividuals, string selection, string crossover, double crossoverProbability, double mutationProbability, int numberOfGenerations){

	int numberOfVariables = 0;	//from readFile
	int numberOfClauses = 0;	//comes from readFile
	int generationCounter = 0;	//counts our generations
	vector<vector<int> > vectorOfClauses = readFile(fileName, &numberOfVariables, &numberOfClauses);	//clauses in vector<vector<int>>
	vector<vector<int> > breedingPool;	//breeding pool after selection occurs
	vector<vector<int> > candidates = genSolutions(numberOfVariables, numberOfIndividuals); //starting candidate solutions
	int maxFitness = 0;
	int indexOfBestThisGenerationFitness;
	double bestFitnessSoFar = 0;
	int indexOfCandidateWithBestFitness;
    int iterationWhereBestIsFound;

	while(generationCounter < numberOfGenerations){
		int thisGenerationFitness = 0;

		for(int i = 0; i < candidates.size();i++){
			int fitness = getFitnessIntegers(vectorOfClauses, candidates[i]);

			if(fitness > thisGenerationFitness){
				thisGenerationFitness = fitness;
				indexOfBestThisGenerationFitness = i;
			}
		}

		if(thisGenerationFitness > maxFitness){
			maxFitness = thisGenerationFitness;
			indexOfCandidateWithBestFitness = indexOfBestThisGenerationFitness;
            iterationWhereBestIsFound = generationCounter;
		}

		if(maxFitness == numberOfClauses){
			break;
		}

		//select from breeding pool
		if(selection == "rs"){
			breedingPool = rankSelection(candidates, vectorOfClauses);
		}
		else if(selection == "ts"){
			breedingPool = tournamentSelection(candidates, vectorOfClauses);

		}
		else if(selection=="bs"){
			//pass numberOfClauses as a parameter so that we can easily go from the fitness number to the proportion / percentage
			breedingPool = boltzmannSelection(candidates, vectorOfClauses, numberOfClauses);
			//boltzmann
		}
		else{
			cout<<"ERROR SELECTION\n"<<endl;
			exit(-1);
		}

		//crossover
		vector<vector<int> > newCandidates = crossoverWrapper(crossover, numberOfIndividuals, breedingPool, crossoverProbability, mutationProbability);
		generationCounter+=1;

		candidates = newCandidates;
	}

	string bestAsString;

	for(int i = 0; i < candidates[indexOfCandidateWithBestFitness].size(); i++){
		bestAsString += to_string(candidates[indexOfCandidateWithBestFitness][i]);
	}
    /*Check this*/
    cout<<"1. filename: "<<fileName<<endl;
    cout<<"2. Number of variables: "<<numberOfVariables<<", Number of clauses: "<<numberOfClauses<<endl;
    cout<<"3. Best assignment number of clauses satisfied:"<<maxFitness<<", Best assignment percentage of clauses satisfied: "
    <<(double)maxFitness/(double)numberOfClauses<<endl;
    cout<<"4. Assignment with best fitness: "<<bestAsString<<endl;
    cout<<"5. Iteration assignment was found: "<<iterationWhereBestIsFound<<endl;
    /*Check this*/

}

vector<int> gaMutate(vector<int> child, double probabilityOfMutation){

	for(int i = 0; i < child.size(); i++){

    double randomProbability = (double)rand()/(double)RAND_MAX;

		if(randomProbability <= probabilityOfMutation){
			if(child[i] == 0){
				child[i] = 1;
			}
			else{
				child[i] = 0;
			}
		}
	}
	return child;
}

vector<vector<int> > tournamentSelection(vector< vector<int> > candidates, vector<vector<int> > vectorOfClauses){

	int counter = 0;

	vector<vector<int> > breedingPopulation;

	while(counter < candidates.size()){

		int firstCandidateIndex = rand() % candidates.size();
		int secondCandidateIndex = rand() % candidates.size();

		int firstCandidateFitness = getFitnessIntegers(vectorOfClauses, candidates[firstCandidateIndex]);
		int secondCandidateFitness = getFitnessIntegers(vectorOfClauses, candidates[secondCandidateIndex]);

		if(firstCandidateFitness > secondCandidateFitness){
			breedingPopulation.push_back(candidates[firstCandidateIndex]);
		}
		else{
			breedingPopulation.push_back(candidates[secondCandidateIndex]);
		}

		counter+=1;
	}

	return breedingPopulation;

}


vector<vector<int> > boltzmannSelection(vector< vector<int> > candidates, vector<vector<int> > vectorOfClauses, int numberOfClauses){

  double denom = 0;
  double scale = 100.0;
  //prevents this from running really slowly
  vector<double> efitness;
  //calculate denominator
  for(int i = 0; i< candidates.size(); i++){

  	//gives us fitness as a proportion

    double fitness = ( (double)getFitnessIntegers(vectorOfClauses, candidates[i]) ) / ((double)numberOfClauses);
    fitness *= scale;
   
    double eToPower = pow(EulerConstant, fitness);
    denom+= eToPower;
    efitness.push_back(eToPower);
  }

  vector<vector<int> > breedingPopulation;

  int counter = 0;

  //randomly select candidates.size() number of individuals for the vector 
  while(counter < candidates.size()) {
    int randomIndex = rand() % candidates.size();

    double randomProbability =  (double)rand()/ (double)RAND_MAX;

    double numerator = efitness[randomIndex];

    if(randomProbability <= numerator/denom){
      breedingPopulation.push_back(candidates[randomIndex]);
      counter+=1;
    }
  }

  return breedingPopulation;
}

/*
For each candidate, create a struct for it with probability of selection, fitness, and index in original candidate vector. This is so we can get the right candidate
after we sort our vector for ranking.

We assign probability to each candidate based upon the summation formula, and then we select candidates until we have filled our breeding
pool

*/
vector<vector<int> > rankSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses){

  vector<candidateFitnessAndPosition> candidateFitness;

  for(int i = 0; i < candidates.size(); i++){
    int fitness = getFitnessIntegers(vectorOfClauses, candidates[i]);
    candidateFitnessAndPosition candidateStruct;
    candidateStruct.fitnessScore = fitness;
    candidateStruct.indexInCandidateVector = i;
    candidateStruct.probabilityForSelection = 0;

    candidateFitness.push_back(candidateStruct);

  }
  //candidateFitness stores the candidateStruct info for each vector within the vector of vectors

  //yields a vector of structs where they are ordered by highest fitness score to lowest fitness score
  sort(candidateFitness.begin(), candidateFitness.end(), descending);

  double sum = 0.0f;
  for(int i =1; i <candidateFitness.size()+1; i++){
  	sum+=i;
  }

  for(int i = 1; i < candidateFitness.size()+1; i++){
  	candidateFitness[i-1].probabilityForSelection = (double) i  /sum;
  }
  vector<vector<int> > breedingPopulation;

  int counter = 0;

  while(counter < candidates.size()){
    int randomIndex = rand() % candidateFitness.size();

    double randomProbability = (double)rand()/ (double)RAND_MAX;

    if(randomProbability <= candidateFitness[randomIndex].probabilityForSelection){
      breedingPopulation.push_back(candidates[candidateFitness[randomIndex].indexInCandidateVector]);
      counter+=1;
    }
  }
  return breedingPopulation;
}

void pbil(string fileName, int numIndividuals, double plr, double nlr, double mutProb, double mutAmount, int numIterations)
{
	//reads in cnf
	int numVar;
	int numClauses;
	vector< vector<int> > cnf = readFile(fileName, &numVar, &numClauses);

	//gens probability vector
	vector<double> probVector;
	double startingValue = 0.5;
	for(int i = 0; i < numVar; i++)
	{
		probVector.push_back(startingValue);
	}

	//vector of ints used to keep track of fitness
	vector<int> fitness;

	int itNum = 0;
	bool updating = true;
    int iterationWhereBestIsFound = 0;
    int maxFitness = 0;
    int indexOfCandidateWithBestFitness = 0;

    vector<vector<bool> > cs;

	while(updating)
	{
		cs = genCS(numIndividuals, probVector);

		//gets fitness for each solution
		for(int i = 0; i < numIndividuals; i++)
		{
            int gottenFitness = getFitness(cnf, cs.at(i));
			fitness.push_back(gottenFitness);

            if(gottenFitness > maxFitness){
                maxFitness = gottenFitness;
                indexOfCandidateWithBestFitness = i;
                iterationWhereBestIsFound = itNum;
            }

		}

		//finds best and worst solution
		int bestIndex = 0;
		int worstIndex = 0;
		for(int i = 1; i < numIndividuals - 1; i = i + 2)
		{
			if(fitness.at(i) < fitness.at(i+1))
			{
				if(fitness.at(i) < fitness.at(worstIndex))
				{
					worstIndex = i;
				}
				if(fitness.at(i+1) > fitness.at(bestIndex))
				{
					bestIndex = i + 1;
				}
			}else{
				if(fitness.at(i+1) < fitness.at(worstIndex))
				{
					worstIndex = i + 1;
				}
				if(fitness.at(i) > fitness.at(bestIndex))
				{
					bestIndex = i;
				}
			}
		}

		//break if the solution is found
		if(fitness.at(bestIndex) == cnf.size())
		{
			break;
            cout << "Found solution." << endl;
		}
		//updating the probVector
		//create bestVect (where true=1)
		//create worstVect (where false=1)
		//needs to be 1s and 0s to update probVector
		vector<int> bestVect;
		vector<int> worstVect;

		for(int i = 0; i < numVar; i++)
		{
			if(cs.at(bestIndex).at(i) == true)
			{
				bestVect.push_back(1);
			}else{
				bestVect.push_back(0);
			}

			if(cs.at(worstIndex).at(i) == true)
			{
				worstVect.push_back(1);
			}else{
				worstVect.push_back(0);
			}
	
            probVector.at(i) = probVector.at(i) * (1.0 - plr) + bestVect.at(i) * plr;
            
			if(bestVect.at(i) != worstVect.at(i))
			{
				probVector.at(i) = probVector.at(i) * (1.0 - nlr) + bestVect.at(i) * nlr;
			}

		}
		probVector = mutate(probVector, mutProb, mutAmount);
		itNum++;
		if(itNum == numIterations)
		{
			updating = false;
		}
	}

    string bestAsString;

    for(int i = 0; i < cs[indexOfCandidateWithBestFitness].size(); i++){
        bestAsString += to_string(cs[indexOfCandidateWithBestFitness][i]);
    }
     /*Check this*/
    cout<<"1. filename: "<<fileName<<endl;
    cout<<"2. Number of variables: "<<numVar<<", Number of clauses: "<<numClauses<<endl;
    cout<<"3. Best assignment number of clauses satisfied:"<<maxFitness<<", Best assignment percentage of clauses satisfied: "
    <<(double)maxFitness/(double)numClauses<<endl;
    cout<<"4. Assignment with best fitness: "<<bestAsString<<endl;
    cout<<"5. Iteration assignment was found: "<<iterationWhereBestIsFound<<endl;
    /*Check this*/
}

//function that generates candidate solutions taking in N and probVector --> vector of vectors
vector< vector<bool> > genCS(int numIndividuals, vector<double> probVector)
{
	vector< vector<bool> > cs;
	for(int i = 0; i < numIndividuals; i++)
	{
		vector<bool> csToAdd;
		cs.push_back(csToAdd);
		for(int j = 0; j < probVector.size(); j++)
		{
			bool csBool;
			cs.at(i).push_back(csBool);
			double prob = (double)rand()/(double)RAND_MAX;
			if(prob <= probVector.at(j))
			{
				cs.at(i).at(j) = true;
			}else{
				cs.at(i).at(j) = false;
			}
		}
	}
	return cs;
}

//mutation function employed by pbil that updates the probability vector
vector<double> mutate(vector<double> pv, double mutProb, double mutAmount)
{
	vector<double> newPV = pv;
	for(int i = 0; i < pv.size(); i++)
	{
		double prob1 = (double)rand()/(double)RAND_MAX;
		int mutateDir;
		if(prob1 <= mutProb)
		{
			double prob2 = (double)rand()/(double)RAND_MAX;
			if(prob2 <= 0.5)
			{
				mutateDir = 1;
			}else{
				mutateDir = 0;
			}
			newPV.at(i) = pv.at(i) * (1.0 - mutAmount) + mutateDir * (mutAmount);
		}
	}
	return newPV;
}

//get fitness function employed by pbil
int getFitness(vector< vector<int> > cnf, vector<bool> cs)
{
	int indexNum;
	int fitValue = 0;

	for(int i = 0; i < cnf.size(); i++)
	{
		for(int j = 0; j < cnf.at(i).size(); j++)
		{
			indexNum = cnf.at(i).at(j);
			if(indexNum < 0)
			{
				indexNum = indexNum / -1;
			}
			indexNum--;
			if((cnf.at(i).at(j) > 0 && cs.at(indexNum)) || (cnf.at(i).at(j) < 0 && !cs.at(indexNum)))
			{
				fitValue++;
				break;
			}
		}
	}
	return fitValue;
}