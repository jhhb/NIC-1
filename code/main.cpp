
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

void martinTests();
void justinTests();
const double EulerConstant = std::exp(1.0);

vector<vector<int> > rankSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses);
vector<vector<int> > tournamentSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses);
vector<vector<int> > boltzmannSelection( vector< vector<int> > candidates, vector< vector<int> > vectorOfClauses, int numberOfClauses);
vector<int> gaMutate(vector<int> child, double probabilityOfMutation);
void runGA(string fileName, int numberOfIndividuals, string selection, string crossover, double crossoverProbability, double mutationProbability, int numberOfGenerations);

bool descending (candidateFitnessAndPosition i, candidateFitnessAndPosition j) { return (i.fitnessScore > j.fitnessScore); }



int main( int argc, const char* argv[] )
{
	//format of a GA run

	runGA("/Users/jamesboyle/Desktop/NIC-1/NIC-1/project1-ga-pbil-for-maxsat 2/maxsat-problems/maxsat-crafted/MAXCUT/SPINGLASS/t3pm3-5555.spn.cnf",
		100, "ts", "1c", 0.7, 0.05, 100);

}


/*While we dont have all the individuals we want, we select two random parents from the breeding pool, cross them using a random position in the array, and add to
new array.

We then return the newCandidates

*/
vector<vector<int> > crossoverWrapper(string crossover, int numberOfIndividuals, vector<vector<int> > breedingPool, double mutationProbability){
	vector<vector<int> > newCandidates;

	//1point crossover
	if(crossover == "1c"){

			int offspringCounter = 0;

			//generate numberOfIndividuals offspring from the breeding pool
			while(offspringCounter < numberOfIndividuals){
				int parent1Index = rand() % breedingPool.size();
				int parent2Index = rand() % breedingPool.size();
				
				//Is this right???
				vector<int> child = onePointCrossover(breedingPool[parent1Index], breedingPool[parent2Index], rand() % breedingPool[parent1Index].size() );

				gaMutate(child, mutationProbability);

				newCandidates.push_back(child);

				offspringCounter+=1;
			}

		}

		//Uniform Crossover
		else{
			int offspringCounter = 0;

			//generate numberOfIndividuals offspring from the breeding pool
			while(offspringCounter < numberOfIndividuals){
				int parent1Index = rand() % breedingPool.size();
				int parent2Index = rand() % breedingPool.size();
				vector<int> child = uniformCrossover(breedingPool[parent1Index], breedingPool[parent2Index]);

				gaMutate(child, mutationProbability);

				newCandidates.push_back(child);

				offspringCounter+=1;
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
	double max = 0;
	double firstBestFitness;
	/* Done just to see what our "best" is at the start so we can compare at end */
	for(int i = 0; i < candidates.size();i++){
		if(getFitnessIntegers(vectorOfClauses, candidates[i]) > max){
			max = getFitnessIntegers(vectorOfClauses, candidates[i]);
		}

	}

	cout<<"First best fitness: "<<endl<<max<<"\n"<<endl;

	firstBestFitness= max;

	while(generationCounter < numberOfGenerations){
		cout<<generationCounter<<"\n"<<endl;
		//select from breeding pool
		if(selection == "rs"){
			breedingPool = rankSelection(candidates, vectorOfClauses);
		}
		else if(selection == "ts"){
			breedingPool = tournamentSelection(candidates, vectorOfClauses);

		}
		else{
			//pass numberOfClauses as a parameter so that we can easily go from the fitness number to the proportion / percentage
			breedingPool = boltzmannSelection(candidates, vectorOfClauses, numberOfClauses);
			//boltzmann
		}

		//crossover
		/*Need to check if we are doing 1c, uc correctly, if uc needs a probability associated with it or a probability for one parent */
		vector<vector<int> > newCandidates = crossoverWrapper(crossover, numberOfIndividuals, breedingPool, mutationProbability);

		generationCounter+=1;

		candidates = newCandidates;
	}

	for(int i = 0; i < candidates.size();i++){
		if(getFitnessIntegers(vectorOfClauses, candidates[i]) > max){
			max = getFitnessIntegers(vectorOfClauses, candidates[i]);
		}

	}
	cout<<"Filename: "<<fileName<<"\n"<<endl;
	cout<<numberOfVariables<<" variables and "<<numberOfClauses<<" clauses\n"<<endl;
	cout<<"Best assignment satisfies "<<numberOfClauses<<" clauses and "<<max / numberOfClauses<<"%"<<" of clauses\n"<<endl;

	cout<<"First best fitness: "<<endl<<firstBestFitness<<"\n"<<endl;

	cout<<"Final best fitness: "<<endl<<max<<"\n"<<endl;

}

vector<int> gaMutate(vector<int> child, double probabilityOfMutation){

	for(int i = 0; i < child.size(); i++){

	//IS THIS RIGHT???

    double randomProbability = (rand() % 101) / 100.0f;

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


/*
We arbitrarily pick a k and an mValue and use them for tournament selection.

We randomly select mValue candidates, get their fitness, save the associated fitness and index into a candidateStruct, and push it back into a candidate struct

Outside the for loop, we sort these, and then take the k with the highest fitness Scores
*/
vector<vector<int> > tournamentSelection(vector< vector<int> > candidates, vector<vector<int> > vectorOfClauses){

/*
HOW DO WE PICK THESE??
*/
  int k = 20;

  int mValue = 50;

  vector<candidateFitnessAndPosition> candidateFitness;

  for(int i = 0; i < mValue; i++){

  	//For candidates of size N, this gives us range of 0 to N-1 indexes
    int randomIndex = rand() % candidates.size();

    int fitness = getFitnessIntegers(vectorOfClauses, candidates[randomIndex]);

    candidateFitnessAndPosition candidateStruct;
    candidateStruct.fitnessScore = fitness;

    candidateStruct.indexInCandidateVector = randomIndex;

    candidateStruct.probabilityForSelection = 0;
    candidateFitness.push_back(candidateStruct);
  }

  sort(candidateFitness.begin(), candidateFitness.end(), descending);

  vector<vector<int> > breedingPopulation;

  for(int i = 0; i < k; i ++){
    breedingPopulation.push_back(candidates[candidateFitness[i].indexInCandidateVector]);
  }

  return breedingPopulation;

}


/* Why the fuck is this so slow? And how can we speed it up? 

First calculate the denominator of the probability calculation, which is just the some of e^(fitness) for each individual.

*/
vector<vector<int> > boltzmannSelection(vector< vector<int> > candidates, vector<vector<int> > vectorOfClauses, int numberOfClauses){


  double denom = 0;

  //calculate denominator
  for(int i = 0; i< candidates.size(); i++){
    double fitness = getFitnessIntegers(vectorOfClauses, candidates[i]) / (double)numberOfClauses;

    double eToPower = pow(EulerConstant, fitness);
    denom+= eToPower;
  }

  vector<vector<int> > breedingPopulation;

  int counter = 0;

  //randomly select candidates.size() number of individuals for the vector 
  while(counter < candidates.size()){
    int randomIndex = rand() % candidates.size();

    //is this right???
    //THIS IS A BUG. NUMERATOR / DENOM is almost never greater than random probability, which is why this shit takes so long 
    double randomProbability =  (rand() % 101) / 100.0f;

    double fitness = getFitnessIntegers(vectorOfClauses, candidates[randomIndex]) / (double)numberOfClauses;

    double numerator = pow(EulerConstant, fitness);
  //  cout<<numerator/denom<<endl;

    if(numerator/denom >= randomProbability){
      breedingPopulation.push_back(candidates[randomIndex]);
      counter+=1;
    }
  }

  return breedingPopulation;

}

/*

For each candidate, create a struct for it with probability of selection, fitness, and index in original candidate vector. This is so we can get the right candidate
after we sort our vector for ranking.

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

  //calculate sum from i = 1 to N
  int iValueInSum = 1;
  double sum = iValueInSum;

  for(int i = 0; i < candidateFitness.size(); i++){
    candidateFitness[i].probabilityForSelection = iValueInSum / sum;
    iValueInSum +=1;
    sum += iValueInSum;
  }

  vector<vector<int> > breedingPopulation;

  int counter = 0;

  while(counter < candidates.size()){
    int randomIndex = rand() % candidateFitness.size();

    double randomProbability = (rand() % 101) / 100.0f;

    if(candidateFitness[randomIndex].probabilityForSelection >= randomProbability){
      breedingPopulation.push_back(candidates[candidateFitness[randomIndex].indexInCandidateVector]);
      counter+=1;
    }

  }

  return breedingPopulation;

}

void martinTests(){

	// need this so that the rand() function works. the beginning of main() is the only place it ever needs to go
  	srand(time(0));

  	// vectors to test crossover on
  	vector<int> mom;
  	vector<int> dad;

  	// length of the solutions used throughout testing. each solution vector has 25 elements (this comes from the problem file)
  	int variables = 25;
  	// size of the desired population for each generation to have (this comes from the command line argument)
  	int genPopulation = 20;

  	// just puting all 0s in mom and all 1s in dad
  	for (int i = 0; i < variables; i++) {

  		mom.push_back(0);
  		dad.push_back(1);
  	}

  	// crosspoint for 1 point c-over is calculated outside the function so that the children each have the genes that the other doesn't
  	// just make sure you switch the order of the two vectors in the two calls (see mom, dad; then dad, mom below)
  	int crossPoint = rand() % variables; 
  	vector<int> baby1 = onePointCrossover(mom, dad, crossPoint);
  	vector<int> baby2 = onePointCrossover(dad, mom, crossPoint);

  	// children in uniform crossover are random each time. it was easier to do it this way and if it gets brought up we can say its
  	// more like how it is in nature.
  	vector<int> baby3 = uniformCrossover(mom, dad);
  	vector<int> baby4 = uniformCrossover(dad, mom);


  	// printing.......

  	cout << "Two Single Point Crossover Babies:" << endl;
  	for (int i = 0; i < variables; i++) {
  		cout << baby1.at(i);
  	}
  	cout << endl;
  	for (int i = 0; i < variables; i++) {
  		cout << baby2.at(i);
  	}
  	cout << endl;
  	cout << "Two Uniform Crossover Babies:" << endl;
  	for (int i = 0; i < variables; i++) {
  		cout << baby3.at(i);
  	}
  	cout << endl;
  	for (int i = 0; i < variables; i++) {
  		cout << baby4.at(i);
  	}
  	cout << endl;
	
	// create population of solutions
  	vector< vector<int> > test;

  	//have vector of vector of ints 
	test = genSolutions(variables, genPopulation);

	// print those solutions
	cout << "Original Generation:" << endl;
	for (int i = 0; i < genPopulation; i++) {
		cout << "Candidate #" << i << ": ";
		for (int j = 0; j < variables; j++) {
  			cout << test.at(i).at(j);
  		}
  		cout << endl;
	}

	//have a generation of candidate solutions
	//

}

void justinTests(){

	vector< vector<int> > testCNF;
	bool testCS[] = {false, false, false, false};

	vector<int> testVector;

	testVector.push_back(1);
	testVector.push_back(2);
	testCNF.push_back(testVector);
	testVector.at(0) = 3;
	testVector.at(1) = 4;
	testCNF.push_back(testVector);
	testVector.at(0) = -1;
	testVector.at(1) = -2;
	testCNF.push_back(testVector);
	testVector.at(0) = -3;
	testVector.at(1) = -4;
	testCNF.push_back(testVector);

	int z;
	z = getFitness(testCNF, testCS);

	cout << "4 variables; 4 clauses. Fitness is " <<  z << endl;
	printVofV(testCNF);

}