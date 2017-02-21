
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

int getFitness(vector< vector<int> > cnf, vector<bool> cs);
void printVofV(vector< vector<int> > cnf);
void pbil(string fileName, int numIndividuals, double plr, double nlr, double mutProb, double mutAmount, int numIterations);
vector< vector<bool> > genCS(int numIndividuals, vector<double> probVector);
vector<double> mutate(vector<double> pv, double mutProb, double mutAmount);

bool descending (candidateFitnessAndPosition i, candidateFitnessAndPosition j) { return (i.fitnessScore > j.fitnessScore); }



int main( int argc, const char* argv[] )
{
	//format of a GA run

	runGA("/Users/mbernard/Computer Science/Nature_Inspired/project1-ga-pbil-for-maxsat/maxsat-problems/maxsat-crafted/MAXCUT/SPINGLASS/t3pm3-5555.spn.cnf",
		100, "rs", "uc", 0.7, 0.01, 20);

	pbil("/Users/mbernard/Computer Science/Nature_Inspired/project1-ga-pbil-for-maxsat/maxsat-problems/maxsat-crafted/MAXCUT/SPINGLASS/t3pm3-5555.spn.cnf",
		100, 0.01, 0.01, 0.01, 0.05, 20);
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
				if(crossoverProbability <= randomProbability){

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
		/*Need to check if we are doing 1c, uc correctly, if uc needs a probability associated with it or a probability for one parent */
		vector<vector<int> > newCandidates = crossoverWrapper(crossover, numberOfIndividuals, breedingPool, crossoverProbability, mutationProbability);
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


/*
We arbitrarily pick a k and an mValue and use them for tournament selection.

We randomly select mValue candidates, get their fitness, save the associated fitness and index into a candidateStruct, and push it back into a candidate struct

Outside the for loop, we sort these, and then take the k with the highest fitness Scores
*/
vector<vector<int> > tournamentSelection(vector< vector<int> > candidates, vector<vector<int> > vectorOfClauses){
/*
Ask majercik about k and mValue 
*/


/*
HOW DO WE PICK THESE??
*/



  int k = 1; // k = 1;

  int mValue = 2; // mValue = 2;

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

    if(numerator/denom <= randomProbability){
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

/* Ostensibly fixes the sum issue */
  double sum;
  for(int i =1; i <candidateFitness.size()+1; i++){
  	sum+=i;
  }

  for(int i = 1; i < candidateFitness.size()+1; i++){
  	candidateFitness[i-1].probabilityForSelection = (double) i  /sum;
  }
  /*Fixes sum */

  vector<vector<int> > breedingPopulation;

  int counter = 0;

  while(counter < candidates.size()){
    int randomIndex = rand() % candidateFitness.size();

    double randomProbability = (rand() % 101) / 100.0f;

    if(candidateFitness[randomIndex].probabilityForSelection <= randomProbability){
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

	double firstFitness = 0;

	while(updating)
	{
		vector< vector<bool> > cs = genCS(numIndividuals, probVector);

		//gets fitness for each solution
		for(int i = 0; i < numIndividuals; i++)
		{
			fitness.push_back(getFitness(cnf, cs.at(i)));
		}

		if(itNum == 0){

			for(int i =0; i<fitness.size(); i++){
				if(fitness[i] > firstFitness){
					firstFitness = fitness[i];
				}
			}
			cout<<"First fitness from PBIL is: "<<firstFitness<<endl;
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
			//add cout
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
				worstVect.push_back(0);
			}else{
				worstVect.push_back(1);
			}

				//MOVE PLR UP
			if(bestVect.at(i) != worstVect.at(i))
			{
				probVector.at(i) = probVector.at(i) * (1.0 - plr) + bestVect.at(i) * plr;
				probVector.at(i) = probVector.at(i) * (1.0 - nlr) + worstVect.at(i) * nlr;
			}

		}

		probVector = mutate(probVector, mutProb, mutAmount);
		itNum++;
		if(itNum == numIterations)
		{
			updating = false;
		}
	}

	double maxFitness = 0;

	for(int i =0; i<fitness.size(); i++){
		if(fitness[i] > maxFitness){
			maxFitness = fitness[i];
		}
	}
	cout<<"Max fitness from PBIL is: "<<maxFitness<<endl;
	//couts
}

//function that generates candidate solutions taking in N and probVector --> vector of vectors
//TESTED: GOOD
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