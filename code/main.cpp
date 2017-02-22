
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
int runGA(string fileName, int numberOfIndividuals, string selection, string crossover, double crossoverProbability, double mutationProbability, int numberOfGenerations);

int getFitness(vector< vector<int> > cnf, vector<bool> cs);
void printVofV(vector< vector<int> > cnf);
int pbil(string fileName, int numIndividuals, double plr, double nlr, double mutProb, double mutAmount, int numIterations);
vector< vector<bool> > genCS(int numIndividuals, vector<double> probVector);
vector<double> mutate(vector<double> pv, double mutProb, double mutAmount);

bool descending (candidateFitnessAndPosition i, candidateFitnessAndPosition j) { return (i.fitnessScore < j.fitnessScore); }

void gaTester(string fileName);
void pbilTester(string fileName);

int main( int argc, const char* argv[] )
{
	//format of a GA run
	string filePath = "/Users/jboyle/Desktop/NIC/NIC-1/project1-ga-pbil-for-maxsat 2/maxsat-problems/maxsat-crafted/MAXCUT/DIMACS_MOD/";


	//string fileName = "t7pm3-9999.spn.cnf";
    string fileName = "San1000.clq.cnf";
	//"/Users/jamesboyle/Desktop/NIC-1/NIC-1/project1-ga-pbil-for-maxsat 2/maxsat-problems/maxsat-crafted/MAXCUT/SPINGLASS/t3pm3-5555.spn.cnf";
	
	string fullPath = filePath + fileName;
	int numberOfIndividuals = 100;
	string selection = "ts";
	string crossover = "1c";
	double crossoverProbability = 0.7;
	double mutationProbability = 0.01;
	double generations = 300;
  //  cout << pbil(fileName, 100, 0.1, 0.075, 0.02, 0.05, 100) << endl;
    gaTester(fullPath);
   //pbilTester(fullPath);


}

void pbilTester(string fileName){
    cout<<fileName<<endl;
    int numVar;
    int numClauses;
    vector< vector<int> > cnf = readFile(fileName, &numVar, &numClauses);
    cout << "numVar=" << numVar << " and numClauses=" << numClauses << endl;
    
    
    int iterationsTested = 10;
    double workingSum = 0.0;
    double averagePerf = 0.0;
    
    //starting values
    int numberOfIndividuals     =100;

    double plr                  =0.05;
    double nlr                  =0.075;
    double mutProb              =0.02;
    double mutAmount            =0.1;

    int numberOfIterations      =500;


    for(int g = 0; g < iterationsTested; g++){

        workingSum = workingSum + pbil(fileName, numberOfIndividuals, plr, nlr, mutProb, mutAmount, numberOfIterations);
        cout<<g<<endl;
    }
    
    averagePerf = (double)workingSum/(double)iterationsTested;
    workingSum = 0.0;
    
    cout << "NumInd=" << numberOfIndividuals << "  plr=" << plr << "  nlr=" << nlr << "  mProb=" << mutProb << "  mAmount=" <<
    mutAmount << "  NumIt=" << numberOfIterations << "  average fitness=" << averagePerf << endl;

    // for(int a = 0; a <= 1; a++){
    //     //mod numIndividuals
    //     if(a==0){numberOfIndividuals=100;}
    //     if(a==1){numberOfIndividuals=200;}
    //     for(int b = 0; b <= 2; b++){
    //         //mod plr
    //         if(b==0){plr=0.1;}
    //         if(b==1){plr=0.05;}
    //         if(b==2){plr=0.2;}
    //         for(int c = 0; c <= 2; c++){
    //             //mod nlr
    //             if(c==0){nlr=0;}
    //             if(c==1){nlr=0.075;}
    //             if(c==2){nlr=0.15;}
    //             for(int d = 0; d <= 1; d++){
    //                 //mod mutProb
    //                 if(d==0){mutProb=0.01;}
    //                 if(d==1){mutProb=0.02;}
    //                 if(d==2){mutProb=0.04;}
    //                 for(int e = 0; e <= 1; e++){
    //                     //mod mutProb
    //                     if(e==0){mutAmount=.05;}
    //                     if(e==1){mutAmount=0.1;}
    //                     for(int f = -1; f <=1; f++){
    //                         //mod numGen
    //                         if(f==-1){numberOfIterations = 200;}
    //                         if(f==0){numberOfIterations = 1000;}
    //                         if(f==1){numberOfIterations = 2000;}
    //                         for(int g = 0; g < iterationsTested; g++){
    //                             workingSum = workingSum + pbil(fileName, numberOfIndividuals, plr, nlr, mutProb, mutAmount, numberOfIterations);
    //                         }
    //                         averagePerf = (double)workingSum/(double)iterationsTested;
    //                         workingSum = 0.0;
    //                         cout << "NumInd=" << numberOfIndividuals << "  plr=" << plr << "  nlr=" << nlr << "  mProb=" << mutProb << "  mAmount=" << mutAmount << "  NumIt=" << numberOfIterations << "  average fitness=" << averagePerf << endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
}

void gaTester(string fileName){
    cout<<"FILENAME: "<<fileName<<endl;
    int numVar;
    int numClauses;
    vector< vector<int> > cnf = readFile(fileName, &numVar, &numClauses);
    cout << "numVar=" << numVar << " and numClauses=" << numClauses << endl;
    
    
    int iterationsTested = 10;
    double workingSum = 0.0;
    double averagePerf = 0.0;
    
    //starting values
    int numberOfIndividuals             =100;
    string selection                    ="bs";
    string crossover                    ="1c";
    double crossoverProbability         =0.7;
    double mutationProbability          =0.01;
    int numberOfGenerations             =500;

    for(int g = 0; g < iterationsTested; g++){
        workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
        cout<<g<<endl;
    }
    
    averagePerf = (double)workingSum/(double)iterationsTested;
    workingSum = 0.0;

    cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << 
    crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" 
    << averagePerf << endl;
    
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    
    // //adjusting one var at a time
    // //number of individuals
    // numberOfIndividuals = 200;
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    // numberOfIndividuals = 100;
    
    // //selection type
    // selection = "rs";
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    // selection = "bs";
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    // selection = "ts";
    
    // //crossover
    // crossover = "uc";
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    // crossover = "1c";
    
    // //crossover probability
    // crossoverProbability = 0.9;
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    // crossoverProbability = 0.7;
    
    // //mutation probability
    // mutationProbability = .05;
    // for(int i = 0; i < iterationsTested; i++){
    //     workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    // }
    // averagePerf = (double)workingSum/(double)iterationsTested;
    // workingSum = 0.0;
    // cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    // mutationProbability = 0.01;
    

    
    // for(int a = 0; a <= 1; a++){
    //     //mod numIndividuals
    //     if(a==0){numberOfIndividuals=100;}
    //     if(a==1){numberOfIndividuals=200;}
    //     for(int b = 0; b <= 2; b++){
    //         //mod selection
    //         if(b==0){selection="ts";}
    //         if(b==1){selection="rs";}
    //         if(b==2){selection="bs";}
    //         for(int c = 0; c <= 1; c++){
    //             //mod crossover
    //             if(c==0){crossover="1c";}
    //             if(c==1){crossover="uc";}
    //             for(int d = 0; d <= 1; d++){
    //                 //mod crossover prob
    //                 if(d==0){crossoverProbability=.7;}
    //                 if(d==1){crossoverProbability=.9;}
    //                 for(int e = 0; e <= 1; e++){
    //                     //mod mutProb
    //                     if(e==0){mutationProbability=.01;}
    //                     if(e==1){mutationProbability=.1;}
    //                     for(int f = -1; f <=1; f++){
    //                         //mod numGen
    //                         if(f==-1){numberOfGenerations = 200;}
    //                         if(f==0){numberOfGenerations = 400;}
    //                         if(f==1){numberOfGenerations = 600;}
    //                         for(int g = 0; g < iterationsTested; g++){
    //                             workingSum = workingSum + runGA(fileName, numberOfIndividuals, selection, crossover, crossoverProbability, mutationProbability, numberOfGenerations);
    //                         }
    //                         averagePerf = (double)workingSum/(double)iterationsTested;
    //                         workingSum = 0.0;
    //                         cout << "NumInd=" << numberOfIndividuals << "  selection=" << selection << "  crossover=" << crossover << "  cProb=" << crossoverProbability << "  mProb=" << mutationProbability << "  NumGen=" << numberOfGenerations << "  average fitness=" << averagePerf << endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
     
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

int runGA(string fileName, int numberOfIndividuals, string selection, string crossover, double crossoverProbability, double mutationProbability, int numberOfGenerations){

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

	/* Done just to see what our "best" is at the start so we can compare at end */

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
		}

		if(maxFitness == numberOfClauses){
			break;
			//time to print
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
    return maxFitness;
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

int pbil(string fileName, int numIndividuals, double plr, double nlr, double mutProb, double mutAmount, int numIterations)
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

	double maxFitness = 0;

	for(int i =0; i<fitness.size(); i++){
		if(fitness[i] > maxFitness){
			maxFitness = fitness[i];
		}
	}
	return maxFitness;
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