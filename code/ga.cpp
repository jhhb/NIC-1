#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "ga.h"

using namespace std;

 // takes two parent vectors and the crosspoint(random value in length of vector) 
  // and returns 1 baby. its 2 to 1 conversion so call it twice to get two chilren
  vector<int> onePointCrossover(vector<int> mom, vector<int> dad, int crossPoint) {

  	vector<int> baby;
  	for (int i = 0; i < mom.size(); i++) {
  		if (i < crossPoint) {
  			baby.push_back(mom.at(i));
  		} else {
  			baby.push_back(dad.at(i));
  		}
  	}
  	return baby;
  }

  // takes two parent vectors and randomly scrambles their genes to make 1 baby
  // call twice for two random kids from these parents
   vector<int> uniformCrossover(vector<int> mom, vector<int> dad){

   	vector<int> baby;
  	for (int i = 0; i < mom.size(); i++) {
  		int crossBool = rand() % 2;
  		if (crossBool == 0) {
  			baby.push_back(mom.at(i));
  		} else {
  			baby.push_back(dad.at(i));
  		}
  	}
  	return baby;
  }

  // takes the desired variable size (from problem file that is read in) and the population size (from command line args.)
  // and generates a population of random boolean vectors that can then be tested for fitness.
  vector< vector<int> > genSolutions(int varSize, int populationSize) {

  	vector< vector<int> > Population;
  	for (int i = 0; i < populationSize; i++) {

  		vector<int> candidate;
  		for (int j = 0; j < varSize; j++) {
  			candidate.push_back(rand() % 2);
  		}
  		Population.push_back(candidate);
  	}
  	return Population;
  }