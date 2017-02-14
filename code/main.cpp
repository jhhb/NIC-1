
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "utility.h"

using namespace std;

int main( int argc, const char* argv[] )
{

	// if(argc < 10){
	// 	cout <<"\nusage: ./a.out filename.spn.cnf individuals selection crossover crossover_probability mutation_probability iterations algorithm\n"<<endl;
	// }

	//read parameters into proper vector 
	


	// vector<string> stringParameterVector;
	// vector<double> doubleParameterVector;

	// for(int i = 2; i <10; i++){
	// 	if(i == 2 || i == 9 || i == 4 || i == 5){
	// 		stringParameterVector.push(argv[i]);
	// 	}
	// 	else{
	// 		doubleParameterVector.push(argv[i]);
	// 	}
	// }

	vector<vector<int> > vectorOfClauses = readFile("/Users/jamesboyle/Desktop/NIC-lab1/project1-ga-pbil-for-maxsat 2/maxsat-problems/maxsat-random/highgirth/3SAT/HG-3SAT-V250-C1000-1.cnf");

	printVofV(vectorOfClauses);

/*
Justin's tests for fitness function

*/


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



	//read file into memory
}
