#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include "utility.h"

using namespace std;

vector< vector<int> > readFile(std::string filename, int *numberOfVariables, int *numberOfClauses){


    ifstream cnfStream; // will be the CNF file
    string line;
    int numVariables;
    int numClauses;

    int counter = 0;

    vector< vector<int> > vectorOfClauses;

    cnfStream.open(filename, ios::in);

    while(cnfStream.good()){

    	while(getline(cnfStream, line)){

    		//gives the index after numVariables and before numClauses
    		int startIndex = 0;
    		vector<int> clauseVector;

    		//if the strings
    		if( strncmp(&line[0],"p", 1) == 0 ) {

    			for(int i = 6; i < line.length(); i++){
    				if( (strncmp(&line[i], " ", 1)) == 0) {
    					startIndex = i;
    					break;
    				}
    			}
    			numVariables = atoi(line.substr(5, startIndex).c_str());
    			numClauses = atoi(line.substr(startIndex+1, 100).c_str());

                *numberOfVariables = numVariables;
                *numberOfClauses = numClauses;

//    			cout<<"Got numVariables = "<< numVariables <<"\n"<<endl;
//    			cout<<"Got numClauses = "<< numClauses <<"\n"<<endl;    			
    		}

    		else if( (strncmp(&line[0], "c", 1)) != 0) {
    			string numberString = "";

    			for(int i = 0; i < line.length(); i++){
    				if( (strncmp(&line[i], " ", 1)) != 0) {
    					numberString += line[i];
    				}
    				else{
    					clauseVector.push_back(atoi(numberString.c_str()));
    					numberString = "";
    				}
    			}

	    		vectorOfClauses.push_back(clauseVector);
    		}
    		
    	}
    }
	return vectorOfClauses;
}

void printVofV(vector< vector<int> > cnf)
{
	for(int i = 0; i < cnf.size(); i++)
	{
		for(int j = 0; j < cnf.at(i).size(); j++)
		{
			cout << "Value at i=" << i << ", j=" << j << " is " << cnf.at(i).at(j) << endl;
		}
		cout<<"\n"<<endl;
	}
}

int getFitness(vector< vector<int> > cnf, bool cs[])
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
            if((cnf.at(i).at(j) > 0 && cs[indexNum]) || (cnf.at(i).at(j) < 0 && !cs[indexNum]))
            {
                fitValue++;
                break;
            }
                
        }
    }
    return fitValue;
}

int getFitnessIntegers(vector<vector<int> > CNF, vector<int> candidate){
      int indexNum;
    
    int fitValue = 0;

    for(int i = 0; i < CNF.size(); i++)
    {
        for(int j = 0; j < CNF.at(i).size(); j++)
        {
            indexNum = CNF.at(i).at(j);
            if(indexNum < 0)
            {
                indexNum = indexNum / -1;
            }
            indexNum--;
            if((CNF.at(i).at(j) > 0 && candidate[indexNum]) || (CNF.at(i).at(j) < 0 && !candidate[indexNum]))
            {
                fitValue++;
                break;
            }
                
        }
    }
    return fitValue;
}