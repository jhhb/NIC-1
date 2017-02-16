#ifndef __GA_H
#define __GA_H
 
 using namespace std;
// This is the content of the .h file, which is where the declarations go
vector<int> onePointCrossover(vector<int> mom, vector<int> dad, int crossPoint);
vector<int> uniformCrossover(vector<int> mom, vector<int> dad);
vector< vector<int> > genSolutions(int varSize, int populationSize);

struct candidateFitnessAndPosition {
  int indexInCandidateVector;
  int fitnessScore;
  double probabilityForSelection;
} ;

// This is the end of the header guard
#endif