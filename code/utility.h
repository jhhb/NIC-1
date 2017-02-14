#ifndef __UTILITY_H
#define __UTILITY_H
 
 using namespace std;
// This is the content of the .h file, which is where the declarations go
vector<vector<int> > readFile(string filename); // function prototype for add.h -- don't forget the semicolon!
void printVofV(vector<vector<int> > cnf);
int getFitness(vector< vector<int> > cnf, bool cs[]);

// This is the end of the header guard
#endif