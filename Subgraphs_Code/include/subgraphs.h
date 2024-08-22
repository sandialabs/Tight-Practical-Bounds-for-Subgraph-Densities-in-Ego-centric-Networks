//subgraphs.h

#ifndef SUBGRAPHS_H
#define SUBGRAPHS_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <ios>
#include <sstream>
#include <set>
#include <queue>
#include <stack>
#include <algorithm> //Set difference
#include <time.h> //Seed random
#include <random>
#include <unordered_map>
#include <unordered_set>

//Used to convert from non-induced counts of three vertex subgraphs to induced counts
const long A3[4][4] = { 
    {1, -1, 1, -1},
    {0, 1, -2, 3},
    {0, 0, 1, -3},
    {0, 0, 0, 1}
};

//Used to convert from non-induced counts of four vertex subgraphs to induced counts
const long A4[11][11] = { 
    {  1,  -1,   1,   1,  -1,  -1,  -1,   1,   1,  -1,   1},
    {  0,   1,  -2,  -2,   3,   3,   3,  -4,  -4,   5,  -6},
    {  0,   0,   1,   0,   0,   0,  -1,   1,   2,  -2,   3},
    {  0,   0,   0,   1,  -3,  -3,  -2,   5,   4,  -8,  12},
    {  0,   0,   0,   0,   1,   0,   0,  -1,   0,   2,  -4},
    {  0,   0,   0,   0,   0,   1,   0,  -1,   0,   2,  -4},
    {  0,   0,   0,   0,   0,   0,   1,  -2,  -4,   6, -12},
    {  0,   0,   0,   0,   0,   0,   0,   1,   0,  -4,  12},
    {  0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   3},
    {  0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -6},
    {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1}
};

std::vector< std::pair<int,int> > getEdgeList(std::ifstream &graphFile, int &n, std::vector<int> &indices);
std::vector< std::vector<double> > subgraphCounts(const std::vector<std::pair<int, int> > &edges, const int n, const std::vector<int> &indices, int cutOff = 10);
std::vector<bool> largestConnectedComponent(const std::vector< std::pair<double, double> > &coordinates, double d);
std::vector< std::vector< std::pair<double, double> > > prune(const std::vector< std::vector<double> > &coordinates, std::vector<double> &d, const double c = 0.95, const double eps = 10E-6);
double area(const std::vector< std::pair<double, double> >  &coordinates, const double d, double &minX, double &maxX, const int numVertices, const int PFAMindex, const int numIt = 10E6);
std::vector<double> getVals(const std::vector< std::vector< std::pair<double, double> > > &prunedCoordinates, const std::vector<double> & d, const int numVertices, const int n, const int numIt = 10E5);
void mainFn(std::string fileName, std::string fileExtension);

#endif