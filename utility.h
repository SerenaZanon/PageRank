#include <iostream>
#include <vector>
#include<omp.h>
using namespace std;

//GENERIC UTILITY

bool compareFirstElement(const pair<int, int>& a, const pair<int, int>& b);

pair<double, double> greatestDifference(const vector<pair<double, double>>& v);

vector<string> retrieveAdHocFolder();

vector<string> retrieveGenericFolder();

vector<string> retrieveAllFolders();

void savingIntoCSV(const string& filePath, const vector<vector<float>>& v);

void savingPairsIntoCSV(const string& filePath, const vector<vector<pair<double, double>>>& v);

//INPUT AND DATA STRUCTURES

pair<int, vector<pair<int, int>>> readingNormalizingEdges(const string& fileName);

vector<pair<int, vector<int>>> identifyingFullRows(const vector<pair<int, int>>& rows);

vector<vector<int>> columnWiseDataStructure(const int& numNodes, const vector<pair<int, vector<int>>>& indexes);

vector<vector<float>> pageRankMatrixDataStructure(const vector<vector<int>>& indexes);

pair<int, vector <vector<int>>> columnInitialization(const string& fileName);

vector<vector<float>> matrixInitialization(const string& fileName);

vector<float> probabilityDistributionVector(int numNodes);

//COMPUTATION

void pageRankMatrix(const vector<vector<float>>& m, const vector<float>& v, vector<float>& res, const float& beta, const float& teleport, const int& numThreads, vector<double>& timeThreads);

void pageRankColumn(const vector<vector<int>>& indexes, const vector<float>& v, vector<float>& res, const float& beta, const float& teleport, const int& numThreads, vector<double>& timeThreads);