#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include <algorithm>
#include<cmath>
#include<map>
#include<omp.h>
#include "utility.h"

using namespace std;

//GENERIC UTILITY

/**
 * @brief Comparison between the first elements of the two given pairs.
 *
 * @param param1 The pair a composed by two integers.
 * @param param2 The pair b composed by two integers.
 * @return A boolean value: true when the first value of the first pair is smaller than the first value of the second pair,
 *                        false otherwise.
 */
bool compareFirstElement(const pair<int, int>& a, const pair<int, int>& b) {
    return a.first < b.first;
}

/**
 * @brief Comparing pairs.
 *
 * This function individualize the pair in which the gap between the second and the first element is bigger.
 *
 * @param param1 A vector of pairs of double.
 * @return A pair of double.
 */
pair<double, double> greatestDifference(const vector<pair<double, double>>& v){

    double currentDiff = 0.0;
    pair<double, double> final;
    for(pair<double, double> couple: v){
        double diff = couple.second - couple.first;
        if(currentDiff - diff){
            currentDiff = diff;
            final = couple;
        }
    }
    return final;
}

/**
 * @brief Retrieving the paths of the datasets contained in 'adHoc' folder.
 *
 * @return A vector of strings containing the paths of the datasets.
 */
vector<string> retrieveAdHocFolder(){

    vector<string> datasets;
    datasets.push_back("datasets/adHoc/adHocDataset_025.txt");
    datasets.push_back("datasets/adHoc/adHocDataset_05.txt");
    datasets.push_back("datasets/adHoc/adHocDataset_075.txt");
    datasets.push_back("datasets/adHoc/adHocDataset_1.txt");

    return datasets;
}

/**
 * @brief Retrieving the paths of the datasets contained in 'generic' folder.
 *
 * @return A vector of strings containing the paths of the datasets.
 */
vector<string> retrieveGenericFolder(){

    vector<string> datasets;
    datasets.push_back("datasets/generic/smallDataset_sparse.txt");
    datasets.push_back("datasets/generic/smallDataset_dense.txt");
    datasets.push_back("datasets/generic/mediumDataset_sparse.txt");
    datasets.push_back("datasets/generic/mediumDataset_dense.txt");
    datasets.push_back("datasets/generic/bigDataset_sparse.txt");
    datasets.push_back("datasets/generic/bigDataset_dense.txt");

    return datasets;
}

/**
 * @brief Retrieving the paths of the datasets contained in both 'adHoc' and 'generic' folders.
 *
 * @return A vector of strings containing the paths of the datasets of both folders.
 */
vector<string> retrieveAllFolders(){

    vector<string> datasets;
    vector<string> smallDatasets = retrieveAdHocFolder();
    vector<string> genericDatasets = retrieveGenericFolder();
    datasets.reserve(smallDatasets.size() + genericDatasets.size()); // preallocate memory
    datasets.insert(datasets.end(), smallDatasets.begin(), smallDatasets.end());
    datasets.insert(datasets.end(), genericDatasets.begin(), genericDatasets.end());

    return datasets;
}

/**
 * @brief Saving information inside a csv file.
 *
 * This function allows to store the given vector of vectors of floats into a csv file.
 * The file is created and if it is open then each internal vector is written in one distinct row.
 * Each element of the internal vector is separated by a comma ','.
 * At the end, the file is closed.
 * If the file is not open then an error shows up.
 *
 * @param param1 The path of the file in which is going to be written the information.
 * @param param2 The vector of vectors of float to be stored.
 */
void savingIntoCSV(const string& filePath, const vector<vector<float>>& v) {

    ofstream file(filePath);

    if (file.is_open()) {
        for (const auto& value : v) {
            for (size_t i = 0; i < value.size(); ++i) {
                file << value[i];
                if (i != value.size() - 1) {
                    file << ",";
                }
                else {
                    file << "\n";
                }
            }
        }
        file.close();
    } else {
        cerr << "Impossible to open the file to save data" << endl;
    }
}

/**
 * @brief Saving information inside a csv file.
 *
 * This function allows to store the given vector of vectors of pairs of two double into a csv file.
 * The file is created and if it is open then each pair is written in one distinct row.
 * Each element of the pair is separated by a comma ','.
 * At the end, the file is closed.
 * If the file is not open then an error shows up.
 *
 * @param param1 The path of the file in which is going to be written the information.
 * @param param2 The vector of vectors of pairs of two double to be stored.
 */
void savingPairsIntoCSV(const string& filePath, const vector<vector<pair<double, double>>>& v) {

    ofstream file(filePath);

    if (file.is_open()) {
        for (const auto& innerVector : v) {
            for (const auto& pair : innerVector) {
                file << pair.first << "," << pair.second << ",";
            }
            file << "\n";
        }
        file.close();
    } else {
        cerr << "Impossible to open the file to save data" << endl;
    }
}

//INPUT AND DATA STRUCTURES


/**
 * @brief Reading and normalizing the ids of the nodes of the given graph.
 *
 * This function tries to open the file at the given path and if it is successful it reads each line of the file.
 * The lines that are not representing an edge, so that are not composed by two integers
 * (the ids of the nodes characterizing the edge) are ignored.
 * Once it finds a line with two integers a mapping is done:
 *      every distinct id encountered in the file is associated with a new id.
 * This is made in order to normalize the names of the nodes and
 * to manage in a better way the further computation.
 * Each pair of new ids representing the just read edge is appended to the vector 'edges'.
 * Finally, the file is closed and the vector 'edges' is ordered by the first value of each pair.
 *
 * @param param1 The path of the file in which are contained the edges.
 * @return Pair composed by an integer which represents the total number of nodes and a vector of pairs of integers.
 * Each of these pairs represents an edge in the graph where
 * - the first value is the id of the node from which the link starts,
 * - the second is the id of the node to which the link arrives.
 */
pair<int, vector<pair<int, int>>> readingNormalizingEdges(const string& filePath) {

    ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        cerr << "Impossible to open the file: " << filePath << endl;
        exit(1);
    }

    map<int, int> renaming;
    vector<pair<int, int>> edges;
    int currentID = 0, startingNode, finalNode;
    string line;

    while (getline(inputFile, line)) {
        stringstream edge(line);
        edge >> startingNode >> finalNode;
        if (edge.fail()) {
            continue;
        }

        if (renaming.find(startingNode) == renaming.end()) {
            renaming[startingNode] = currentID;
            currentID++;
        }

        if (renaming.find(finalNode) == renaming.end()) {
            renaming[finalNode] = currentID;
            currentID++;
        }

        edges.push_back(make_pair(renaming[startingNode], renaming[finalNode]));
    }

    inputFile.close();
    if(edges.empty()){
        throw runtime_error("The file isn't containing any edge to process!");
    }
    sort(edges.begin(), edges.end(), compareFirstElement);
    return make_pair(currentID, edges);

}

/**
 * @brief Identifying the outgoing links for each node.
 *
 * The goal of this function is to group together all the ids of the destination nodes
 * for all nodes that have outgoing links.
 * This is achieved by saving for each node that has at least one outgoing link,
 * a vector of ints containing the ids of the nodes that are reached from it.
 *
 * @param param1 A pair composed by an integer which represents the total number of nodes and a vector of pairs of integers.
 * Each of these pairs represents an edge in the graph where
 * - the first value is the id of the node from which the link starts,
 * - the second is the id of the node to which the link arrives.
 * @return A vector of pairs where:
 * - the first value is the id of the node from which the link starts,
 * - the second is a vector containing all the ids of the nodes reached by the node
 * represented by the first value of the pair.
 * Let's notice that this function (for the purpose of the assignment) will always be called
 * after 'readingNormalizingEdges' function, so the elements of the resulting vector are ordered by the first value of the pair.
 */
vector<pair<int, vector<int>>> identifyingFullRows(const vector<pair<int, int>>& rows) {

    vector<pair<int, vector<int>>> indexes;
    for(pair<int, int> row : rows){
        if(indexes.empty() || indexes.at(indexes.size() - 1).first != row.first){
            pair<int, vector<int>> newPair;
            newPair.first = row.first;
            indexes.push_back(newPair);
        }
        indexes.at(indexes.size() - 1).second.push_back(row.second);
    }
    return indexes;
}

/**
 * @brief Creating the column wise data structure.
 *
 * The core of the function it's the for loop over the total number of nodes.
 * An if checks if the first value of the current pair is equal to 'numProgressive':
 * - If the answer is positive than the associated vector is appended to the final vector;
 * - Otherwise an empty vector is appended.
 * Just to remind, if a node does not have an associated vector it means that it is a dead end, so it doesn't link to anybody.
 *
 * @param param1 An integer representing the total number of nodes of the graph.
 * @param param2 A vector of pairs where:
 * - the first value is the id of the node from which the link starts,
 * - the second is a vector containing all the ids of the nodes reached by the node
 * represented by the first value of the pair.
 * @return A vector of vectors of ints where:
 * - The i-th element of the external vector corresponds to the i-th node;
 * - The i-th internal vector contains the ids of the nodes reached by the i-th node, if any.
 */
vector<vector<int>> columnWiseDataStructure(const int& numNodes, const vector<pair<int, vector<int>>>& indexes){

    vector<vector<int>> rowIndexes;
    int k = 0;
    for(int numProgressive = 0; numProgressive < numNodes; numProgressive++){
        if(k < indexes.size() && numProgressive == indexes.at(k).first){
            rowIndexes.push_back(indexes.at(k).second);
            k++;
        }
        else{
            vector<int> emptyVector;
            rowIndexes.push_back(emptyVector);
        }
    }
    return rowIndexes;
}

/**
 * @brief Creating the Page Rank Matrix data structure.
 *
 * The initial point is a matrix composed by only 0s whose dimension is equal to the size of the given parameter 'indexes',
 * which means the number of nodes of the input graph.
 * For all the columns of the matrix, a check is done on the corresponding element of the input parameter:
 * - If i-th internal vector of 'indexes' is NOT empty then the all the rows identified by the integers presented in this internal vector
 * assume a value equal to 1 over the size of the i-th internal vector;
 * - Otherwise, all the rows will assume a value equal to 1 over the size of the matrix, so the number of the nodes.
 * Finally, the matrix is returned.
 *
 * @param param1 A vector of vectors of ints where:
 * - The i-th element of the external vector corresponds to the i-th node;
 * - The i-th internal vector contains the ids of the nodes reached by the i-th node, if any.
 * @return A vector of vectors of floats representing the Page Rank Matrix.
 */
vector<vector<float>> pageRankMatrixDataStructure(const vector<vector<int>>& indexes){

    vector<vector<float>> matrix (indexes.size(), vector<float>(indexes.size(), 0.0));
    for(int i = 0; i < matrix.size(); i++){
        if(!indexes.at(i).empty()){
            for(int j = 0; j < indexes.at(i).size(); j++){
                matrix.at(i).at(indexes.at(i).at(j)) = 1.0 / (float) indexes.at(i).size();
            }
        }
        else{
            for(int j = 0; j < matrix.size(); j++){
                matrix.at(i).at(j) = 1.0 / (float) matrix.size();
            }
        }
    }
    return matrix;
}

/**
 * @brief Initialization of the column wise data structure.
 *
 * This function calls three functions:
 *  - 'readingNormalizingEdges' in order to normalize the ids of the nodes and identify their exact cardinality;
 *  - 'identifyingFullRows' in order to understand which is the destination node of each link;
 *  - 'creatingDataStructure' in order to build the column wise data structure.
 *
 * @param param1 The path of the file containing the edges of the graph.
 * @return A pair composed by an integer and a vector of vectors of integers:
 * - the integer represents the total number of nodes composing the graph;
 * - the vector of vectors of integer represents the column wise data structure.
 */
pair<int, vector <vector<int>>> columnInitialization(const string& filePath){

    cout << "Reading file " << filePath << endl;
    pair < int, vector < pair < int, int>>> edges = readingNormalizingEdges(filePath);

    cout << "Linearizing columns" << endl;
    vector<pair<int, vector<int>>> linearColumns = identifyingFullRows(edges.second);

    cout << "Creating column wise approach" << endl;
    vector<vector<int>> columnWise = columnWiseDataStructure(edges.first, linearColumns);

    return make_pair(edges.first, columnWise);
}

/**
 * @brief Initialization of the adjacency matrix data structure.
 *
 * This function calls two functions:
 * - 'columnInitialization' in order to understand from which to which node the every link goes;
 * - 'pageRankMatrixDataStructure' in order to fill the Page Rank matrix with the correct values.
 *
 * @param param1 The path of the file containing the edges of the graph.
 * @return A vector of vectors of float representing the Page Rank matrix.
 */
vector<vector<float>> matrixInitialization(const string& filePath){

    pair<int, vector<vector<int>>> dataStructure = columnInitialization(filePath); //column
    cout << "Creating the Page Rank matrix" << endl;
    vector<vector<float>> matrix = pageRankMatrixDataStructure(dataStructure.second);

    return matrix;
}

/**
 * @brief Creating the probability distribution vector.
 *
 * Each element of the vector is equal to 1 over the number of nodes.
 *
 * @param param1 An integer representing the number of nodes of the graph.
 * @return A vector whose elements follow a uniform probability distribution.
 */
vector<float> probabilityDistributionVector(int numNodes){
    return vector<float>(numNodes, (float) (1.0/numNodes));
}

//COMPUTATION

/**
 * @brief Computing the Google Metrics for the matrix data structure.
 *
 * Initially, the number of threads are set and then the actual computation starts:
 * - The scalar product between two vectors are made (the product is made by scanning each column of the matrix);
 * - The 'beta' and 'teleport' constants are applied.
 * If this function is run in parallel then even the time of each thread is taken and saved inside the 'timeThreads' vector.
 *
 * @param param1 A vector of vectors of floats representing the Page Rank Matrix.
 * @param param2 A vector of floats which initially is the uniform probability distribution vector.
 * @param param3 A vector of floats which will contain the result of the computation.
 * @param param4 A float which is a constant.
 * @param param5 A float which is a constant.
 * @param param6 An int which specifies the number of threads with which the computation must be done.
 * @param param7 A vector of doubles which will contain the execution times of the threads.
 */
void pageRankMatrix(const vector<vector<float>>& m, const vector<float>& v, vector<float>& res, const float& beta, const float& teleport, const int& numThreads, vector<double>& timeThreads){

    omp_set_num_threads(numThreads);
    #pragma omp parallel
    {
        double startTime = omp_get_wtime();
        #pragma omp for ordered
        for (int col = 0; col < m.size(); col++) {
            for (int row = 0; row < m.size(); row++) {
                res.at(row) += v.at(col) * m.at(col).at(row);
            }
        }
        double endTime = omp_get_wtime();
        double executionTime = endTime - startTime;
        #pragma omp critical
        {
            timeThreads.push_back(executionTime);
        }
    }

    omp_set_num_threads(numThreads);
    #pragma omp parallel for
    for(int i = 0; i < res.size(); i++){
        res.at(i) = (beta * res.at(i)) + teleport;
    }
}

/**
 * @brief Computing the Google Metrics for the column wise data structure.
 *
 * Initially, the number of threads are set and then the actual computation starts:
 * - For each node is checked if it has some outgoing links:
 *      - If it has them, then the classic scalar product is made;
 *      - Otherwise, it is a dead end and its contribution is saved in a support variable.
 * - The contribution of the dead ends is added and the 'beta' and 'teleport' constants are applied.
 * If this function is run in parallel then even the time of each thread is taken and saved inside the 'timeThreads' vector.
 *
 * @param param1 A vector of vectors of ints where:
 * - The i-th element of the external vector corresponds to the i-th node;
 * - The i-th internal vector contains the ids of the nodes reached by the i-th node, if any.
 * @param param2 A vector of floats which initially is the uniform probability distribution vector.
 * @param param3 A vector of floats which will contain the result of the computation.
 * @param param4 A float which is a constant.
 * @param param5 A float which is a constant.
 * @param param6 An int which specifies the number of threads with which the computation must be done.
 * @param param7 A vector of doubles which will contain the execution times of the threads.
 */
void pageRankColumn(const vector<vector<int>>& indexes, const vector<float>& v, vector<float>& res, const float& beta, const float& teleport, const int& numThreads, vector<double>& timeThreads){

    float deadEnds = 0.0;
    omp_set_num_threads(numThreads);
    #pragma omp parallel
    {
        double startTime = omp_get_wtime();
        #pragma omp for ordered
        for (int col = 0; col < indexes.size(); col++) {
            if (!indexes.at(col).empty()) {
                size_t length = indexes.at(col).size();
                for (int row = 0; row < length; row++) {
                    #pragma omp atomic
                    res.at(indexes.at(col).at(row)) += v.at(col) * (float) (1.0 / (float) length);
                }
            } else {
                deadEnds += v.at(col) * (float) (1.0 / (float) res.size());
            }
        }
        double endTime = omp_get_wtime();
        double executionTime = endTime - startTime;
        #pragma omp critical
        {
            timeThreads.push_back(executionTime);
        }
    }
    omp_set_num_threads(numThreads);
    #pragma omp parallel for
    for(int i = 0; i < res.size(); i++){
        res.at(i) += deadEnds;
        res.at(i) = (beta * res.at(i)) + teleport;
    }

}
