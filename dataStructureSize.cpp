#include <algorithm>
#include<chrono>
#include <cmath>
#include "utility.h"


using namespace std;

int main(){

    /**
     * Retrieving the datasets' paths.
     * */
    vector<string> datasetNamesColumns = retrieveAllFolders();
    vector<string> datasetNameMatrix = retrieveAdHocFolder();

    /**
     * In 'spaceMatrix' there will be the space occupied by the datasets based on the Page Rank Matrix Data Structure.
     * In 'spaceColumn' there will be the space occupied by the datasets based on the Column Wise Data Structure.
     * */
    vector<size_t> spaceMatrix;
    vector<size_t> spaceColumn;

    /**
     * For each dataset the total space is computed.
     * It is given by the sum of the space of
     * - The external vector
     * - The sum of the space occupied by the internal vectors
     * - The sum of the space occupied by the floats.
     * */
    for(string dataset: datasetNameMatrix) {
        vector<vector<float>> matrix = matrixInitialization(dataset);
        size_t totalSpace = sizeof(matrix);
        totalSpace += matrix.size() * sizeof(vector<float>);
        for (vector<float> v: matrix) {
            totalSpace += sizeof(float) * v.size();
        }
        spaceMatrix.push_back(totalSpace);
    }

    /**
     * For each dataset the total space is computed.
     * It is given by the sum of the space of
     * - The external vector
     * - The sum of the space occupied by the internal vectors
     * - The sum of the space occupied by the ints.
     * */
    for(string dataset: datasetNamesColumns) {
        pair<int, vector<vector<int>>> dataStructure = columnInitialization(dataset);
        size_t totalSpace = sizeof(dataStructure.second);
        totalSpace += dataStructure.second.size() * sizeof(vector<int>);
        for (vector<int> v: dataStructure.second) {
            totalSpace += sizeof(int) * v.size();
        }
        spaceColumn.push_back(totalSpace);
    }

    /**
     * Visualizing results
     * */
    for(int i = 0; i < datasetNameMatrix.size(); i++){
        cout<<"dataset: "<<datasetNameMatrix.at(i)<<", spaceMatrix: "<<spaceMatrix.at(i)<<endl;
    }

    for(int i = 0; i < datasetNamesColumns.size(); i++){
        cout<<"dataset: "<<datasetNamesColumns.at(i)<<", spaceColumn: "<<spaceColumn.at(i)<<endl;
    }

    return 0;
}