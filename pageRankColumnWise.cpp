#include <chrono>
#include <algorithm>
#include "utility.h"


int main(){

    /**
     * Retrieving the datasets' paths.
     * */
    vector<string> datasetNames = retrieveAllFolders();
    //vector<string> datasetNames;
    //datasetNames.push_back("datasets/generic/smallDataset_sparse.txt");
    /**
     * 'allTimePoints' is a vector of vectors of float where:
     * - Each internal vector correspond to a different dataset;
     * - Each internal vector will contain the mean of the time required for each execution according to the number of threads used
     * The size of the internal vector depends on the number of threads that can be run (in this case, 20).
     * */
    vector<vector<float>> allTimePoints;
    /**
     * 'timeThreads' is a vector of vectors of pairs of double where:
     * - Each internal vector correspond to a different dataset;
     * - Each internal vector will contain the minimum and the maximum time required by a thread during the computation.
     * The size of the internal vector depends on the number of threads that can be run (in this case, 20).
     * */
    vector<vector<pair<double, double>>> timeThreads;

    int maxNumThread = 1;
    int maxNumRound = 1;

    for(const string& dataset: datasetNames) {

        /**
         * Initialization of the Column Wise data structure.
         * */
        auto startInit = std::chrono::high_resolution_clock::now();
        pair<int, vector<vector<int>>> dataStructure = columnInitialization(dataset);
        auto endInit = std::chrono::high_resolution_clock::now();
        auto elapsedInit = std::chrono::duration_cast<std::chrono::milliseconds>(endInit - startInit);
        cout<<"Construction Time: "<<elapsedInit.count()<<" ms"<<endl;

        /**
         * 'timePoints' is a vector of float that will contain the mean time required to compute the instructions.
         * The size of this vector depends on the number of threads that can be run (in this case, 20).
         * */
        vector<float> timePoints;
        /**
         * 'bigGapThreads' is a vector of pair of double.
         * Each pair will contain as first element the minimum time required by a thread and
         * as second element the maximum time required by a thread.
         * The size of this vector depends on the number of threads that can be run (in this case, 20).
         * */
        vector<pair<double, double>> bigGapThreads;

        for (int thread = 1; thread <= maxNumThread; thread++) {

            cout<<"N. threads: "<<thread<<endl;
            float meanTime = 0.0;
            /**
             * 'timeThreadsRounds' is a vector of pair of double.
             * Each pair will contain as first element the minimum time required by a thread and
             * as second element the maximum time required by a thread.
             * The size of this vector depends on the number of rounds that can be run (in this case, 5).
             * */
            vector<pair<double, double>> timeThreadsRounds;

            for (int round = 0; round < maxNumRound; round++) {

                cout<<"N. round:"<<round<<endl;
                /**
                 * 'threads' is a vector of double.
                 * It will contain the time required by each thread.
                 * */
                vector<double> threads;
                vector<float> v = probabilityDistributionVector(dataStructure.first);
                float beta = 0.85;
                float teleport_constant = (float) (1.0 - beta) / (float) dataStructure.first;
                vector<float> res(dataStructure.first, 0.0);

                int iter = 0;
                /**
                 * Auxiliary variables.
                 * */
                double minTime;
                double maxTime;

                /**
                 * Iterating for 50 times the 'pageRankMatrix' function.
                 * */
                auto startComp = std::chrono::high_resolution_clock::now();
                while(iter < 50) {
                    res.assign(res.size(), 0.0);
                    pageRankColumn(dataStructure.second, v, res, beta, teleport_constant, thread, threads);
                    v = res;
                    iter++;
                }
                auto endComp = std::chrono::high_resolution_clock::now();
                auto elapsedComp = std::chrono::duration_cast<std::chrono::milliseconds>(endComp - startComp);
                meanTime += elapsedComp.count();

                for(int i = 0; i < 10; i++){
                    cout<<v.at(i)<<endl;
                }

                /**
                 * Identifying the minimum and the maximum time required by a thread.
                 * */
                minTime = *min_element(threads.begin(), threads.end());
                maxTime = *max_element(threads.begin(), threads.end());
                timeThreadsRounds.emplace_back(minTime, maxTime);
            }
            /**
             * Appending the mean of the total time for the number of rounds.
             * */
            timePoints.push_back(meanTime/(float) maxNumRound);
            /**
             * Appending the pair whose gap is the biggest in the 'timeThreadRounds' vector.
             * */
            bigGapThreads.push_back(greatestDifference(timeThreadsRounds));
        }
        /**
         * Appending the vectors for the specific dataset.
         * */
        allTimePoints.push_back(timePoints);
        timeThreads.push_back(bigGapThreads);
    }
    /**
     * Saving data.
     * */
    savingIntoCSV("datasets/executionTimesColumn.csv", allTimePoints);
    savingPairsIntoCSV("datasets/threadTimesColumn.csv", timeThreads);

    return 0;
}