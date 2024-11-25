The project is characterized by diferrent files:

C++ FILES
- utility.h and utilityImplementation.cpp contain respectively the declaration and the implementation of all the functions necessary to compute the Page Rank.
- pageRankMatrix.cpp and pageRankColumnWise.cpp contain repsectively the main for the Page Rank Matrix Data Structure and the Column Wise Data Structure. In both there is:
	- The while that loops 50 times over the associate pageRank function;
	- A for loop that repeats the computation 5 times, in order to have a better representation of the execution time (a mean it's taken);
	- A for loop over all the possible number of threads;
	- A for loop over all the possible datasets on which the function can be tested.
- dataStructureSize.cpp allows to measure the size in byte of the Data Structure for all the datasets.

PYTHON FILE
- DenseGraphsCreation.py contains the procedure with which the tested datasets were made.

JUPYTER NOTEBOOK FILE
- analysis.ipynb contains the instructions to read the csv files created by the cpp files in order to visualize the Speed Up and the execution time of the threads.

TXT FILE
- constructionTime.txt contains the construction times retrieved by the execution of the cpp files and the space occupied by the Data Structures.

datasets FOLDER
This folder contains:
- 'adHoc' folder which contains the datasets that were tested on both algorithms;
- 'generic' folder which contains the datasets that were tested only on the Column Wise algorithm;
- csv files obtained by running the cpp files;
- chosenDatasets.txt contains information about the datasets tested.


The project can be runned with the following commands:
1) g++ -O3 -fopenmp -std=c++11 pageRankMatrix.cpp utilityImplementation.cpp -o matrix
2) ./matrix

1) g++ -O3 -fopenmp -std=c++11 pageRankColumnWise.cpp utilityImplementation.cpp -o column
2) ./column

1) g++ -O3 -fopenmp -std=c++11 dataStructureSize.cpp utilityImplementation.cpp -o dataStructure
2) ./dataStructure

!!!
How was the project run
!!!
In order to make things going faster and to try an higher number of threads, the results of the project were obtained by running it on the DAIS cluster with the following commands:
1) g++ -O3 -fopenmp -std=c++11 pageRankMatrix.cpp utilityImplementation.cpp -o matrix
2) srun -v -l -c 20 ./matrix 20

1) g++ -O3 -fopenmp -std=c++11 pageRankColumnWise.cpp utilityImplementation.cpp -o column
2) srun -v -l -c 20 ./column 20