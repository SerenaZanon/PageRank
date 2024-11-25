# University course: LEARNING WITH MASSIVE DATA
# Assignment: Multi-threaded implementation of PageRank

PageRank measures the importance of a node i in a graph G as the weighted sum of the importance of its neighbours.

Let v be a vector storing the importance of each node:
- v[i] = importance of the i-th node, and elements of v sum up to 1
- The contribution of each neighbor j is normalized by its out-degree o(j)
  $v^{t+1}[i] = \sum_{j \arrow i} \frac{v^t[j]}{o(j)}$

Note that the above update rule can be rewritten as a matrix-vector multiplication:

$v^{t+1} = Mv^t$ with 
