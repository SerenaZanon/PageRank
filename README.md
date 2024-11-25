# University course: LEARNING WITH MASSIVE DATA
# Assignment: Multi-threaded implementation of PageRank

PageRank measures the importance of a node i in a graph G as the weighted sum of the importance of its neighbours.

Let v be a vector storing the importance of each node:
- v[i] = importance of the i-th node, and elements of v sum up to 1
- The contribution of each neighbor j is normalized by its out-degree o(j)
  $v^{t+1}[i] = \sum_{j \rightarrow i} \frac{v^t[j]}{o(j)}$

Note that the above update rule can be rewritten as a matrix-vector multiplication:

$v^{t+1} = Mv^t$ with 
$
M[i, j] =
\left\{
\begin{array}{ll}
\frac{1}{o(j)} & \text{if } j \rightarrow i \\
0 & \text{otherwise}
\end{array}
\right.
$

The above update rule, defines an iterative process:
- we start from a random vector v0, ( v0 sums up to 1, it is a probability distribution)
- the update rule is applied for several iterations (<=50) until convergence
- a.k.a. random surfer model

Does it converge?

Only if the original graph G is irreducible (all states are reachable from any other state) and aperiodic (no “cycles” of fixed length)

The Web Graph does not satisfy this criteria (it has dead ends and cycles)!

Solution (leading to the so called Google Matrix):
- Teleportation and dead ends removal

$v^{t+1} = \beta Mv^t + (\ - \beta)\frac{1}{n}$ with 
$
M[i, j] =
\left\{
\begin{array}{ll}
\frac{1}{o(j)} & \text{if } j \rightarrow i \\
\frac{1}{n} $ \text{if } o(j) = 0 \\
0 & \text{otherwise}
\end{array}
\right.
$

- with probability β we follow M, with probability (1-β) we jump to a random node
- dead end nodes link to all nodes of the graph
