# University course: LEARNING WITH MASSIVE DATA
# Assignment: Multi-threaded implementation of PageRank

PageRank measures the importance of a node i in a graph G as the weighted sum of the importance of its neighbours.

Let v be a vector storing the importance of each node:
- v[i] = importance of the i-th node, and elements of v sum up to 1
- The contribution of each neighbor j is normalized by its out-degree o(j)

  $v^{t+1}[i] = \sum_{j \rightarrow i} \frac{v^t[j]}{o(j)}$

Note that the above update rule can be rewritten as a matrix-vector multiplication:

$$
v^{t+1} = M v^t
$$

with 

$$
    \[
        M[i, j] =
        \begin{cases}
            \text{$\frac{1}{o(j)}$ if $j \xrightarrow{} i$} \\
            \text{$\frac{1}{n}$ if o(j) = 0} \\
            \text{0 otherwise}
        \end{cases}
    \]
$$



The above update rule, defines an iterative process:
- we start from a random vector v0, ( v0 sums up to 1, it is a probability distribution)
- the update rule is applied for several iterations (<=50) until convergence
- a.k.a. random surfer model

Does it converge?

Only if the original graph G is irreducible (all states are reachable from any other state) and aperiodic (no “cycles” of fixed length)

The Web Graph does not satisfy this criteria (it has dead ends and cycles)!

Solution (leading to the so called Google Matrix):
- Teleportation and dead ends removal

$v^{t+1} = \beta Mv^t + (1 - \beta)\frac{1}{n}$ with 
$$
M[i, j] =
\begin{cases} 
\frac{1}{o(j)} & \text{se } j \to i \\ 
\frac{1}{n} & \text{se } o(j) = 0 \\ 
0 & \text{altrimenti}
\end{cases}
$$

- with probability β we follow M, with probability (1-β) we jump to a random node
- dead end nodes link to all nodes of the graph

For space issues, the txt file are omitted in the repository.
