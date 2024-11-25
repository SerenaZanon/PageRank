import networkx as nx

def generate_directed_graph_with_density(num_nodes, density):
    num_edges = int(density * num_nodes * (num_nodes - 1))  # Calculate the number of edges based on density
    print(num_edges)
    G = nx.gnm_random_graph(num_nodes, num_edges, directed=True)
    return G

def write_edges_to_file(graph, filename):
    with open(filename, 'w') as file:
        for edge in graph.edges():
            file.write(f"{edge[0]} {edge[1]}\n")

num_nodes = 7000
desired_density = 0.25  # Specify the desired density of the graph

graph = generate_directed_graph_with_density(num_nodes, desired_density)
print("Generated graph with density", nx.density(graph))

output_file = "smallDataset_025.txt"
write_edges_to_file(graph, output_file)
print(f"Edges written to {output_file}")
