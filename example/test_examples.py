# Undirected graph
id_ = "ud_graph"
directed = False
vertex_num = 6
vertex_list = [f"V{i}" for i in range(vertex_num)]
adj_matrix = [
    [0, 10, 0, 6, 0, 0], 
    [10, 0, 5, 0, 3, 1],
    [0, 5, 0, 0, 9, 0], 
    [6, 0, 0, 0, 15, 0], 
    [0, 3, 9, 15, 0, 0], 
    [0, 1, 0, 0, 0, 0]
]

# Directed acyclic graph, used for topological sorting.
id_ = "dag"
directed = True
vertex_num = 15
vertex_list = [f"V{i}" for i in range(vertex_num)]
adj_matrix = [
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] 
]

# Activity on edge (AOE) network, used for critical path search.
# AOE network is a DAG with only one vertex with zero in-degree and only one vertex with zero out-degree.
id_ = "aoe_network"
directed = True
vertex_num = 9
vertex_list = [f"V{i}" for i in range(vertex_num)]
adj_matrix = [
    [0, 6, 4, 5, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 7, 5, 0],
    [0, 0, 0, 0, 0, 0, 0, 4, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 2],
    [0, 0, 0, 0, 0, 0, 0, 0, 4],
    [0, 0, 0, 0, 0, 0, 0, 0, 0]
]


if __name__ == "__main__":
    g = Graph(id_, directed, vertex_num, vertex_list, adj_matrix)
    g.visualize()
