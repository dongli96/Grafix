import graphviz
from typing import *

import graph


class Graph(object):
    
    def __init__(self, *args, **kw):
        if "graph_obj" not in kw:
            self.graph_ = graph.Graph(*args, **kw)
        else:
            self.graph_ = kw["graph_obj"]

    @property
    def id(self):
        return self.graph_.get_id()

    @property
    def is_directed(self):
        return self.graph_.is_directed()

    @property
    def vertex_list(self) -> List[str]:
        return self.graph_.get_vertex_list()

    @property
    def adj_matrix(self) -> List[List[int]]:
        return self.graph_.get_adj_matrix()

    @property
    def edgeset_array(self) -> List[List[int, int, int]]:
        edgeset_array_ = []
        edgeset_array = self.graph_.get_edgeset_array()
        vertex_list = self.vertex_list
        for edge in edgeset_array:
            edgeset_array_.append([vertex_list[edge[0]], vertex_list[edge[1]], edge[2]])
        return edgeset_array_

    @property
    def adj_list(self) -> List[List]:
        adj_list = self.graph_.get_adj_list()
        adj_list_ = []
        for node in adj_list:
            node_ = [self.vertex_list[node.get_idx()], node.get_next()]
            adj_list_.append(node_)
        return adj_list_

    @property
    def in_degree(self) -> List[int]:
        return self.graph_.get_in_degree()

    @property
    def out_degree(self) -> List[int]:
        return self.graph_.get_out_degree()

    @property
    def is_connected(self):
        return self.graph_.is_connected()

    @property
    def is_strongly_connected(self):
        return self.graph_.is_strongly_connected()

    @property
    def is_cyclic(self):
        return self.graph_.is_cyclic()
    
    def add_vertex(self, add_v):
        self.graph_.add_vertex(add_v)

    def add_edge(self, add_e):
        self.graph_.add_edge(add_e)

    def del_vertex(self, del_v):
        self.graph_.delete_vertex(del_v)

    def del_edge(self, del_e):
        self.graph_.delete_edge(del_e)

    def info(self):
        self.graph_.info()

    def random_init(self):
        self.graph_.random_init()

    def traverse(self, algo="dfs") -> List[List[str]]:
        """ Traverses the current graph.

        Args:
            algo (str, optional): Algorithm used. Defaults to "dfs",
            which is depth first search, another is "bfs" (breath 
            first search).

        Returns:
            Traversal sequence of all vertexes. 
        """
        return self.graph_.traverse(algo)

    def shortest_path(self, start_idx, algo="f"):
        """ If algo == "d" (Dijkstra), prints shortest paths from 
            the specified vertexes to other vertexes. 
            If algo == "f (Floyd), prints the shortest paths from 
            all vertexes to all other vertexes. 

        Args:
            start_idx (uint): start vertex for Prim algorithm.
            algo (str, optional): Algorithm used. Defaults to "f".
        """
        self.graph_.shortest_path(start_idx, algo)

    def ms_tree(self, algo="p"):
        """ Generates a `Py::Graph` object encapsulating the minimum
            spanning tree of the current graph.

        Args:
            algo (str, optional): Alogorithm used. Defaults to "p", 
            which is 'Prim', another choice is "k" (Kruskal).

        Returns:
            Graph: Minimum spanning tree of the current graph.
        """
        return Graph(graph_obj=self.graph_.minimum_spanning_tree(algo))

    def topo_sort(self) -> List[int]:
        return self.graph_.topological_sort()

    def critical_path(self):
        return Graph(graph_obj=self.graph_.critical_path())

    def visualize(self, output_name=None):
        graph_ = self.graph_
        adj_list = graph_.get_vertex_list()
        vertex_num = len(adj_list)
        adj_matrix = graph_.get_adj_matrix()
        if self.is_directed:
            dot = graphviz.Digraph()
            for i in range(vertex_num):
                for j in range(vertex_num):
                    degree = adj_matrix[i][j]
                    if degree:
                        dot.edge(adj_list[i], adj_list[j], str(degree))
        else:
            dot = graphviz.Graph()
            for i in range(vertex_num):
                for j in range(i + 1, vertex_num):
                    degree = adj_matrix[i][j]
                    if degree:
                        dot.edge(adj_list[j], adj_list[i], str(degree))
        for node_id in adj_list:
            dot.node(node_id)
        if output_name is None:
            output_name = graph_.get_id()
        dot.render(output_name, format="png", cleanup=True)
