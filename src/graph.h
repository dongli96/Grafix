#include <string>
#include <vector>
using namespace std;

#ifndef GRAPH_H
#define GRAPH_H

using lst_uint = vector<unsigned>;
using lst_str = vector<string>;
using matrix_uint = vector<lst_uint>;
using matrix_str = vector<lst_str>;

struct AdjacencyListNode {
    unsigned idx;
    AdjacencyListNode* next;
    unsigned get_idx() { return idx; }
    AdjacencyListNode* get_next() { return next; }
};
using adj_lst = vector<AdjacencyListNode*>;

class Graph {
private:
    string id = "graph1";
    bool directed = false;
    unsigned vertex_num = 0;
    lst_str vertex_list;
    matrix_uint adjacency_matrix;
public:
    Graph() = default;
    Graph(const string& id_, const bool& directed_, const unsigned& vertex_num_): id(id_), directed(directed_), vertex_num(vertex_num_) {};
    Graph(const string& id_, const bool& directed_, const unsigned& vertex_num_, const lst_str& vertex_list_, const matrix_uint& adjacency_matrix_); 
    virtual ~Graph() {};
    void info();
    string get_id();
    lst_str get_vertex_list();
    matrix_uint get_adj_matrix();
    void add_vertex(lst_str add_v);
    void add_edge(matrix_uint add_e);
    void delete_vertex(lst_uint del_v);
    void delete_edge(matrix_uint del_e);
    matrix_uint get_edgeset_array();
    adj_lst get_adj_list();
    lst_uint get_in_degree();
    lst_uint get_out_degree();
    bool is_directed();
    bool is_cyclic();
    bool is_connected();
    bool is_strongly_connected();
    void random_init();
    matrix_uint traverse(const string& algo="dfs");
    void shortest_path(const unsigned& start_idx, const string& algo="d");
    Graph minimum_spanning_tree(const string& algo="p");
    lst_uint topological_sort();
    Graph critical_path();
};

#endif
