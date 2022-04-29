#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <stack>
#include <stdexcept>
#include <unordered_set>

#include "graph.h"

#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"

namespace py = pybind11;

void print(const lst_uint& obj) {
    for (auto& elem: obj) {
        cout << elem << " ";
    }
    cout << endl;
}

void print(const lst_str& obj) {
    for (auto& elem: obj) {
        cout << elem << " ";
    }
    cout << endl;
}

void print(const matrix_uint& obj) {
    for (auto& line: obj) {
        for (auto& elem: line) {
            cout << elem << " ";
        }
        cout << endl;
    }
}

bool compare_by_last_element(lst_uint a, lst_uint b) {
    auto elem_a = *(a.end() - 1);
    auto elem_b = *(b.end() - 1);
    if (elem_a < elem_b) {
        return true;
    }
    else {
        return false;
    }
}

Graph::Graph(const string& id_, const bool& directed_, const unsigned& vertex_num_, const lst_str& vertex_list_, const matrix_uint& adjacency_matrix_) {
    if (vertex_list_.size() != vertex_num_) {
        throw invalid_argument("The length of \"vertex_list_\" not matched with \"vertex_num_\".");
    }
    if (adjacency_matrix_.size() != vertex_num_) {
        throw invalid_argument("The length of \"adjacency_matrix_\" not matched with \"vertex_num_\".");
    }
    else {
        for (auto& line: adjacency_matrix_) {
            if (line.size() != vertex_num_) {
                throw invalid_argument("The shape of \"adjacency_matrix_\" not suitable, should be: \"vertex_num_\" * \"vertex_num_\".");
            }
        }
    }
    if (!directed_) {
        for (unsigned i = 0; i < vertex_num_; ++i) {
            if (adjacency_matrix_[i][i]) {
                throw invalid_argument("Self weight should always be zero.");
            }
        }
        for (unsigned i = 0; i < vertex_num_; ++i) {
            for (unsigned j = i + 1; j < vertex_num_; ++j) {
                if (adjacency_matrix_[i][j] != adjacency_matrix_[j][i]) {
                    throw invalid_argument("\"adjacency_matrix_\" implies a directed graph.");
                }
            }
        }
    }
    id = id_;
    directed = directed_;
    vertex_num = vertex_num_;
    vertex_list = vertex_list_;
    adjacency_matrix = adjacency_matrix_;
}

void Graph::info() {
    cout << "ID: " << id << endl;
    cout << "Directed: " << directed << endl;
    cout << "Number of Vertex: " << vertex_num << endl;
    cout << "Adjacency list: ";
    print(vertex_list);
    cout << "Adjacency matrix:" << endl;
    print(adjacency_matrix);
}

string Graph::get_id() {
    return id;
}

lst_str Graph::get_vertex_list() {
    return vertex_list;
}

matrix_uint Graph::get_adj_matrix() {
    return adjacency_matrix;
}

// O(V^2)
void Graph::add_vertex(lst_str add_v) {
    vertex_num += add_v.size();
    for (auto& v: add_v) {
        vertex_list.push_back(v);
    }
    matrix_uint new_adj_mat(vertex_num);
    for (unsigned i = 0; i < vertex_num; ++i) {
        new_adj_mat[i].resize(vertex_num);
        if (i < adjacency_matrix.size()) {
            for (unsigned j = 0; j < vertex_num; ++j) {
                if (j < adjacency_matrix.size()) {
                    new_adj_mat[i][j] = adjacency_matrix[i][j];
                }
                else {
                    new_adj_mat[i][j] = 0;
                }
            } 
        }
        else {
            for (auto& w: new_adj_mat[i]) {
                w = 0;
            }
        }
    }
    adjacency_matrix = new_adj_mat;
}

// O(A)
void Graph::add_edge(matrix_uint add_e) {
    for (auto& e: add_e) {
        auto begin = e[0];
        auto end = e[1];
        if (begin >= vertex_list.size() || end >= vertex_list.size()) {
            throw invalid_argument("Cannot add edge between nonexistent vertexes.");
        }
        if (directed) {
            adjacency_matrix[begin][end] = e[2];
        }
        else {
            adjacency_matrix[begin][end] = adjacency_matrix[end][begin] = e[2];
        }
    }
}

// O(V^2)
void Graph::delete_vertex(lst_uint del_v) {
    vertex_num -= del_v.size();
    lst_str new_vertex_list;
    for (unsigned i = 0; i < vertex_list.size(); ++i) {
        if (find(del_v.begin(), del_v.end(), i) == del_v.end()) {
            new_vertex_list.push_back(vertex_list[i]);
        }
    }
    matrix_uint new_adj_mat;
    for (unsigned i = 0; i < adjacency_matrix.size(); ++i) {
        if (find(del_v.begin(), del_v.end(), i) == del_v.end()) {
            lst_uint line;
            for (unsigned j = 0; j < adjacency_matrix.size(); ++j) {
                if(find(del_v.begin(), del_v.end(), j) == del_v.end()) {
                    line.push_back(adjacency_matrix[i][j]);
                }
            }
            new_adj_mat.push_back(line);
        }
    }
    vertex_list = new_vertex_list;
    adjacency_matrix = new_adj_mat;
}

// O(A)
void Graph::delete_edge(matrix_uint del_e) {
    for (auto& e: del_e) {
        auto begin = e[0];
        auto end = e[1];
        if (begin >= vertex_list.size() || end >= vertex_list.size()) {
            throw invalid_argument("Cannot delete edge between nonexistent vertexes.");
        }
        if (directed) {
            adjacency_matrix[begin][end] = 0;
        }
        else {
            adjacency_matrix[begin][end] = adjacency_matrix[end][begin] = 0;
        }
    }
}

// O(V^2)
matrix_uint Graph::get_edgeset_array() {
    matrix_uint edgeset_array; 
    if (directed) {
        for (unsigned i = 0; i < vertex_num; ++i) {
            for (unsigned j = 0; j < vertex_num; ++j) {
                if (adjacency_matrix[i][j]) {
                    lst_uint edge = {i, j, adjacency_matrix[i][j]};
                    edgeset_array.push_back(edge);
                }
            }
        }
    } 
    else {
        for (unsigned i = 0; i < vertex_num; ++i) {
            for (unsigned j = i + 1; j < vertex_num; ++j) {
                if (adjacency_matrix[i][j]) {
                    lst_uint edge = {i, j, adjacency_matrix[i][j]};
                    edgeset_array.push_back(edge);
                }
            }
        }
    }
    return edgeset_array;
}

adj_lst Graph::get_adj_list() {
    adj_lst adj_list;
    for (unsigned i = 0; i < vertex_num; ++i) {
        AdjacencyListNode* p_head = new AdjacencyListNode;
        p_head->idx = i;
        AdjacencyListNode* last = p_head;   
        for (unsigned j = 0; j < vertex_num; ++j) {
            if (adjacency_matrix[i][j]) {
                AdjacencyListNode* p_node = new AdjacencyListNode;
                p_node->idx = j;
                last->next = p_node;
                last = p_node;
            }
        }
        last->next = 0;
        adj_list.push_back(p_head);
    }
    return adj_list;
}

lst_uint Graph::get_in_degree() {
    lst_uint in_degrees(vertex_num, 0);
    for (unsigned j = 0; j < vertex_num; ++j) {
        unsigned count = 0;
        for (unsigned i = 0; i < vertex_num; ++i) {
            if (adjacency_matrix[i][j] > 0) {
                ++count;
            }
        }
        in_degrees[j] = count;
    }
    return in_degrees; 
}

lst_uint Graph::get_out_degree() {
    lst_uint out_degrees(vertex_num, 0);
    for (unsigned i = 0; i < vertex_num; ++i) {
        unsigned count = 0;
        for (unsigned j = 0; j < vertex_num; ++j) {
            if (adjacency_matrix[i][j] > 0) {
                ++count;
            }
        }
        out_degrees[i] = count;
    }
    return out_degrees;
}

bool Graph::is_directed() {
    return directed;
}

bool Graph::is_connected() {
    if (!directed) {
        unordered_set<unsigned> visited;
        lst_uint vertexes;
        visited.insert(0);
        queue<unsigned> Q;
        Q.push(0);
        while (!Q.empty()) {
            auto front = Q.front();
            vertexes.push_back(front);
            Q.pop();
            for (unsigned j = 0; j < vertex_num; ++j) {
                auto degree = adjacency_matrix[front][j];
                if (degree && visited.find(j) == visited.end()) {
                    Q.push(j);
                    visited.insert(j);
                }
            }
        }
        if (vertexes.size() != vertex_num) {
            return false;
        }
        else {
            return true;
        }
    }
    else {
        matrix_uint temp = adjacency_matrix;
        for (unsigned i = 0; i < vertex_num; ++i) {
            for (unsigned j = i + 1; j < vertex_num; ++j) {
                if (temp[i][j]) {
                    temp[i][j] = temp[j][i] = 1;
                }
            }
        }
        auto g = Graph();
        g.directed = false;
        g.vertex_list = vertex_list;
        g.adjacency_matrix = temp;
        g.is_connected();
    }
}

bool Graph::is_strongly_connected() {
    if (!directed) {
        throw invalid_argument("Strong connectivity only applies to directed graphs.");
    }
    else {
        if (traverse("dfs").size() > 1) {
            return false;
        }
        else {
            return true;
        }
    }
}

bool Graph::is_cyclic() {
    unordered_set<unsigned> visited;
    if (!directed) {
        auto connect_components = traverse("dfs");
        for (auto& connect_component: connect_components) {
            Graph g1(*this);
            lst_uint diff;
            for (unsigned i = 0; i < vertex_num; ++i) {
                if (find(connect_component.begin(), connect_component.end(), i) == connect_component.end()) {
                    diff.push_back(i);
                }
            }
            g1.delete_vertex(diff);
            if (g1.get_edgeset_array().size() > g1.vertex_num - 1) {
                return true;
            }
        }
        return false;
    }
    else {
        while (visited.size() != vertex_num) {
            unordered_set<unsigned> loop_visited;
            unsigned start_idx;
            for (unsigned i = 0; i < vertex_num; ++i) {
                if (visited.find(i) == visited.end()) {
                    start_idx = i;
                    break;
                }
            }
            stack<unsigned> S;
            S.push(start_idx);
            visited.insert(start_idx);
            loop_visited.insert(start_idx);
            while (!S.empty()) {
                auto top = S.top();
                bool push = true;
                unsigned j = 0;
                for (; j < vertex_num; ++j) {
                    auto degree = adjacency_matrix[top][j];
                    if (degree) {
                        if (visited.find(j) == visited.end()) {
                            S.push(j);
                            visited.insert(j);
                            push = false;
                            break;
                        }
                        else if (loop_visited.find(j) != loop_visited.end()) {
                            return true;
                        }
                    }
                }
                if (push) {
                    S.pop();
                }
            }
        }
        return false;
    }
}

void Graph::random_init() {
    for (unsigned num = 1; num < vertex_num + 1; ++num) {
        vertex_list.push_back(to_string(num)); //!
    }
    adjacency_matrix.resize(vertex_num);
    for (auto& line: adjacency_matrix) {
        line.resize(vertex_num);
    }
    if (directed) {
        for (unsigned i = 0; i < vertex_num; ++i) {
            for (unsigned j = 0; j < vertex_num; ++j) {
                adjacency_matrix[i][j] = static_cast<unsigned>(rand()); //!
            }
        }
        for (unsigned i = 0; i < vertex_num; ++i) {
            adjacency_matrix[i][i] = 0;
            
        }
    }
    else {
        for (unsigned i = 0; i < vertex_num; ++i) {
            adjacency_matrix[i][i] = 0;
        }
        for (unsigned i = 0; i < vertex_num; ++i) {
            for (unsigned j = 0; j < i; ++j) {
                adjacency_matrix[i][j] = adjacency_matrix[j][i] = static_cast<unsigned>(rand());
            }
        }
    }
}

matrix_uint Graph::traverse(const string& algo) {
    if (algo != "bfs" && algo != "dfs") {
        throw invalid_argument("Expected \"dfs\" or \"bfs\".");
    }
    matrix_uint result; // stores the final result
    unordered_set<unsigned> visited; // stores all visited vertex indices
    if (algo == "bfs") {
        while (visited.size() != vertex_num) {
            // a new loop means traversal of a new connected component
            unsigned start_idx;
            lst_uint vertexes; // stores all vertex indices in a connected component in the traversal order
            for (unsigned i = 0; i < vertex_num; ++i) {
                if (visited.find(i) == visited.end()) {
                    start_idx = i;
                    break;
                }
            }
            queue<unsigned> Q;
            Q.push(start_idx);
            visited.insert(start_idx);
            while (!Q.empty()) {
                auto front = Q.front();
                vertexes.push_back(front);
                Q.pop();
                for (unsigned j = 0; j < vertex_num; ++j) {
                    auto degree = adjacency_matrix[front][j];
                    if (degree && visited.find(j) == visited.end()) {
                        Q.push(j);
                        visited.insert(j);
                    }
                }
            }
            result.resize(result.size() + 1);
            for (auto& idx: vertexes) {
                result[result.size() - 1].push_back(idx);
            }
        }
        return result;
    }
    else {
        while (visited.size() != vertex_num) {
            unsigned start_idx;
            lst_uint vertexes; //stores all vertex indices in a connected component in the traversal order
            for (unsigned i = 0; i < vertex_num; ++i) {
                if (visited.find(i) == visited.end()) {
                    start_idx = i;
                    break;
                }
            }
            stack<unsigned> S;
            S.push(start_idx);
            visited.insert(start_idx);
            vertexes.push_back(start_idx);
            if (!directed) {
                while (!S.empty()) {
                    auto top = S.top();
                    bool push = true;
                    for (unsigned j = 0; j < vertex_num; ++j) {
                        auto degree = adjacency_matrix[top][j];
                        if (degree && visited.find(j) == visited.end()) {
                            S.push(j);
                            visited.insert(j);
                            vertexes.push_back(j);
                            push = false;
                            break;
                        }
                    }
                    if (push) {
                        S.pop();
                    }
                }
                result.resize(result.size() + 1);
                for (auto& idx: vertexes) {
                    result[result.size() - 1].push_back(idx);
                }
            }
            else {
                while (!S.empty()) {
                    auto top = S.top();
                    bool push = true;
                    for (unsigned j = 0; j < vertex_num; ++j) {
                        auto degree = adjacency_matrix[top][j];
                        if (degree && visited.find(j) == visited.end()) {
                            S.push(j);
                            visited.insert(j);
                            vertexes.push_back(j);
                            push = false;
                            break;
                        }
                    }
                    if (push) {
                        S.pop();
                        if (vertexes.size()) {
                            result.resize(result.size() + 1);
                            for (auto& idx: vertexes) {
                                result[result.size() - 1].push_back(idx);
                            }
                            vertexes = lst_uint();
                        }
                    }
                }
            }
        }
        return result;
    }
}

void Graph::shortest_path(const unsigned& start_idx, const string& algo) {
    if (algo != "d" && algo != "f") {
        throw invalid_argument("Only \"d\" or \"f\" is expected.");
    }
    matrix_uint adj_mat_ = adjacency_matrix;
    for (unsigned i = 0; i < vertex_num; ++i) {
        for (unsigned j = 0; j < vertex_num; ++j) {
            if (!adjacency_matrix[i][j]) {
                adj_mat_[i][j] = UINT_MAX;
            }
        }
    }
    if (algo == "d") {
        lst_uint lowcost(adjacency_matrix[start_idx]);
        for (unsigned i = 0; i < vertex_num; ++i) {
            if (lowcost[i] < UINT_MAX) {
                for (unsigned j = 0; j < vertex_num; ++j) {
                    if (adjacency_matrix[i][j] < UINT_MAX) {
                        if (lowcost[i] + adjacency_matrix[i][j] < lowcost[j]) {
                            lowcost[j] = lowcost[i] + adjacency_matrix[i][j];
                        }
                    }
                }
            }
        }
        print(lowcost);
    }
    else if (algo == "f") {
        matrix_uint lowcost(adjacency_matrix);
        for (unsigned i = 0; i < vertex_num; ++i) {
            for (unsigned j = 0; j < vertex_num; ++j) {
                for (unsigned k = 0; k < vertex_num; ++k) {
                    if (lowcost[k][j] < UINT_MAX && lowcost[i][k] < UINT_MAX) {
                        if (lowcost[k][j] + lowcost[i][k] < lowcost[i][j])
                            lowcost[i][j] = lowcost[k][j] + lowcost[i][k];
                    }
                }
            }
        }
        print(lowcost);
    }
}

Graph Graph::minimum_spanning_tree(const string& algo) {
    if (algo != "p" && algo != "k") {
        throw invalid_argument("Only \"prim\" or \"kruskal\" is expected.");
    }
    if (directed || !is_connected()) {
        throw invalid_argument("Only connected undirected graphs have spanning trees.");
    }
    if (algo == "p") {
        lst_uint v_new;
        matrix_uint e_new(vertex_num);
        for (auto& line: e_new) {
            line.resize(vertex_num);
        }
        v_new.push_back(0);
        while (v_new.size() < vertex_list.size()) {
            unsigned min_weight = UINT_MAX;
            lst_uint obj_vertex = {0, 0};
            for (auto& i: v_new) {
                for (unsigned j = 0; j < vertex_num; ++j) {
                    if (find(v_new.begin(), v_new.end(), j) == v_new.end()) {
                        if (adjacency_matrix[i][j] && adjacency_matrix[i][j] < min_weight) {
                            min_weight = adjacency_matrix[i][j];
                            obj_vertex = {i, j};
                        }
                    }
                }
            }
            if (obj_vertex[1]) {
                v_new.push_back(obj_vertex[1]);
                e_new[obj_vertex[0]][obj_vertex[1]] = e_new[obj_vertex[1]][obj_vertex[0]] = min_weight;
            }
        }
        return Graph("mst_prim_1", false, vertex_num, vertex_list, e_new);
    }
    else {
        auto g = Graph("mst_kruskal_1", false, vertex_num);
        g.vertex_list = vertex_list;
        g.adjacency_matrix.resize(vertex_num);
        for (auto& line: g.adjacency_matrix) {
            line.resize(vertex_num);
        } 
        //Kruskal algo is based on edgeset array
        auto edgeset_array = get_edgeset_array();
        sort(edgeset_array.begin(), edgeset_array.end(), compare_by_last_element);
        for (auto& line: edgeset_array) {
            auto begin = line[0];
            auto end = line[1];
            g.adjacency_matrix[begin][end] = line[2];
            if (g.is_cyclic()) {
                g.adjacency_matrix[begin][end] = 0;
            }
        }
        return g;
    }
}

lst_uint Graph::topological_sort() {
    if (!directed) {
        throw invalid_argument("Topological sorting only applies to directed graph.");
    }
    if (is_cyclic()) {
        throw invalid_argument("Topological sorting only applies to acyclic graph.");
    }
    lst_uint result;
    lst_uint in_degrees = get_in_degree();
    auto adj_list = get_adj_list();
    stack<unsigned> S;
    for (unsigned i = 0; i < vertex_num; ++i) {
        if (!in_degrees[i]) {
            S.push(i);
        }
    }
    while (!S.empty()) {
        auto top = S.top();
        S.pop();
        result.push_back(top);
        auto p = adj_list[top]->next;
        while (p) {
            --in_degrees[p->idx];
            if (!in_degrees[p->idx]) {
                S.push(p->idx);
            }
            p = p->next;
        }
    }
    return result;
}

Graph Graph::critical_path() {
    auto in_degrees = get_in_degree();
    if (count(in_degrees.begin(), in_degrees.end(), 0) != 1) {
        throw invalid_argument("Not a AOE network.");
    }
    auto out_degrees = get_out_degree();
    if (count(out_degrees.begin(), out_degrees.end(), 0) != 1) {
        throw invalid_argument("Not a AOE network.");
    }
    lst_uint topo_seq = topological_sort();
    lst_uint etv(vertex_num, 0); //earliest_times_of_vertex
    for (auto& idx: topo_seq) {
        for (unsigned j = 1; j < vertex_num; ++j) {
            if (adjacency_matrix[idx][j]) {
                if (etv[idx] + adjacency_matrix[idx][j] > etv[j]) {
                    etv[j] = etv[idx] + adjacency_matrix[idx][j];
                }
            }
        }
    }
    lst_uint ltv(vertex_num, *(etv.end() - 1)); //latest_times_of_vertex
    for (int i = vertex_num - 2; i > 0; --i) {
        auto idx = topo_seq[i];
        for (unsigned j = 0; j < vertex_num; ++j) {
            if (adjacency_matrix[idx][j]) {
                if (ltv[j] - adjacency_matrix[idx][j] < ltv[idx]) {
                    ltv[idx] = ltv[j] - adjacency_matrix[idx][j];
                }
            }
        }
    }
    lst_uint critical_path_vertexes;
    for (auto& i: topo_seq) {
        for (unsigned j = 1; j < vertex_num; ++j) {
            if (adjacency_matrix[i][j]) {
                auto ete = etv[i];
                auto lte = ltv[j] - adjacency_matrix[i][j];
                if (ete == lte) {
                    critical_path_vertexes.push_back(i);
                    critical_path_vertexes.push_back(j);
                } 
            }
        }
    }
    lst_uint diff;
    for (unsigned i = 0; i < vertex_num; ++i) {
        if (find(critical_path_vertexes.begin(), critical_path_vertexes.end(), i) == critical_path_vertexes.end()) {
            diff.push_back(i);
        }
    }
    Graph g1(*this);
    g1.delete_vertex(diff);
    return g1;
}

PYBIND11_MODULE(graph, m) {
    
    m.def("compare_by_last_element", &compare_by_last_element);
    
    py::class_<AdjacencyListNode>(m, "AdjacencyListNode")
    .def("get_idx", &AdjacencyListNode::get_idx)
    .def("get_next", &AdjacencyListNode::get_next);

    py::class_<Graph>(m, "Graph")
    .def(py::init<>())
    .def(py::init<const string&, const bool&, const unsigned&>())
    .def(py::init<const string&, const bool&, const unsigned&, const vector<string>&, const vector<vector<unsigned>>&>())
    .def("info", &Graph::info)
    .def("get_id", &Graph::get_id)
    .def("get_vertex_list", &Graph::get_vertex_list)
    .def("get_adj_matrix", &Graph::get_adj_matrix)
    .def("add_vertex", &Graph::add_vertex)
    .def("add_edge", &Graph::add_edge)
    .def("delete_vertex", &Graph::delete_vertex)
    .def("delete_edge", &Graph::delete_edge)
    .def("get_edgeset_array", &Graph::get_edgeset_array)
    .def("get_adj_list", &Graph::get_adj_list)
    .def("get_in_degree", &Graph::get_in_degree)
    .def("get_out_degree", &Graph::get_out_degree)
    .def("is_directed", &Graph::is_directed)
    .def("is_connected", &Graph::is_connected)
    .def("is_strongly_connected", &Graph::is_strongly_connected)
    .def("is_cyclic", &Graph::is_cyclic)
    .def("random_init", &Graph::random_init)
    .def("traverse", &Graph::traverse)
    .def("shortest_path", &Graph::shortest_path)
    .def("minimum_spanning_tree", &Graph::minimum_spanning_tree)
    .def("topological_sort", &Graph::topological_sort)
    .def("critical_path", &Graph::critical_path);
    
}
