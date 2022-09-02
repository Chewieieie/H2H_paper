#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <string>
#include "stdlib.h"
#include <time.h>

using namespace std;



struct Graph {
    // basic structure
    int vertex_num;
    int **g;

    // 使用容器实现
    unordered_map<int, unordered_map<int, int>> dis_map;

    // used for marked removed vertex
    set<int> removed_vertex_set;

};

struct TreeNode {
    int root_node;  // head of the tree node
    vector<int> vertex_list;
    vector<TreeNode*> children;
};


// X(v) is a star that contains not only the vertices but also the edges from v to u in N(v, H), H is the graph.
// see pp.716
struct StarTreeNode {
    int root_node;
    vector<int> vertex_list;
    vector<StarTreeNode*> children;
    StarTreeNode* parent;
    unordered_map<int, int> edges; // distance from root_node to the u in vertex_list
    vector<int> phi;

    vector<int> pos;  // pos array, see pp.717
    vector<int> anc;  // anc array, see pp.714
    vector<int> dist; // dist array, see pp.717
};


struct DPGraph {
    int vertex_num;
    vector<int> vertex_list;
    unordered_map<int, unordered_map<int, int>> dist_map;
};


// 新建一个树节点（bag）
TreeNode* new_tree_node(vector<int> *vertex_list) {
    auto t = new TreeNode;
    t->root_node = (*vertex_list)[0];
    for (int i : (*vertex_list)) {
        t->vertex_list.push_back(i);
    }

    return t;
}


// 新建一个星型树节点
StarTreeNode* new_star_tree_node_new(vector<int> *vertex_list, Graph* H) {
    auto t = new StarTreeNode;
    t->root_node = (*vertex_list)[0];
    for (int i : (*vertex_list)) {
        t->vertex_list.push_back(i);
    }

    for (int i=1; i < t->vertex_list.size(); i++) {

        int dis = H->dis_map.count(t->root_node) > 0 ? H->dis_map[t->root_node][t->vertex_list[i]] : 0;
        if (dis > 0) {
            t->edges.insert({t->vertex_list[i], dis});
        }
    }

    t->parent = nullptr;
    return t;
}


// 新建一个星型树节点
StarTreeNode* new_star_tree_node(vector<int> *vertex_list, Graph* H) {
    auto t = new StarTreeNode;
    t->root_node = (*vertex_list)[0];
    for (int i : (*vertex_list)) {
        t->vertex_list.push_back(i);
    }

    for (int i=1; i < t->vertex_list.size(); i++) {
        int dis = H->g[t->root_node][t->vertex_list[i]];
        if (dis > 0) {
            t->edges.insert({t->vertex_list[i], dis});
        }
    }

    t->parent = nullptr;
    return t;
}

void print_tree_node(TreeNode *p) {

    printf("node bag: ");
    for (auto item : p->vertex_list) {
        printf("%d ", item+1);
    }
    printf("\n");

    printf("children: ");
    for (auto item : p->children) {
        printf("%d ", item->root_node+1);
    }
    printf("\n");

}

void print_star_tree_node(StarTreeNode *p) {
    printf("root node: %d\n", p->root_node+1);
    printf("node bag: ");
    for (auto item : p->vertex_list) {
        printf("%d ", item+1);
    }
    printf("\n");

    printf("children: ");
    for (auto item : p->children) {
        printf("%d ", item->root_node+1);
    }
    printf("\n");

    if (p->parent == nullptr){
        printf("parent: null \n");
    } else{
        printf("parent: %d \n", p->parent->root_node+1);
    }

    printf("phi: ");
    for (auto item : p->phi){
        printf("%d ", item);
    }
    printf("\n");


    printf("anc: ");
    for (auto item : p->anc){
        printf("%d ", item + 1);
    }
    printf("\n");


    printf("dist: ");
    for (auto item : p->dist){
        printf("%d ", item);
    }
    printf("\n");
    printf("******************************* \n");
}


void print_graph_new(Graph *G) {
    printf("total vertex num: %d \n", G->vertex_num);
    printf("******************************* \n");

    for (int i=0; i < G->vertex_num; i++) {
        bool flag = false;
        if (G->dis_map.count(i) > 0) {
            for (int j=0; j < G->vertex_num; j++) {

                if (G->dis_map[i].count(j) > 0) {
                    if (flag) {
                        printf("        |-");
                    } else{
                        printf("vertex: %d", i+1);
                    }

                    printf("-----%d-----> %d\n", G->dis_map[i][j], j+1);
                    flag = true;
                }
            }

            if (!flag) {
                printf("vertex: %d\n", i+1);
            }

            printf("******************************* \n");

        }
    }

}


void print_graph(Graph *G) {
    printf("total vertex num: %d \n", G->vertex_num);
    printf("******************************* \n");

    for (int i=0; i < G->vertex_num; i++) {
        bool flag = false;
        for (int j=0; j < G->vertex_num; j++) {
            if (G->g[i][j] > 0) {
                if (flag) {
                    printf("        |-");
                } else{
                    printf("vertex: %d", i+1);
                }

                printf("-----%d-----> %d\n", G->g[i][j], j+1);
                flag = true;
            }
        }

        if (!flag) {
            printf("vertex: %d\n", i+1);
        }

        printf("******************************* \n");
    }


}


// 获取图中度最小的顶点
int get_vertex_with_smallest_degree_new(Graph* H) {
    int num = H->vertex_num;
    int min_degree = num;
    int min_vertex = -1;

    for (int i = 0; i < num; i++) {
        if (H->removed_vertex_set.count(i) > 0) {
            continue;
        }

        int temp = 0;

        temp += int(H->dis_map.at(i).size());

        if (temp < min_degree) {
            min_vertex = i;
            min_degree = temp;
        }

    }
    return min_vertex;
}


// 获取图中度最小的顶点
int get_vertex_with_smallest_degree(Graph* H) {
    int num = H->vertex_num;
    int min_degree = num;
    int min_vertex = -1;

    for (int i=0; i < num; i++) {

        if (H->removed_vertex_set.count(i) > 0) {
            continue;
        }

        int temp = 0;
        for (int j=0; j < num; j++) {
            if (H->g[i][j] > 0){
                temp++;
            }
        }

        if (temp < min_degree) {
            min_vertex = i;
            min_degree = temp;
        }

    }

    return min_vertex;
}

// 按照事先的排序获取图的下一个处理节点
int get_vertex_rank(int index) {
    const int rank[20] = {12, 9, 10, 3, 11, 4, 5, 7, 8, 6, 2, 1, 15, 19, 20, 17, 16, 18, 13, 14};

    return rank[index] - 1;
}


// 获取一个顶点在图中的所有邻居
vector<int>* get_vertex_neighbors_in_graph_new(int v, Graph *G) {
    auto list = new vector<int>;

    if (G->dis_map.count(v) <= 0) {
        return list;
    }

    for (auto kv : G->dis_map[v]) {
        list->push_back(kv.second);
    }

    return list;
}

// 获取一个顶点在图中的所有邻居
vector<int>* get_vertex_neighbors_in_graph(int v, Graph *G) {
    auto list = new vector<int>;
    for (int i=0; i < G->vertex_num; i++) {
        if (G->g[v][i] > 0) {
            list->push_back(i);
        }
    }
    return list;
}

// 图的深拷贝
Graph* deep_copy_graph_new(Graph *G) {
    auto* H = new Graph;
    int vertex_num = G->vertex_num;

    unordered_map<int, unordered_map<int, int>> dis_map;
    int counter = 0;

    for (const auto& kv : G->dis_map) {
        dis_map[kv.first] = unordered_map<int, int>();
        for (auto item : kv.second) {
            dis_map[kv.first][item.first] = item.second;
        }
    }

    H->vertex_num = vertex_num;
    H->dis_map = dis_map;

    return H;
}


// 图的深拷贝
Graph* deep_copy_graph(Graph *G) {
    Graph* H = new Graph;
    int vertex_num = G->vertex_num;
    
    
    int **arr = new int*[vertex_num];
    for (int i=0; i < vertex_num; i++) {
        arr[i] = new int[vertex_num];
    }
    
    for (int i=0; i<vertex_num; i++) {
        for (int j=0; j<vertex_num; j++) {
            arr[i][j] = G->g[i][j];
        }
    }

    H->vertex_num = vertex_num;
    H->g = arr;
    
    return H;
}

// 连接图中的顶点v的所有邻居顶点
void connect_vertex_neighbors_in_graph_new(int v, Graph *G) {
    const int DEFAULT_DISTANCE = 999;
    vector<int> *neighbors = get_vertex_neighbors_in_graph_new(v, G);

    for (int i=0; i < neighbors->size(); i++) {
        for (int j=0; j < neighbors->size(); j++) {

            if (i == j) {
                continue;
            }

            int &from = neighbors->at(i);
            int &to = neighbors->at(j);

            if (G->dis_map.count(from) <= 0) {
                G->dis_map[from] = unordered_map<int, int>();
            }

            if (G->dis_map.count(to) <= 0) {
                G->dis_map[to] = unordered_map<int, int>();
            }

            G->dis_map[from][to] = DEFAULT_DISTANCE;
            G->dis_map[to][from] = DEFAULT_DISTANCE;

        }
    }

}


// 连接图中的顶点v的所有邻居顶点
void connect_vertex_neighbors_in_graph(int v, Graph *G) {
    const int DEFAULT_DISTANCE = 999;
    vector<int> *neighbors = get_vertex_neighbors_in_graph(v, G);
    
    for (int i=0; i < neighbors->size(); i++) {
        for (int j=0; j < neighbors->size(); j++) {
            
            if (i == j) {
                continue;
            }

            int &from = neighbors->at(i);
            int &to = neighbors->at(j);

            if (G->g[from][to] <= 0) {
                G->g[from][to] = DEFAULT_DISTANCE;
                G->g[to][from] = DEFAULT_DISTANCE;
            }

        }
    }

}


// 移除图中的一个顶点以及其相邻的边
void remove_vertex_and_adjacent_edges_in_graph_new(int v, Graph *G) {
    G->removed_vertex_set.insert(v);

    for (int i=0; i < G->vertex_num; i++) {

        if (G->dis_map.count(v) > 0) {
            if (G->dis_map[v].count(i) > 0) {
                G->dis_map[v].erase(i);
            }

            if (G->dis_map[v].empty()) {
                G->dis_map.erase(v);
            }
        }

        if (G->dis_map.count(i) > 0) {
            if (G->dis_map[i].count(v) > 0) {
                G->dis_map[i].erase(v);
            }

            if (G->dis_map[i].empty()) {
                G->dis_map.erase(i);
            }
        }

    }

}


// 移除图中的一个顶点以及其相邻的边
void remove_vertex_and_adjacent_edges_in_graph(int v, Graph *G) {
    G->removed_vertex_set.insert(v);

    for (int i=0; i < G->vertex_num; i++) {
        if (G->g[v][i] > 0) {
            G->g[v][i] = 0;
            G->g[i][v] = 0;
        }
    }

}


void dp_vertex_elimination_new(int v, Graph* G, vector<int>* neighbors) {
    for (int i=0; i < neighbors->size(); i++) {
        for (int j=0; j < neighbors->size(); j++) {

            if (i == j) {
                continue;
            }

            int n_i = neighbors->at(i);
            int n_j = neighbors->at(j);

            if (G->dis_map.count(n_i) <= 0 || G->dis_map[n_i][n_j] <= 0) {

                int i_to_v = G->dis_map.count(n_i) > 0 ? G->dis_map[n_i][v] : 0;
                int v_to_j = G->dis_map.count(v) > 0 ? G->dis_map[v][n_j] : 0;

                if(G->dis_map.count(n_i) <= 0) {
                    G->dis_map[n_i] = unordered_map<int, int>();
                }
                if (G->dis_map.count(n_j) <= 0) {
                    G->dis_map[n_j] = unordered_map<int, int>();
                }

                G->dis_map[n_i][n_j] = i_to_v + v_to_j;
                G->dis_map[n_j][n_i] = G->dis_map[n_i][n_j];
            } else {

                if (G->dis_map.count(n_i) > 0 && G->dis_map[n_i][n_j] > 0) {
                    int v_to_j = G->dis_map.count(v) > 0 ? G->dis_map[v][n_j] : 0;
                    if (G->dis_map[n_i][n_j] > G->dis_map[n_i][v] + v_to_j) {
                        if(G->dis_map.count(n_i) <= 0) {
                            G->dis_map[n_i] = unordered_map<int, int>();
                        }
                        if (G->dis_map.count(n_j) <= 0) {
                            G->dis_map[n_j] = unordered_map<int, int>();
                        }

                        G->dis_map[n_i][n_j] = G->dis_map[n_i][v] + v_to_j;
                        G->dis_map[n_j][n_i] = G->dis_map[n_i][n_j];
                    }
                }
            }

        }
    }

    // remove v from the graph
    G->removed_vertex_set.insert(v);
    G->dis_map.erase(v);
    for (int i=0; i < G->vertex_num; i++) {
        if (G->dis_map.count(i) > 0) {
            G->dis_map[i].erase(v);
        }
    }
}

void dp_vertex_elimination(int v, Graph* G) {
    vector<int> *neighbors = get_vertex_neighbors_in_graph(v, G);


    for (int i=0; i < neighbors->size(); i++) {
        for (int j=0; j < neighbors->size(); j++) {

            if (i == j) {
                continue;
            }

            int n_i = neighbors->at(i);
            int n_j = neighbors->at(j);

            // the is no edge connect i and j, then we insert edge i -> j = i -> v + v -> j
            if (G->g[n_i][n_j] <= 0) {
                G->g[n_i][n_j] = G->g[n_i][v] + G->g[v][n_j];
                G->g[n_j][n_i] = G->g[n_i][v] + G->g[v][n_j];
            } else if (G->g[n_i][n_j] > 0 && G->g[n_i][n_j] > G->g[n_i][v] + G->g[v][n_j]) {
                G->g[n_i][n_j] = G->g[n_i][v] + G->g[v][n_j];
                G->g[n_j][n_i] = G->g[n_i][v] + G->g[v][n_j];
            }
        }
    }

    // remove v from the graph
    G->removed_vertex_set.insert(v);
    for (int i=0; i < G->vertex_num; i++) {
        if (G->g[v][i] > 0) {
            G->g[v][i] = 0;
            G->g[i][v] = 0;
        }
    }

}

// 获取顶点列表中拥有最小消去序的顶点
int get_vertex_with_smallest_elimination_order(const int& v, const vector<int>& vertex_list, const unordered_map<int, int>& pi) {
    int min_order = INT32_MAX;
    int res = -1;

    for (auto item : vertex_list) {
        if (item != v && pi.count(item) > 0 && pi.at(item) < min_order) {
            min_order = pi.at(item);
            res = item;
        }
    }

    return res;
}

Graph* get_test_road_network_graph_new() {
    const unsigned int vertex_num = 20;

    unordered_map<int, unordered_map<int, int>> dis_map;

    for (int i=0; i < vertex_num; i++) {
        dis_map[i] = unordered_map<int, int>();
    }

    // 1
    dis_map[0][1] = 2;
    dis_map[0][2] = 1;
    dis_map[0][12] = 2;

    // 2
    dis_map[1][0] = 2;
    dis_map[1][4] = 1;
    dis_map[1][13] = 2;

    // 3
    dis_map[2][0] = 1;
    dis_map[2][3] = 1;
    dis_map[2][4] = 2;

    // 4
    dis_map[3][2] = 1;
    dis_map[3][4] = 1;
    dis_map[3][5] = 2;

    // 5
    dis_map[4][1] = 1;
    dis_map[4][2] = 2;
    dis_map[4][3] = 1;
    dis_map[4][5] = 3;

    //6
    dis_map[5][3] = 2;
    dis_map[5][4] = 3;
    dis_map[5][7] = 4;
    dis_map[5][8] = 2;

    // 7
    dis_map[6][7] = 2;
    dis_map[6][9] = 2;
    dis_map[6][13] = 2;

    // 8
    dis_map[7][6] = 2;
    dis_map[7][10] = 2;
    dis_map[7][13] = 2;

    // 9
    dis_map[8][5] = 2;
    dis_map[8][11] = 2;

    // 10
    dis_map[9][6] = 2;
    dis_map[9][10] = 2;

    // 11
    dis_map[10][7] = 2;
    dis_map[10][9] = 2;
    dis_map[10][11] = 1;

    // 12
    dis_map[11][8] = 2;
    dis_map[11][10] = 1;

    // 13
    dis_map[12][0] = 2;
    dis_map[12][13] = 1;
    dis_map[12][17] = 1;

    // 14
    dis_map[13][1] = 2;
    dis_map[13][12] = 1;
    dis_map[13][6] = 2;
    dis_map[13][7] = 2;
    dis_map[13][14] = 1;
    dis_map[13][15] = 3;

    // 15
    dis_map[14][13] = 1;
    dis_map[14][15] = 2;

    // 16
    dis_map[15][14] = 2;
    dis_map[15][13] = 3;
    dis_map[15][16] = 2;
    dis_map[15][19] = 2;

    // 17
    dis_map[16][15] = 2;
    dis_map[16][17] = 2;
    dis_map[16][18] = 2;

    // 18
    dis_map[17][12] = 1;
    dis_map[17][16] = 2;

    // 19
    dis_map[18][16] = 2;
    dis_map[18][19] = 1;

    // 20
    dis_map[19][15] = 2;
    dis_map[19][18] = 1;

    auto *G = new Graph;
    
    G->vertex_num = vertex_num;
    G->dis_map = dis_map;
    
    return G;
}



Graph* get_test_road_network_graph() {
    const unsigned int vertex_num = 20;
    int **arr = new int*[vertex_num];
    for (int i=0; i<vertex_num; i++) {
        arr[i] = new int[vertex_num];
        for (int j=0; j<vertex_num; j++) {
            arr[i][j] = 0;
        }
    }

    Graph *G = new Graph;

    // 1
    arr[0][1] = 2;
    arr[0][2] = 1;
    arr[0][12] = 2;

    // 2
    arr[1][0] = 2;
    arr[1][4] = 1;
    arr[1][13] = 2;

    // 3
    arr[2][0] = 1;
    arr[2][3] = 1;
    arr[2][4] = 2;

    // 4
    arr[3][2] = 1;
    arr[3][4] = 1;
    arr[3][5] = 2;

    // 5
    arr[4][1] = 1;
    arr[4][2] = 2;
    arr[4][3] = 1;
    arr[4][5] = 3;

    //6
    arr[5][3] = 2;
    arr[5][4] = 3;
    arr[5][7] = 4;
    arr[5][8] = 2;

    // 7
    arr[6][7] = 2;
    arr[6][9] = 2;
    arr[6][13] = 2;

    // 8
    arr[7][6] = 2;
    arr[7][10] = 2;
    arr[7][13] = 2;

    // 9
    arr[8][5] = 2;
    arr[8][11] = 2;

    // 10
    arr[9][6] = 2;
    arr[9][10] = 2;

    // 11
    arr[10][7] = 2;
    arr[10][9] = 2;
    arr[10][11] = 1;

    // 12
    arr[11][8] = 2;
    arr[11][10] = 1;

    // 13
    arr[12][0] = 2;
    arr[12][13] = 1;
    arr[12][17] = 1;

    // 14
    arr[13][1] = 2;
    arr[13][12] = 1;
    arr[13][6] = 2;
    arr[13][7] = 2;
    arr[13][14] = 1;
    arr[13][15] = 3;

    // 15
    arr[14][13] = 1;
    arr[14][15] = 2;

    // 16
    arr[15][14] = 2;
    arr[15][13] = 3;
    arr[15][16] = 2;
    arr[15][19] = 2;

    // 17
    arr[16][15] = 2;
    arr[16][17] = 2;
    arr[16][18] = 2;

    // 18
    arr[17][12] = 1;
    arr[17][16] = 2;

    // 19
    arr[18][16] = 2;
    arr[18][19] = 1;

    // 20
    arr[19][15] = 2;
    arr[19][18] = 1;

    G->vertex_num = vertex_num;
    G->g = arr;

    return G;
}

void sort_start_tree_node(StarTreeNode* t, unordered_map<int, int>* pi) {
    
    for (int i=0; i < t->vertex_list.size(); i++) {
        for (int j=i+1; j < t->vertex_list.size(); j++) {
            int pi_i = pi->at(i);
            int pi_j = pi->at(j);
            
            if (pi_i < pi_j) {
                int temp = t->vertex_list[i];
                t->vertex_list[i] = t->vertex_list[j];
                t->vertex_list[j] = temp;
            }
            
        }
    }
    
}

void compute_star_tree_node_phi(StarTreeNode* t) {
    for (int i=0; i<t->vertex_list.size(); i++) {
        if (t->root_node == t->vertex_list[i]) {
            t->phi.push_back(0);
        } else{
            t->phi.push_back(t->edges.at(t->vertex_list[i]));
        }

    }
}

void compute_anc_and_pos_array(StarTreeNode* node) {
    vector<int> anc;
    vector<int> pos;

    unordered_map<int, int> temp_map; //k: node index, v: index in anc.

    StarTreeNode *temp = node;
    while (temp != nullptr) {
        anc.insert(anc.begin(), temp->root_node);
        temp = temp->parent;
    }

    for (int i=0; i < anc.size(); i++) {
        temp_map.insert({anc[i], i});
    }

    for (auto item : node->vertex_list) {
        pos.push_back(temp_map.at(item));
    }

    node->anc = anc;
    node->pos = pos;

}


DPGraph* construct_dp_graph(StarTreeNode *node, unordered_map<int, StarTreeNode*> tree_node_map) {
    auto *dp_graph = new DPGraph;

    int vertex_num = node->anc.size();
    dp_graph->vertex_num = vertex_num;

    // 初始化dp_graph中所有的顶点
    for (int i=0; i < vertex_num; i++) {
        auto tree_node = tree_node_map.at(node->anc[i]);
        int v = tree_node->root_node;
        unordered_map<int, int> mp;
        dp_graph->dist_map.insert({v, mp});
        dp_graph->vertex_list.push_back(v);
    }


    for (int i=vertex_num-1; i>=0; i--) {
        auto tree_node = tree_node_map.at(node->anc[i]);
        int v = tree_node->root_node;

        for (auto item : tree_node->edges) {
            dp_graph->dist_map.at(v).insert({item.first, item.second});
            dp_graph->dist_map.at(item.first).insert({v, item.second});

        }

    }
    
    return dp_graph;

}

// 获取树的根
int find_tree_root(unordered_map<int, StarTreeNode*> tree_node_map) {
    for (auto &kv:tree_node_map) {
        StarTreeNode *node = kv.second;
        if (node->parent == nullptr) {
            return node->root_node;
        }
    }

    return -1;
}

// new Algorithm 3 DPTreeDecomposition(G(V, E))
unordered_map<int, StarTreeNode*> dp_tree_decomposition_new(Graph* G) {
    auto H = deep_copy_graph_new(G);

    unordered_map<int, StarTreeNode*> star_tree_node_map;

    // the elimination order
    unordered_map<int, int> pi;

    for (int i=0; i < H->vertex_num; i++) {

        if (i % 50 == 0) {
            printf("dp_tree_decomposition progress: %d \n", i);
            fflush(stdout);
        }

        // todo: 优化到O(1)
        int v = get_vertex_with_smallest_degree_new(H);

        vector<int> *neighbors = get_vertex_neighbors_in_graph_new(v, H);
        neighbors->insert(neighbors->begin(), v);

        StarTreeNode *Xv = new_star_tree_node_new(neighbors, H);
        star_tree_node_map.insert({Xv->root_node, Xv});

        // todo: 优化
        dp_vertex_elimination_new(v, H, neighbors);

        pi.insert({v, i});
    }


    for (int i=0; i < H->vertex_num; i++) {
        if (star_tree_node_map.count(i) > 0 && star_tree_node_map.at(i)->vertex_list.size() > 1) {
            StarTreeNode *X_v = star_tree_node_map.at(i);

            int u = get_vertex_with_smallest_elimination_order(i, X_v->vertex_list, pi);
            StarTreeNode *X_u = star_tree_node_map.at(u);
            X_u->children.push_back(X_v);
            X_v->parent = X_u;
        }
    }

    // sort each start tree node
    for (int i=0; i < H->vertex_num; i++) {
        sort_start_tree_node(star_tree_node_map.at(i), &pi);
        compute_star_tree_node_phi(star_tree_node_map.at(i));
    }

    return star_tree_node_map;
}

// Algorithm 3 DPTreeDecomposition(G(V, E))
unordered_map<int, StarTreeNode*> dp_tree_decomposition(Graph* G) {
    auto H = deep_copy_graph(G);
    unordered_map<int, StarTreeNode*> star_tree_node_map;

    // the elimination order
    unordered_map<int, int> pi;

    for (int i=0; i < H->vertex_num; i++) {
        int v = get_vertex_with_smallest_degree(H);
        vector<int> *neighbors = get_vertex_neighbors_in_graph(v, H);
        neighbors->insert(neighbors->begin(), v);

        StarTreeNode *Xv = new_star_tree_node(neighbors, H);
        star_tree_node_map.insert({Xv->root_node, Xv});

        dp_vertex_elimination(v, H);

        pi.insert({v, i});
    }

    for (int i=0; i < H->vertex_num; i++) {
        if (star_tree_node_map.count(i) > 0 && star_tree_node_map.at(i)->vertex_list.size() > 1) {
            StarTreeNode *X_v = star_tree_node_map.at(i);

            int u = get_vertex_with_smallest_elimination_order(i, X_v->vertex_list, pi);
            StarTreeNode *X_u = star_tree_node_map.at(u);
            X_u->children.push_back(X_v);
            X_v->parent = X_u;
        }
    }

    // sort each start tree node
    for (int i=0; i < H->vertex_num; i++) {
        sort_start_tree_node(star_tree_node_map.at(i), &pi);
        compute_star_tree_node_phi(star_tree_node_map.at(i));
    }

    return star_tree_node_map;
}

// New Algorithm 6 TreeDecomposition(G(V, E)) new
unordered_map<int, TreeNode*> tree_decomposition_new(Graph* G) {
    auto H = deep_copy_graph_new(G);
    unordered_map<int, TreeNode *> tree_node_map;

    // the elimination order
    unordered_map<int, int> pi;

    for (int i = 0; i < H->vertex_num; i++) {
        int v = get_vertex_with_smallest_degree_new(H);
        vector<int> *neighbors = get_vertex_neighbors_in_graph_new(v, H);
        neighbors->insert(neighbors->begin(), v);

        // create a node X(v) in TG
        TreeNode *Xv = new_tree_node(neighbors);
        tree_node_map.insert({Xv->root_node, Xv});

        // Add edges to H to make every pair of veritces in N(v, H) connected to each other in H;
        connect_vertex_neighbors_in_graph_new(v, H);

        // remove v and its adjacent edges from H
        remove_vertex_and_adjacent_edges_in_graph_new(v, H);

        pi.insert({v, i});
    }

    for (int i=0; i < H->vertex_num; i++) {
        if (tree_node_map.count(i) > 0 && tree_node_map.at(i)->vertex_list.size() > 1) {
            TreeNode *X_v = tree_node_map.at(i);
            
            // no need to modify
            int u = get_vertex_with_smallest_elimination_order(i, X_v->vertex_list, pi);
            TreeNode *X_u = tree_node_map.at(u);
            X_u->children.push_back(X_v);
        }
    }

    return tree_node_map;
    
}


// Algorithm 6 TreeDecomposition(G(V, E))
unordered_map<int, TreeNode*> tree_decomposition(Graph* G) {
    auto H = deep_copy_graph(G);
    unordered_map<int, TreeNode*> tree_node_map;

    // the elimination order
    unordered_map<int, int> pi;

    for (int i=0; i < H->vertex_num; i++) {
        int v = get_vertex_with_smallest_degree(H);
        // int v = get_vertex_rank(i);

        vector<int> *neighbors = get_vertex_neighbors_in_graph(v, H);
        neighbors->insert(neighbors->begin(), v);
        
        // create a node X(v) in TG
        TreeNode *Xv = new_tree_node(neighbors);

        tree_node_map.insert({Xv->root_node, Xv});


        // Add edges to H to make every pair of veritces in N(v, H) connected to each other in H;
        connect_vertex_neighbors_in_graph(v, H);

        // remove v and its adjacent edges from H
        remove_vertex_and_adjacent_edges_in_graph(v, H);

        pi.insert({v, i});
    }


    for (int i=0; i < H->vertex_num; i++) {
        if (tree_node_map.count(i) > 0 && tree_node_map.at(i)->vertex_list.size() > 1) {
            TreeNode *X_v = tree_node_map.at(i);

            int u = get_vertex_with_smallest_elimination_order(i, X_v->vertex_list, pi);
            TreeNode *X_u = tree_node_map.at(u);
            X_u->children.push_back(X_v);
        }
    }

    return tree_node_map;
}

unordered_map<int, int> dijkstra_for_dp_graph(DPGraph *G, int src) {
    unordered_map<int, int> dist;
    set<int> S;
    set<int> Q;

    for (auto v : G->vertex_list) {
        if (src == v) {
            dist.insert({v, 0});
        } else {
            // 如果 src -> v之间存在路径
            if (G->dist_map.count(src) > 0 && G->dist_map.at(src).count(v) > 0) {
                dist.insert({v, G->dist_map.at(src).at(v)});
            } else {
                dist.insert({v, 88888});
            }

        }

        if (src != v) {
            Q.insert(v);
        }
    }

    while(!Q.empty()) {

        // 找到最小距离的那个点
        int min_dist = INT32_MAX;
        int u = -1;
        for (auto item : Q) {
            if (dist.at(item) < min_dist) {
                min_dist = dist.at(item);
                u = item;
            }
        }
        Q.erase(u);
        S.insert(u);

        // 获取u相邻的节点
        set<int> neighbour;

        for (auto item : G->vertex_list) {
            if (Q.count(item) > 0 && G->dist_map.at(u).count(item) > 0) {
                neighbour.insert(item);
            }
        }

        for (auto item : neighbour) {
            if (dist[u] + G->dist_map[u][item] < dist[item]) {
                dist[item] = dist[u] + G->dist_map[u][item];

            }
        }

    }

    return dist;
}


unordered_map<int, int> dijkstra(Graph *G, int src) {

    unordered_map<int, int> dist;

    set<int> S;
    set<int> Q;

    S.insert(src);
    for (int i=0; i<G->vertex_num; i++) {
        // init dist
        if (src == i) {
            dist.insert({i, 0});
        } else{
            if (G->g[src][i] <= 0) {
                dist.insert({i, INT32_MAX});
            } else {
                dist.insert({i, G->g[src][i]});
            }
        }

        // init Q
        if (src != i) {
            Q.insert(i);
        }
    }

    while (!Q.empty()) {
        // find the vertex x in Q, such that s -> x is min in dist
        int min_dist = INT32_MAX;
        int u = -1;
        for (auto item : Q) {
            if (dist.at(item) < min_dist) {
                min_dist = dist.at(item);
                u = item;
            }
        }
        Q.erase(u);
        S.insert(u);

        // 获取u相邻的节点
        set<int> neighbour;
        for (int i=0; i < G->vertex_num; i++) {
            if (G->g[u][i] > 0 && Q.count(i)) {
                neighbour.insert(i);
            }
        }

        for (auto item : neighbour) {
            if (dist[u] + G->g[u][item] < dist[item]) {
                dist[item] = dist[u] + G->g[u][item];
            }
        }
    }

    return dist;
}



// New Algorithm 4 H2H-Naive(G)
unordered_map<int, StarTreeNode*> H2H_naive_new(Graph *G) {
    unordered_map<int, StarTreeNode*> star_tree_node_map = dp_tree_decomposition_new(G);


    for (int i=0; i<G->vertex_num; i++) {

        if (i % 10 == 0) {
            printf("build H2H navie index progress: %d", i);
        }

        StarTreeNode* node = star_tree_node_map.at(i);
        compute_anc_and_pos_array(node);
        auto dp_graph = construct_dp_graph(node, star_tree_node_map);

        const unordered_map<int, int> &min_dist = dijkstra_for_dp_graph(dp_graph, i);

        for (int j=0; j < node->anc.size(); j++) {
            node->dist.push_back(min_dist.at(node->anc[j]));
        }
    }
    
    return star_tree_node_map;
    
}


// Algorithm 4 H2H-Naive(G)
unordered_map<int, StarTreeNode*> H2H_naive(Graph *G) {
    unordered_map<int, StarTreeNode*> star_tree_node_map = dp_tree_decomposition(G);

    int tree_root = find_tree_root(star_tree_node_map);

    for (int i=0; i<G->vertex_num; i++) {
        StarTreeNode* node = star_tree_node_map.at(i);
        compute_anc_and_pos_array(node);
        auto dp_graph = construct_dp_graph(node, star_tree_node_map);

        const unordered_map<int, int> &min_dist = dijkstra_for_dp_graph(dp_graph, i);

        for (int j=0; j < node->anc.size(); j++) {
            node->dist.push_back(min_dist.at(node->anc[j]));
        }
    }

    return star_tree_node_map;
    
}

// Algorithm 5 H2H(G)
void H2H(Graph *G) {



}


void test_algorithm_6_tree_decomposition() {
    Graph *graph = get_test_road_network_graph();
    print_graph(graph);
    const unordered_map<int, TreeNode *> &d = tree_decomposition(graph);

    for (int i=0; i<20; i++) {
        TreeNode *p = d.at(i);
        print_tree_node(p);
        printf("xxxxxxxxxxxxxxxxxxxxxxx \n");
    }
}

void test_dijkstra() {
    const unsigned int vertex_num = 20;
    int **arr = new int*[vertex_num];
    for (int i=0; i<vertex_num; i++) {
        arr[i] = new int[vertex_num];
        for (int j=0; j<vertex_num; j++) {
            arr[i][j] = 9999999;
        }
    }


    Graph *G = new Graph;

    // s
    arr[0][0] = 0;
    arr[0][1] = 5;
    arr[0][3] = 10;

    // y
    arr[1][1] = 0;
    arr[1][3] = 3;
    arr[1][4] = 9;
    arr[1][2] = 2;

    // z
    arr[2][2] = 0;
    arr[2][4] = 6;
    arr[2][0] = 7;

    // t
    arr[3][3] = 0;
    arr[3][1] = 2;
    arr[3][4] = 1;

    // x
    arr[4][4] = 0;
    arr[4][2] = 4;

    G->g = arr;
    G->vertex_num = 5;



    const unordered_map<int, int> &map = dijkstra(G, 0);

    for (auto kv : map) {
        printf("s ----> %d: %d \n", kv.first, kv.second);
    }



}

void test_algorithm_3_dp_tree_decomposition() {
    Graph *graph = get_test_road_network_graph();
    const unordered_map<int, StarTreeNode*> &d = dp_tree_decomposition(graph);

    for (int i=0; i<20; i++) {
        StarTreeNode *p = d.at(i);
        print_star_tree_node(p);
        printf("xxxxxxxxxxxxxxxxxxxxxxx \n");
    }
}

void test_algorithm_4_H2H_naive() {
    printf("^^^^^^^^^^^^test_algorithm_4_H2H_naive^^^^^^^^^^^^ \n");
    Graph *graph = get_test_road_network_graph();
    const unordered_map<int, StarTreeNode *> &map = H2H_naive(graph);

    for (auto kv : map) {
        print_star_tree_node(kv.second);
    }
    printf("^^^^^^^^^^^^test_algorithm_4_H2H_naive^^^^^^^^^^^^ \n");

}


void test_algorithm_4_H2H_naive_new() {
    printf("^^^^^^^^^^^^test_algorithm_4_H2H_naive_new^^^^^^^^^^^^ \n");

    Graph *graph = get_test_road_network_graph_new();
    const unordered_map<int, StarTreeNode *> &map = H2H_naive_new(graph);

    for (auto kv : map) {
        print_star_tree_node(kv.second);
    }

    printf("^^^^^^^^^^^^test_algorithm_4_H2H_naive_new^^^^^^^^^^^^ \n");

}


vector<int> split_one_line(string s){
    int pre = s.find(' ');
    int post = s.find(' ', pre+1);
    vector<int> res;

    while(post > 0) {
        int v = atoi(s.substr(pre, post - pre + 1).c_str());
        res.push_back(v);

        pre = post;
        post = s.find(' ', pre+1);
    }

    res.push_back(atoi(s.substr(pre+1, s.length()-pre).c_str()));
    return res;

}

Graph* read_data_from_file() {
    ifstream infile;
    infile.open("/Users/chenhao/cpp-proj/Hierarchy-2-Hop-Labeling/USA-road-d.NY.gr", ios::in);

    string buffer;
    auto* res = new Graph;


    res->vertex_num = 264346;

    int counter = 0;
    int edge_num = 0;
    while (getline(infile, buffer)) {
        counter++;
        if(counter <= 8) {
            // 前8行都是一些无用的信息，这里直接跳过不处理
            continue;
        }
        const vector<int> &nums = split_one_line(buffer);

        if(res->dis_map.count(nums[0]-1) <= 0) {
            res->dis_map[nums[0]-1] = unordered_map<int, int>();
        }
        edge_num++;
        res->dis_map[nums[0]-1][nums[1]-1] = nums[2];
    }

    printf("total edeges are: %d \n", edge_num);
    printf("total vertexs are: %lu \n", res->dis_map.size());

    return res;
}

void benchmark_build_H2H_naive_index() {
    Graph *pGraph = read_data_from_file();
    H2H_naive_new(pGraph);
}



Graph process_data(){
    ifstream infile;
    infile.open("/Users/chenhao/cpp-proj/Hierarchy-2-Hop-Labeling/USA-road-d.NY.gr", ios::in);

    string buffer;

    Graph res;
    const unsigned int vertex_num = 1000;
    res.vertex_num = vertex_num;
    int **arr = new int*[vertex_num];
    for (int i=0; i<vertex_num; i++) {
        arr[i] = new int[vertex_num];
        for (int j=0; j<vertex_num; j++) {
            arr[i][j] = 0;
        }
    }

    res.g = arr;
    
    int counter = 0;
    while (getline(infile, buffer)) {
        counter++;
        if(counter <= 8) {
            continue;
        }
        const vector<int> &nums = split_one_line(buffer);
        printf("address in process_data: %p \n", &nums);

        if (counter == 9) {
            break;
        }

        arr[nums[0]-1][nums[1]-1] = nums[2];

        if(counter > 2016) {
            cout << "current line is" << nums[0] << " " << nums[1] << " " << nums[2];
            break;
        }
        
    }
    
    return res;

}


int main() {
    benchmark_build_H2H_naive_index();
}
