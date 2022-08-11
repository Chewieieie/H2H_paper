#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>

using namespace std;



struct Graph {
    // basic structure
    int vertex_num;
    int **g;

    // used for marked removed vertex
    set<int> removed_vertex_set;

};

struct TreeNode {
    int root_node;
    vector<int> vertex_list;
    vector<TreeNode*> children;
};


// X(v) is a star that contains not only the vertices but also the edges from v to u in N(v, H), H is the graph.
// see pp.716
struct StarTreeNode {
    int root_node;
    vector<int> vertex_list;
    vector<StarTreeNode*> children;
    unordered_map<int, int> edges; // distance from root_node to the u in vertex_list
    vector<int> phi;

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

    printf("phi: ");
    for (auto item : p->phi){
        printf("%d ", item);
    }
    printf("\n");
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
        printf("******************************* \n");
    }


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
void remove_vertex_and_adjacent_edges_in_graph(int v, Graph *G) {
    G->removed_vertex_set.insert(v);

    for (int i=0; i < G->vertex_num; i++) {
        if (G->g[v][i] > 0) {
            G->g[v][i] = 0;
            G->g[i][v] = 0;
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

            // the is no edge connect i and j, then we insert edge i -> j = i -> v + v -> j
            if (G->g[i][j] <= 0) {
                G->g[i][j] = G->g[i][v] + G->g[v][j];
                G->g[j][i] = G->g[i][v] + G->g[v][j];
            } else if (G->g[i][j] > 0 && G->g[i][j] > G->g[i][v] + G->g[v][j]) {
                G->g[i][j] = G->g[i][v] + G->g[v][j];
                G->g[j][i] = G->g[i][v] + G->g[v][j];
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
        }
    }

    // sort each start tree node
    for (int i=0; i < H->vertex_num; i++) {
        sort_start_tree_node(star_tree_node_map.at(i), &pi);
        compute_star_tree_node_phi(star_tree_node_map.at(i));
    }

    return star_tree_node_map;
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

void test_algorithm_3_dp_tree_decomposition() {
    Graph *graph = get_test_road_network_graph();
    const unordered_map<int, StarTreeNode*> &d = dp_tree_decomposition(graph);

    for (int i=0; i<20; i++) {
        StarTreeNode *p = d.at(i);
        print_star_tree_node(p);
        printf("xxxxxxxxxxxxxxxxxxxxxxx \n");
    }
}


int main() {
    Graph *graph = get_test_road_network_graph();
    const unordered_map<int, StarTreeNode*> &d = dp_tree_decomposition(graph);

    for (int i=0; i<20; i++) {
        StarTreeNode *p = d.at(i);
        print_star_tree_node(p);
        printf("xxxxxxxxxxxxxxxxxxxxxxx \n");
    }

}
