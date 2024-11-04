#include <iostream>
#include <vector>
#include <queue>

using namespace std;

/*
 * 拓扑排序算法：
 * 拓扑排序是针对有向无环图(Directed Acyclic Graph, DAG)的排序算法，它能够将图中的顶点排序，使得对于每一条 **有向边(u, v)**，顶点u在顶点v之前。
 * 卡恩(Kahn)算法是图的拓扑排序(Topological sorting)算法，它是基于队列实现的，类似于《宽度优先搜索(BFS)》->
 * Kahn 算法的核心思想是找到图中 没有前驱（入度为0）的节点，并依次将其加入拓扑排序的结果中。然后，移除这个节点及其边，更新其后续节点的入度，重复这一过程，直到所有节点都被处理完毕。
 * 需要的数据结构：1. 邻接表；2. 入度表；
 */
class Graph {
private:
    int vertices;  // 顶点数量
    vector<vector<int> > adjList;  // 邻接表
    vector<int> inDegree;  // 入度表

public:
    // 构造函数
    Graph(int v) : vertices(v), adjList(v), inDegree(v, 0) {}

    // 添加边，相当于遍历一次图获得入度
    void addEdge(int u, int v) {
        adjList[u].push_back(v);
        inDegree[v]++;  // 更新入度
    }

    // 检查是否有环
    bool hasCycle() {
        vector<int> tempInDegree = inDegree;  // 拷贝入度表
        queue<int> q;

        // 将所有入度为 0 的节点加入队列
        for (int i = 0; i < vertices; ++i) {
            if (tempInDegree[i] == 0) {
                q.push(i);
            }
        }

        int count = 0;  // 记录已排序的节点数

        while (!q.empty()) {
            int node = q.front();
            q.pop();
            count++;

            // 遍历相邻节点并减少入度
            for (int neighbor : adjList[node]) {
                if (--tempInDegree[neighbor] == 0) {
                    q.push(neighbor);
                }
            }
        }

        // 如果排序节点数等于图的顶点数，则无环；否则有环
        return count != vertices;
    }

    // 执行拓扑排序
    // 这个函数其实有点冗余，可以和上述函数合并
    vector<int> topologicalSort() {
        vector<int> result;
        if (hasCycle()) {
            cout << "Graph contains a cycle. Topological sorting is not possible." << endl;
            return result;  // 返回空结果
        }

        vector<int> tempInDegree = inDegree;  // 使用临时入度表
        queue<int> q;

        // 将入度为 0 的节点加入队列
        for (int i = 0; i < vertices; ++i) {
            if (tempInDegree[i] == 0) {
                q.push(i);
            }
        }

        while (!q.empty()) {
            int node = q.front();
            q.pop();
            result.push_back(node);

            // 减少相邻节点的入度
            for (int neighbor : adjList[node]) {
                if (--tempInDegree[neighbor] == 0) {
                    q.push(neighbor);
                }
            }
        }

        return result;
    }

    // 打印拓扑排序结果
    // 因为这里 vector 不能直接打印
    void printTopologicalSort() {
        vector<int> sortedOrder = topologicalSort();

        if (!sortedOrder.empty()) {
            cout << "Topological Sort: ";
            for (int node : sortedOrder) {
                cout << node << " ";
            }
            cout << endl;
        }
    }
};

int main() {
    int vertices, edges;
    cout << "Enter the number of vertices: ";
    cin >> vertices;
    cout << "Enter the number of edges: ";
    cin >> edges;

    Graph graph(vertices);

    cout << "Enter the edges (u v):" << endl;
    for (int i = 0; i < edges; ++i) {
        int u, v;
        cin >> u >> v;
        graph.addEdge(u, v);
    }

    // 验证是否有环并打印拓扑排序
    graph.printTopologicalSort();

    return 0;
}
