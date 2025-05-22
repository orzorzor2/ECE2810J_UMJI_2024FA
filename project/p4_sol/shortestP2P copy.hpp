#include<iostream>
#include<list>
#include<vector>
#include<climits>

#define INF INT_MAX

using namespace std;

class ShortestP2P {
  public:
      ShortestP2P() {}

      void readGraph() {
          unsigned int numVertices;
          cin >> numVertices;
          V = numVertices;
          
          // 优化1：一次性分配内存，避免多次重新分配
          graph.resize(V, vector<int>(V, INF));
          
          for (unsigned int i = 0; i < V; ++i) {
              graph[i][i] = 0;
          }
          
          unsigned int numEdges;
          cin >> numEdges;
          
          for (unsigned int i = 0; i < numEdges; ++i) {
              unsigned int start, end;
              int weight;
              cin >> start >> end >> weight;
              graph[start][end] = weight;
              
              // 优化2：提前检测自环是否为负
              if (start == end && weight < 0) {
                  cout << "Invalid graph. Exiting." << endl;
                  exit(0);
              }
          }
          
          // 优化3：合并Floyd-Warshall和负环检测
          if (!computeFloydWarshall()) {
              cout << "Invalid graph. Exiting." << endl;
              exit(0);
          }
      }

      void distance(unsigned int A, unsigned int B) {
          if (A == B) {
              cout << 0 << endl;
              return;
          }
          
          if (graph[A][B] == INF) {
              cout << "INF" << endl;
              return;
          }
          
          cout << graph[A][B] << endl;
      }

  private:
      vector<vector<int>> graph;
      unsigned int V;

      bool computeFloydWarshall() {
          // 优化4：在一次遍历中完成最短路径计算和负环检测
          for (unsigned int k = 0; k < V; ++k) {
              // 优化5：立即检测经过新顶点k后是否产生负环
              if (graph[k][k] < 0) return false;
              
              for (unsigned int i = 0; i < V; ++i) {
                  if (graph[i][k] == INF) continue;
                  
                  for (unsigned int j = 0; j < V; ++j) {
                      if (graph[k][j] == INF) continue;
                      
                      // 优化6：避免整数溢出
                      long long newDist = (long long)graph[i][k] + graph[k][j];
                      if (newDist < INT_MIN) return false;
                      if (newDist < graph[i][j]) {
                          graph[i][j] = (int)newDist;
                          // 优化7：如果更新了一个顶点到自身的距离为负值，说明有负环
                          if (i == j && graph[i][i] < 0) return false;
                      }
                  }
              }
          }
          
          // 优化8：最终检查，确保没有遗漏的负环
          for (unsigned int i = 0; i < V; ++i) {
              if (graph[i][i] < 0) return false;
          }
          
          return true;
      }
};