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
              
              if (start == end && weight < 0) {
                  cout << "Invalid graph. Exiting." << endl;
                  exit(0);
              }
          }
          
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
          for (unsigned int k = 0; k < V; ++k) {
              if (graph[k][k] < 0) return false;
              for (unsigned int i = 0; i < V; ++i) {
                  if (graph[i][k] == INF) continue;
                  
                  for (unsigned int j = 0; j < V; ++j) {
                      if (graph[k][j] == INF) continue;
                      long long newDist = (long long)graph[i][k] + graph[k][j];
                      if (newDist < INT_MIN) return false;
                      if (newDist < graph[i][j]) {
                          graph[i][j] = (int)newDist;
                          if (i == j && graph[i][i] < 0) return false;
                      }
                  }
              }
          }
          for (unsigned int i = 0; i < V; ++i) {
              if (graph[i][i] < 0) return false;
          }
          
          return true;
      }
};