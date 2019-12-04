#include "utils.hpp"
#include "globais.hpp"

extern int v_;
extern int a_;
extern vector<int> origem;
extern vector<int> destino;

// DSU (Disjoint-Set-Union, or Union-Find)  
// creates a new set consist of only the element v 
void makeSet (int v, vector<int>& r, vector<int>& p) { 
    p[v] = v; r[v] = 0; 
} 
// returns representative of the set in which v lies 
int findSet (int v, vector<int>& r, vector<int>& p) {
    return v == p[v] ? v : p[v] = findSet(p[v], r, p);
}

// unite the sets in which a and b lie
int uniteSets (int a, int b, vector<int>& r, vector<int>& p) {
    a = findSet(a, r, p);
    b = findSet(b, r, p);
    if (a != b) {
        if (r[a] < r[b]) {
            swap(a, b);
        }
        p[b] = a; 
        if (r[a] == r[b]) {
            r[a]++;
        } 
        return 1;
    } else {
      return 0;
    }
}

/* Verifica se as variaveis estao coerentes e correspondem a uma arvore
  retorna 1 quando ha inconsistencia nas variaveis z
  retorna 2 quando o numero de arestas eh diferente de |v| - 1
  retorna 3 quando ha inconsistencia nas variaveis y
  retorna 4 quando o grafo nao é conexo
  retorna 0 quando o as variaveis estao consistentes e correspondem a uma arvore */
int testOk(IloNumArray& xstar, IloNumArray& ystar, IloNumArray& zstar) {
  int edge_counter = 0;
  vector<vector<int> > g(v_, vector<int>());
  for (int i = 0; i < a_; i++) {
    if (xstar[i] > 1 - EPSILON) {
      edge_counter++;
      if (zstar[2 * i] + zstar[2 * i + 1] <= 1 - EPSILON) return 1;
      if (zstar[2 * i] + zstar[2 * i + 1] >= 1 + EPSILON) return 1;
      g[origem[i]].push_back(destino[i]);
      g[destino[i]].push_back(origem[i]);
    }
  }
  if (edge_counter != v_ - 1) return 2;

  int y_counter = 0;
  int branch_counter = 0;
  for (int i = 0; i < v_; i++) {
    if (ystar[i] > 1 - EPSILON) y_counter++;
    if (g[i].size() > 2) {
      branch_counter++;
    }
  }
  if (y_counter != branch_counter) return 3;

  vector<bool> used(v_, false);
  queue<int> myq;
  myq.push(0);
  used[0] = true;
  while (!myq.empty()) {
    int p = myq.front(); myq.pop();
    for (int i = 0; i < g[p].size(); i++) {
      int to = g[p][i];
      if (used[to]) continue;
      used[to] = true;
      myq.push(to);
    }
  }
  for (int i = 0; i < v_; i++) {
    if (!used[i]) return 4;
  }
  return 0;
}