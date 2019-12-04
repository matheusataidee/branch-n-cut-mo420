#ifndef HYBRID_H_
#define HYBRID_H_

#include <iostream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>

#include "ilcplex/ilocplex.h"
#include "utils.hpp"

using namespace std;

extern int v_;
extern int a_;
extern vector<int> origem;
extern vector<int> destino;
extern vector<vector<int> > g;
extern map<pair<int, int>, int> edge_to_index;

extern map<int, IloBoolVar*> y; /* dicionario de ponteiros para as variaveis Y */
extern map<std::string, IloBoolVar*> x; /* dicionario de ponteiros para as variaveis X */  
extern vector<int> arestas_ponte; /* armazena as pontes encontradas no pre processamento */
extern vector<int> vertices_corte; /* armazena os vertices de corte encontrados no pre processamento */

/* declaracao das necessidades do pre processamento */
extern vector<vector<int>> listaAdj;
extern int tempo_global;
extern vector<bool> visitado;
extern vector<int> tempo_entrada; /* Tempo que a árvore DFS chega em v */
extern vector<int> menor_retorno; /* Menor tempo de entrada que v alcança na árvore DFS */

// void pre_processamento();
// void dfs_cortes(int v, int pai = -1);

extern IloBoolVarArray* y_global;     
extern IloBoolVarArray* x_global;     

/* ---- conjunto de variaveis ainda nao utilizadas devidamente (estavam no exemplo) ---- */
	
	extern int totcuts; /* contador do total de cortes por nó da arvore B&B */
	extern int itersep; /* contador do total de lacos de separacao de cortes por nó da arvore B&B */
	extern int MAX_NODE_DEPTH_FOR_SEP; /* profundidade maxima de um no da arvore de B&B em que sera feita separacao */
	extern double objval_relax; /* valor otimo da primeira relaxacao */
	extern double objval_node1; /* valor da relaxacao linear no final do 1o nó */

	extern double zstar; /* melhor limitante primal encontrado */
	extern double melhor_limitante_dual; /* melhor limitante dual encontrado */
	extern int incumbent_node; /* no do melhor limitante primal encontrado */

	extern int contador_sec;
	extern int contador_d18;
	extern int contador_d19;
	extern int contador_d34;
	extern int total_no_exp;

	extern IloInt node_depth; /* profundidade do nó */

/* Corte baseado em min-cut do grafo
  O algoritmo utilizado é o Karger's Algorithm, onde
  repete-se o seguinte procedimento diversas vezes:
    - Contrair duas arestas aleatorias repetidamente
      ate obter apenas dois vertices
    - Soma o peso de todas arestas com um extremo em
      cada um dos dois vertices resultantes
    - Verifica se essa soma de pesos eh menor que um
    - Caso seja, o corte eh adicionado*/
ILOUSERCUTCALLBACK1(CorteMinimo, IloBoolVarArray, x) {
  cout << "Entrei no callback" << endl;

   /* Recupera ambiente do cplex */
  IloEnv env = getEnv();

  /* Pega a solução do LP. */
  IloNumArray val(env);
  getValues(val, x);

  vector<int> r(v_);
  vector<int> p(v_);
  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < v_; j++) makeSet(j, r, p);
    int n_componentes = v_;
    /* Contracao de arestas ate so restar dois vertices */
    while (n_componentes > 2) {
      int aresta = rand() % a_;
      n_componentes -= uniteSets(origem[aresta], destino[aresta], r, p);
    }

    double cut = 0;
    for (int j = 0; j < a_; j++) {
      /* Somar arestas ligando os dois vertices resultantes */
      if (findSet(origem[j], r, p) != findSet(destino[j], r, p))
        cut += val[j];
    }

    /* Adicionando corte */
    if (cut < 1) {
      cout << "Cut pequeno encontrado na iteraçao: " << i << " cut: " << cut << endl;
      for (int j = 0; j < v_; j++) {
        cout << j << ": " << findSet(j, r, p) << endl;
      }
      IloExpr cut(env);
      for (int j = 0; j < a_; j++) {
        if (findSet(origem[j], r, p) != findSet(destino[j], r, p)) {
          cut += x[j];
        }
      }
      add(cut >= 1);
      break;
    }
  }
}

ILOLAZYCONSTRAINTCALLBACK1(LazyConstraints, IloBoolVarArray, x) {
  cout << "Entrou no Lazy constraints" << endl;
  /* Recupera ambiente do cplex */
  IloEnv env = getEnv();

  /* Pega a solução do LP. */
  IloNumArray val(env);
  getValues(val, x);

  vector<set<int> > components;
  vector<int> used(v_, false);
  for (int i = 0; i < v_; i++) {
    if (!used[i]) {
      set<int> comp;
      queue<int> myq;
      myq.push(i);
      while (!myq.empty()) {
        int p = myq.front();  myq.pop();
        comp.insert(p);
        for (int i = 0; i < g[p].size(); i++) {
          int to = g[p][i];
          if (!used[to] && val[edge_to_index[{p, to}] / 2] == 1) {
            used[to] = true;
            myq.push(to);
          }
        }
      }
      components.push_back(comp);
    }
  }
  for (int i = 0; i < components.size(); i++) {
    IloExpr cut(env);
    for (int e = 0; e < a_; e++) {
      if (components[i].count(origem[e]) && components[i].count(destino[e])) {
        cut += x[e];
      }
    }
    double rhs = components[i].size() - 1.0;
    add(cut <= rhs);
  }
}

ILOHEURISTICCALLBACK3(HeuristicaPrimal, IloBoolVarArray, x, IloBoolVarArray, y, IloBoolVarArray,  z) {
  cout << "heuristica primal executada!" << endl;
  IloNumArray val_x(getEnv());
  IloNumArray val_y(getEnv());
  IloNumArray val_z(getEnv());
  getValues(val_x, x);
  getValues(val_y, y);
  getValues(val_z, z);

  cout << "testOk() pre heuristica: " << testOk(val_x, val_y, val_z) << endl;
  /* r e p sao vetores uteis para o algoritmo union-find */
  vector<int> r(v_);
  vector<int> p(v_);
  vector<int> deg(v_, 0); /* grau dos vertices da solucao heuristica */
  vector<bool> used(v_, false); /* util para marcar vars z durante bfs*/
  vector<RegAux> arestas_ordenadas(a_);
  for (int i = 0; i < v_; i++) makeSet(i, r, p);
  for (int i = 0; i < a_; i++) {
    arestas_ordenadas[i].valor = val_x[i];
    arestas_ordenadas[i].indice = i;
  }
  sort(arestas_ordenadas.begin(), arestas_ordenadas.end());
  /* Inicialmente cada vertice eh uma componente */
  int n_components = v_;
  for (int i = 0; i < a_; i++) {
    int id = arestas_ordenadas[i].indice;
    /* Caso em que a aresta conecta duas componentes */
    if (uniteSets(origem[id], destino[id], r, p)) {
      n_components--;
      deg[origem[id]]++;  deg[destino[id]]++;
      val_x[id] = 1.0;
    } else {
      val_x[id] = 0.0;
    }
  }
  double custo = 0;
  /* Demarca variaveis y correspondentes a vertices de ramificacoes */
  for (int i = 0; i < v_; i++) {
    if (deg[i] >= 3) {
      val_y[i] = 1.0;
      custo += 1;
    } else {
      val_y[i] = 0.0;
    }
  }

  /* BFS partindo de root (vertice 0) para marcar variaveis z corretamente
      Eh necessario pois necessita do sentido das arestas, que eh o mesmo
      sentido em que se percorre o grafo por uma bfs partindo de root.*/
  queue<int> myq;
  used[0] = true;
  myq.push(0);
  while (!myq.empty()) {
    int node = myq.front(); myq.pop();
    for (int i = 0; i < a_; i++) {
      if (origem[i] != node && destino[i] != node) continue;
      int to = origem[i] + destino[i] - node;
      if (used[to]) continue;
      used[to] = true;
      myq.push(to);
      if (origem[i] == node) {
        val_z[2 * i] = 1.0;
        val_z[2 * i + 1] = 0.0;
      } else {
        val_z[2 * i] = 0.0;
        val_z[2 * i + 1] = 1.0;
      }
    }
  }
  /* Mescla variaveis nas estruturas correspondentes */
  IloNumVarArray vars(getEnv());
  IloNumArray vals(getEnv());
  for (int i = 0; i < a_; i++)      vars.add(x[i]);
  for (int i = 0; i < v_; i++)      vars.add(y[i]);
  for (int i = 0; i < 2 * a_; i++)  vars.add(z[i]);

  for (int i = 0; i < a_; i++)      vals.add(val_x[i]);
  for (int i = 0; i < v_; i++)      vals.add(val_y[i]);
  for (int i = 0; i < 2 * a_; i++)  vals.add(val_z[i]);

  setSolution(vars, vals, custo);
  cout << "testOk() pos heuristica: " << testOk(val_x, val_y, val_z) << endl;
}

#endif  // HYBRID_H_