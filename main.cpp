#ifndef HGLOBAIS
#include "globais.hpp"
#endif

#include <map>
#include <string>
#include <fstream>
#include <queue>
#include <set>
#include <cstdlib>

#include "ilcplex/ilocplex.h"

using namespace std;

/* variaveis globais */

/* dados da mochila */
int v_;
int a_;              /* número de itens */
vector<int> v;   /* custos */
vector<int> origem;
vector<int> destino;
vector<int> w;      /* pesos */
vector<vector<int> > g(v_, vector<int>()); /* grafo de entrada em lista de adjacencias */
map<pair<int, int>, int> edge_to_index; /* mapeia aresta no indice correspondente */
int W;              /* capacidade */
/* contador do total de cortes por nó da arvore B&B */
int totcuts=0;
/* contador do total de lacos de separacao de cortes por nó da arvore B&B */
int itersep=0;
/* profundidade maxima de um no da arvore de B&B em que sera feita separacao */
int MAX_NODE_DEPTH_FOR_SEP=1000000;
/* - valor otimo da primeira relaxacao */
double objval_relax;
/* - valor da relaxacao linear no final do 1o nó */
double objval_node1;
/* Profundidade do nó */
IloInt node_depth;

enum ModelType {
  ILP,
  HYBRID
};

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

struct RegAux {
  double valor;
  int indice;
  friend bool operator<(const RegAux &a, const RegAux &b) {
    return (a.valor > b.valor);
  }
};

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

int main(int argc, char * argv[]) {

  /* variaveis auxiliares */
  /* - diz se usara ou nao a heuristica primal */
  ModelType model_type;
  double timelimit;
  bool use_primal_heur;

  if (argc < 5 || argc > 5) {
    cout << "Usage: " << string(argv[0]) << " <model> <time-limit> <heur-primal> <arq-input>" << endl;
    cout << "\t<model>: s for ILP, h for HYBRID" << endl;
    cout << "\t<time-limit>: Time limit in seconds" << endl;
    cout << "\t<heur-primal>: 1 for using, 0 for not using" << endl;
    cout << "\t<arq-input>: path to input instance" << endl;
    return -1;
  }
  model_type = (string(argv[1]) == "s") ? ILP : HYBRID;
  timelimit = atof(argv[2]);
  use_primal_heur = (string(argv[3]) == "1");
  ifstream fin(argv[4]);

  /* ambiente do cplex */
  IloEnv env;

  /* le dados de entrada do problema da AGMR.*/
  fin >> v_ >> a_;

  int root = 0; // No raiz, escolhido arbitrariamente
  g = vector<vector<int> >(v_, vector<int>());
  v.resize(v_);
  origem.resize(a_);
  destino.resize(a_);

  for (int i = 0; i < a_; i++) {
    fin >> origem[i] >> destino[i];
    origem[i]--;  destino[i]--; // Tornando a indexação 0-based
    g[origem[i]].push_back(destino[i]);
    g[destino[i]].push_back(origem[i]);
    edge_to_index[{origem[i], destino[i]}] = 2 * i;
    edge_to_index[{destino[i], origem[i]}] = 2 * i + 1;
  }

  /* objeto que representa o modelo */
  IloModel model(env);

  IloBoolVarArray x(env, a_);
  IloBoolVarArray y(env, v_);
  IloBoolVarArray z(env, 2 * a_);
  map<pair<int, int>, IloBoolVar*> x_map;
  map<pair<int, int>, IloBoolVar*> z_map;
  map<int, IloBoolVar*> y_map;
  for(int i = 0; i < v_; i++) y_map[i] = &y[i];
  for(int i = 0; i < a_; i++){
    x_map[{origem[i], destino[i]}] = &x[i];
    z_map[{origem[i], destino[i]}] = &z[2 * i];
    z_map[{destino[i], origem[i]}] = &z[2 * i + 1];
  }

  /* funcao objetivo (1) */
  IloExpr obj(env);
  for(int i = 0; i < v_; i++)
    obj += y[i];
  model.add(IloMinimize(env, obj)); 

  /* restricao (3) */
  IloExpr constr_3(env);
  for(int i = 0; i < a_; i++)
    constr_3 += x[i];
  model.add(constr_3 == v_ - 1);

  /* restricao (4) */
  for(int i = 0; i < v_; i++){
    IloExpr constr_4(env);
    for (int e = 0; e < g[i].size(); e++) {
      int to = g[i][e];
      constr_4 += x[edge_to_index[{i, to}] / 2];
    }
    model.add(constr_4 - 2 <= ((int)g[i].size() - 2) * y[i]);
  }

  /* restricao (27) */
  for (int i = 0; i < v_; i++) {
    if (i == root) continue;
    IloExpr constr_27(env);
    for (int e = 0; e < g[i].size(); e++) {
      int to = g[i][e];
      /* Somando arestas de entrada */
      constr_27 += z[edge_to_index[{to, i}]];
    }
    model.add(constr_27 == 1);
  }

  /* restricao (28) */
  for (int i = 0; i < v_; i++) {
    if (i == root) continue;
    IloExpr constr_28(env);
    for (int e = 0; e < g[i].size(); e++) {
      int to = g[i][e];
      /* Somando arestas de saida */
      constr_28 += z[edge_to_index[{i, to}]];
    }
    model.add(constr_28 - 1 <= ((int)g[i].size() - 2) * y[i]);
  }

  /* restricao (29) */
  IloExpr constr_29(env);
  for (int e = 0; e < g[root].size(); e++) {
    int to = g[root][e];
    constr_29 += z[edge_to_index[{root, to}]];
  }
  model.add(constr_29 - 2 <= ((int)g[root].size() - 2) * y[root]);

  /* restricao (35) */
  IloExpr constr_35(env);
  for (int i = 0; i < a_; i++) {
    IloExpr constr_35(env);
    constr_35 += z[2 * i];
    constr_35 += z[2 * i + 1];
    model.add(constr_35 == x[i]);
  }

  /* cria objeto do cplex */
  IloCplex cplex(env);

  /* carrega o modelo */
  cplex.extract(model);

  /*  atribui valores aos diferentes parametros de controle do CPLEX */
	cplex.setParam(IloCplex::Param::TimeLimit, timelimit);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve, CPX_OFF);
	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::FPHeur, -1);
	cplex.setParam(IloCplex::Param::Preprocessing::Linear, 0);
	cplex.setParam(IloCplex::Param::MIP::Strategy::Search, CPX_MIPSEARCH_TRADITIONAL);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1.0);
  
  /* salva um arquivo ".lp" com o LP original */
  cplex.exportModel("LP.lp");


  cplex.use(HeuristicaPrimal(env, x, y, z));
	cplex.use(LazyConstraints(env, x));
  cplex.use(CorteMinimo(env, x));

  cplex.solve();

  cout << "solucao otima: " << cplex.getObjValue() << endl;

  IloNumArray ystar(env);
  cplex.getValues(ystar, y);
  cout << "valores de y:" << endl;
  for(int i = 0; i < v_; i++){
    cout << "y[" << i << "]: " << ystar[i] << endl;
  }

  IloNumArray xstar(env);
  cplex.getValues(xstar, x);
  cout << "valores de x:" << endl;
  for(int i = 0; i < a_; i++){
    cout << "x[" << origem[i] << "," << destino[i] << "]: " << xstar[i] << endl;
  }

  IloNumArray zstar(env);
  cplex.getValues(zstar, z);
  cout << "valores de z:" << endl;
  for (int i = 0; i < a_; i++) {
    cout << "(" << origem[i] << ", " << destino[i] << ") " << zstar[2 * i] << endl;
    cout << "(" << destino[i] << ", " << origem[i] << ") " << zstar[2 * i + 1] << endl;
  }

  cout << "test return code: " << testOk(xstar, ystar, zstar) << endl;

  return 0;
}
