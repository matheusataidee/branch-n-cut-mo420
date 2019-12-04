#ifndef HGLOBAIS
#include "globais.hpp"
#endif

#include <map>
#include <string>
#include <fstream>
#include <queue>
#include <set>
#include <cstdlib>

#include "utils.hpp"
#include "hybrid.hpp"

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
