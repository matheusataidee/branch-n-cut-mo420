#ifndef HGLOBAIS
#include "globais.hpp"
#endif

#include <map>
#include <string>

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
  string input_path;

  if (argc < 5 || argc > 5) {
    cout << "Usage: " << string(argv[0]) << " <model> <time-limit> <heur-primal> <arq-input>" << endl;
    cout << "\t<model>: s for ILP, h for HYBRID" << endl;
    cout << "\t<time-limit>: Time limit in seconds" << endl;
    cout << "\t<heur-primal>: 1 for using, 0 for not using" << endl;
    cout << "\t<arq-input>: path to input instance" << endl;
    return -1;
  } else {
    model_type = (string(argv[1]) == "s") ? ILP : HYBRID;
    timelimit = atof(argv[2]);
    use_primal_heur = (string(argv[3]) == "1");
    input_path = string(argv[4]);
  }

  /* ambiente do cplex */
  IloEnv env;

  /* inicializa valores de variaveis globais */
  totcuts=0;   itersep=0;

  cout << "[Entrar com os dados da instancia:]" << endl;

  /* le dados de entrada do problema da Mochila 0-1.*/
  scanf("%d %d",&v_,&a_);
  v.resize(v_);
  origem.resize(a_);
  destino.resize(a_);

  for(int i=0;i<v_;i++) scanf("%d",&v[i]);
  for(int i=0;i<a_;i++) scanf("%d %d",&origem[i],&destino[i]);

  /* objeto que representa o modelo */
  IloModel model(env);

  map<int, IloBoolVar*> y;
  IloBoolVarArray y_temp(env, v_);
  for(int i = 0; i < v_; i++){
    y[i] = &y_temp[i];
  }

  map<std::string, IloBoolVar*> x;
  IloBoolVarArray x_temp(env, a_);
  for(int i = 0; i < a_; i++){
    x["" + to_string(origem[i]) + "-" + to_string(destino[i]) + ""] = &x_temp[i];
  }

  /* funcao objetivo (1) */
  IloExpr obj(env);
  for(int i = 0; i < v_; i++)
    obj += *y[i];
  model.add(IloMinimize(env, obj)); 

  /* restricao (3) */
  IloExpr constr_3(env);
  for(int i = 0; i < a_; i++)
    constr_3 += *x[to_string(origem[i]) + "-" + to_string(destino[i])];
  model.add(constr_3 == v_ - 1);

  /* restricao (4) */
  for(int i = 0; i < v_; i++){
    IloExpr constr_4(env);
    for(int e = 0; e < a_; e++){
      if(i == origem[e]) constr_4 += *x[to_string(i) + "-" + to_string(destino[e])];
      if(i == destino[e]) constr_4 += *x[to_string(origem[e]) + "-" + to_string(i)];
    }
    model.add(constr_4 - 2 - 100 * *y[i] <= 0);
  }

  /* cria objeto do cplex */
  IloCplex cplex(env);

  /* carrega o modelo */
  cplex.extract(model);

  /* salva um arquivo ".lp" com o LP original */
  cplex.exportModel("LP.lp");

  cplex.solve();

  cout << "solucao otima: " << cplex.getObjValue() << endl;

  IloNumArray ystar(env);
  cplex.getValues(ystar, y_temp);
  cout << "valores de y:" << endl;
  for(int i = 0; i < v_; i++){
    cout << "y[" << i << "]: " << ystar[i] << endl;
  }

  IloNumArray xstar(env);
  cplex.getValues(xstar, x_temp);
  cout << "valores de x:" << endl;
  for(int i = 0; i < a_; i++){
    cout << "x[" << origem[i] << "," << destino[i] << "]: " << xstar[i] << endl;
  }

  /* cria variaveis do modelo */
/*  IloBoolVarArray y(env, v_); */

/*  IloBoolVar x(env);  */

  return 0;
}
