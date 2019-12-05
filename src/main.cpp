#ifndef HGLOBAIS
#include "globais.hpp"
#endif

#include <map>
#include <string>
#include <fstream>
#include <queue>
#include <set>
#include <cstdlib>
#include <chrono>

#include "utils.hpp"
#include "hybrid.hpp"
#include "ilp.hpp"

#include "ilcplex/ilocplex.h"

using namespace std;

/* variaveis globais */

/* dados da mochila */
int v_;
int a_;              /* número de itens */
vector<int> v;   /* custos */
vector<int> origem;
vector<int> destino;
vector<vector<int> > g(v_, vector<int>()); /* grafo de entrada em lista de adjacencias */
map<pair<int, int>, int> edge_to_index; /* mapeia aresta no indice correspondente */
map<int, IloBoolVar*> y; /* dicionario de ponteiros para as variaveis Y */
map<std::string, IloBoolVar*> x; /* dicionario de ponteiros para as variaveis X */  
vector<int> arestas_ponte; /* armazena as pontes encontradas no pre processamento */
vector<int> vertices_corte; /* armazena os vertices de corte encontrados no pre processamento */

/* declaracao das necessidades do pre processamento */
vector<vector<int>> listaAdj;
int tempo_global;
vector<bool> visitado;
vector<int> tempo_entrada; /* Tempo que a árvore DFS chega em v */
vector<int> menor_retorno; /* Menor tempo de entrada que v alcança na árvore DFS */

// void pre_processamento();
// void dfs_cortes(int v, int pai = -1);

IloBoolVarArray* y_global;     
IloBoolVarArray* x_global;     

/* ---- conjunto de variaveis ainda nao utilizadas devidamente (estavam no exemplo) ---- */
	
	int totcuts = 0; /* contador do total de cortes por nó da arvore B&B */
	int itersep = 0; /* contador do total de lacos de separacao de cortes por nó da arvore B&B */
	int MAX_NODE_DEPTH_FOR_SEP = 1000000; /* profundidade maxima de um no da arvore de B&B em que sera feita separacao */
	double objval_relax; /* valor otimo da primeira relaxacao */
	double objval_node1; /* valor da relaxacao linear no final do 1o nó */

	double zstar; /* melhor limitante primal encontrado */
	double melhor_limitante_dual; /* melhor limitante dual encontrado */
	int incumbent_node = -1; /* no do melhor limitante primal encontrado */

	int contador_sec = 0;
	int contador_d18 = 0;
	int contador_d19 = 0;
	int contador_d34 = 0;
	int total_no_exp = 0;

	IloInt node_depth; /* profundidade do nó */

int main(int argc, char * argv[]) {
  
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

/*  variaveis auxiliares */
	ModelType model_type;
	double timelimit;
	string input_path;

/*  informa se usara ou nao a heuristica primal */
	bool use_primal_heur;

/*  define o ambiente do cplex */
	IloEnv env;

	if(argc < 5 || argc > 5){
		cout << "Usage: " << string(argv[0]) << " <model> <time-limit> <heur-primal> <arq-input>" << endl;
		cout << "\t<model>: s for ILP, h for HYBRID" << endl;
		cout << "\t<time-limit>: Time limit in seconds" << endl;
		cout << "\t<heur-primal>: 1 for using, 0 for not using" << endl;
		cout << "\t<arq-input>: path to input instance" << endl;
		return -1;
	}else{
		model_type = (string(argv[1]) == "s") ? ILP : HYBRID;
		timelimit = atof(argv[2]);
		use_primal_heur = (string(argv[3]) == "1");
		input_path = string(argv[4]);
	}

  if (model_type == ILP) {

        cout << input_path << endl;

    /*  inicializa valores de variaveis globais de corte e separacao */
      totcuts = 0;   
      itersep = 0;

    //	cout << "[Entrar com os dados da instancia: (copiar e colar todos os valores)]" << endl;

      ifstream file(input_path);

      int descarta;

    /*  le dados de entrada do problema */
    //	scanf("%d %d", &v_, &a_);
      file >> v_ >> a_ >> descarta;
      v.resize(v_);
      origem.resize(a_);
      destino.resize(a_);
      listaAdj.resize(v_);

    //	iota(v.begin(), v.end(), 0);

    //	for(int i = 0; i < v_; i++) scanf("%d", &v[i]);
      for(int i = 0; i < v_; i++) v[i] = i;

    //	for(int i = 0; i < a_; i++) scanf("%d %d", &origem[i], &destino[i]);

      int orig, dest;
      for(int i = 0; i < a_; i++){
    //		scanf("%d %d", &orig, &dest);
        file >> orig >> dest >> descarta;

        origem[i] = orig - 1;
        destino[i] = dest - 1;

        listaAdj[orig - 1].push_back(dest - 1);
        listaAdj[dest - 1].push_back(orig - 1);
      }

      file.close();

        pre_processamento();

    //	for(int value : arestas_ponte) cout << "aresta ponte: " << value + 14 << " - (" << origem[value] << ", " << destino[value] << ")" << endl;
    //	for(int value : vertices_corte) cout << "vertices corte: " << value << endl;

      if(model_type) return -1;

    /*  objeto que representa o modelo */
      IloModel model(env);

      IloBoolVarArray y_temp(env, v_);
      y_global = &y_temp;
      for(int i = 0; i < v_; i++){
        y[i] = &y_temp[i];
      }

      IloBoolVarArray x_temp(env, a_);
      x_global = &x_temp;
      for(int i = 0; i < a_; i++){
        x["" + to_string(origem[i]) + "-" + to_string(destino[i]) + ""] = &x_temp[i];
      }

    /*  define a funcao objetivo (1) */
      IloExpr obj(env);
      for(int i = 0; i < v_; i++)
        obj += *y[i];
      model.add(IloMinimize(env, obj)); 

    /*  cria restricao (3) */
      IloExpr constr_3(env);
      for(int i = 0; i < a_; i++)
        constr_3 += *x[to_string(origem[i]) + "-" + to_string(destino[i])];
      model.add(constr_3 == v_ - 1);

    /*  cria familia de restricoes (4) */
      for(int i = 0; i < v_; i++){
        IloExpr constr_4(env);
        for(int e = 0; e < a_; e++){
          if(i == origem[e]) constr_4 += *x[to_string(i) + "-" + to_string(destino[e])];
          if(i == destino[e]) constr_4 += *x[to_string(origem[e]) + "-" + to_string(i)];
        }
        model.add(constr_4 - 2 - 7 * *y[i] <= 0); /* ajustar o valor 100 para d(v) */
      } 

    //	for(int value : arestas_ponte) model.add(x_temp[value] == 1);

    /*  cria objeto do cplex */
      IloCplex cplex(env);

    /*  carrega o modelo */
      cplex.extract(model);

    /*  silencia o cplex no terminal */
    //	cplex.setOut(env.getNullStream());

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

    /*  salva um arquivo ".lp" com o LP original */
      cplex.exportModel("LP.lp");

      if(use_primal_heur == 1) cplex.use(HeuristicaPrimalILP(env, y_temp, x_temp));
      cplex.use(LazyConstraintsILP(env, x_temp));
      cplex.use(CortesILP(env, x_temp));

      bool achou_solucao;

      if(cplex.solve()){

        incumbent_node = cplex.getIncumbentNode();
        zstar = cplex.getObjValue();
        melhor_limitante_dual = cplex.getBestObjValue();
        total_no_exp = cplex.getNnodes();

      //	cout << "solucao otima: " << zstar << endl << endl;

        IloNumArray ystar(env);
        cplex.getValues(ystar, y_temp);
      /*	cout << "valores de y:" << endl;
        for(int i = 0; i < v_; i++){
          if(ystar[i] > 0) cout << "\ty[" << i << "]: " << ystar[i] << endl;
        } */

        IloNumArray xstar(env);
        cplex.getValues(xstar, x_temp);
      /*	cout << endl << "valores de x:" << endl;
        for(int i = 0; i < a_; i++){
          if(xstar[i] > 0) cout << "\tx[" << origem[i] << "," << destino[i] << "]: " << xstar[i] << endl;
        } */

        chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_time = chrono::duration_cast<chrono::duration<double>>(now - start);

        imprime_solucao(model_type, timelimit, use_primal_heur, input_path, xstar, diff_time.count());

      }else{
        printf("Main: programa terminou sem achar solucao inteira !\n");

        incumbent_node = -1;
        zstar = -1;
        melhor_limitante_dual = cplex.getBestObjValue();
        total_no_exp = cplex.getNnodes();

        chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_time = chrono::duration_cast<chrono::duration<double>>(now - start);

        
        imprime_solucao(model_type, timelimit, use_primal_heur, input_path, diff_time.count());
      }

  } else if (model_type == HYBRID) {

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
    listaAdj.resize(v_);

    for (int i = 0; i < a_; i++) {
      fin >> origem[i] >> destino[i];
      origem[i]--;  destino[i]--; // Tornando a indexação 0-based
      g[origem[i]].push_back(destino[i]);
      g[destino[i]].push_back(origem[i]);
      listaAdj[origem[i]].push_back(destino[i]);
      listaAdj[destino[i]].push_back(origem[i]);
      edge_to_index[{origem[i], destino[i]}] = 2 * i;
      edge_to_index[{destino[i], origem[i]}] = 2 * i + 1;
    }

    pre_processamento();

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
    for (int i = 0; i < a_; i++) {
      IloExpr constr_35(env);
      constr_35 += z[2 * i];
      constr_35 += z[2 * i + 1];
      model.add(constr_35 == x[i]);
    }

    /* restricao (7) */
    for (int i = 0; i < v_; i++) {
      IloExpr constr_7(env);
      for (int j = 0; j < g[i].size(); j++) {
        int to = g[i][j];
        constr_7 += x[edge_to_index[{i, to}] / 2];
      }
      model.add(constr_7 - 1 >= 2 * y[i]);
    }

    /* restricao de pontes */
    for (int ponte : arestas_ponte) {
      IloExpr constr_ponte(env);
      constr_ponte += x[ponte];
      model.add(constr_ponte == 1);
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


    if(use_primal_heur) cplex.use(HeuristicaPrimalHybrid(env, x, y, z));
    cplex.use(LazyConstraintsHybrid(env, x));
    cplex.use(CortesHybrid(env, x, y, z));

    if(cplex.solve()){

        incumbent_node = cplex.getIncumbentNode();
        zstar = cplex.getObjValue();
        melhor_limitante_dual = cplex.getBestObjValue();
        total_no_exp = cplex.getNnodes();

        IloNumArray ystar(env);
        IloNumArray xstar(env);
        IloNumArray zstar_resp(env);
        cplex.getValues(ystar, y);
        cplex.getValues(xstar, x);
        cplex.getValues(zstar_resp, z);

        chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_time = chrono::duration_cast<chrono::duration<double>>(now - start);

        imprime_solucao(model_type, timelimit, use_primal_heur, input_path, xstar, diff_time.count());

        cout << "Limitante dual: " << melhor_limitante_dual << endl;
        cout << "Primal: " << zstar << endl;
        cout << "test return code: " << testOk(xstar, ystar, zstar_resp) << endl;

    } else {
        printf("Main: programa terminou sem achar solucao inteira !\n");

        incumbent_node = -1;
        zstar = -1;
        melhor_limitante_dual = cplex.getBestObjValue();
        total_no_exp = cplex.getNnodes();

        chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_time = chrono::duration_cast<chrono::duration<double>>(now - start);

        
        imprime_solucao(model_type, timelimit, use_primal_heur, input_path, diff_time.count());
    }

  }

  return 0;
}
