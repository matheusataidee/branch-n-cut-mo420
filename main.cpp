#ifndef HGLOBAIS
#include "globais.hpp"
#endif

#include <map>
#include <string>

#include "ilcplex/ilocplex.h"

using namespace std;

int v_; /* armazena o número de vertices */
int a_; /* armazena o número de arestas */
vector<int> v; /* armazena nome dos vertices */
vector<int> origem; /* armazena a origem das arestas */
vector<int> destino; /* armazena o destino das arestas */
map<int, IloBoolVar*> y; /* dicionario de ponteiros para as variaveis Y */
map<std::string, IloBoolVar*> x; /* dicionario de ponteiros para as variaveis X */  

IloBoolVarArray* y_global;     
IloBoolVarArray* x_global;     

/* ---- conjunto de variaveis ainda nao utilizadas devidamente (estavam no exemplo) ---- */

	/* contador do total de cortes por nó da arvore B&B */
	int totcuts = 0;
	/* contador do total de lacos de separacao de cortes por nó da arvore B&B */
	int itersep = 0;
	/* profundidade maxima de um no da arvore de B&B em que sera feita separacao */
	int MAX_NODE_DEPTH_FOR_SEP = 1000000;
	/* - valor otimo da primeira relaxacao */
	double objval_relax;
	/* - valor da relaxacao linear no final do 1o nó */
	double objval_node1;
	/* Profundidade do nó */
	IloInt node_depth;

/* ------------------------------------------------------------------------------------- */

enum ModelType {
  ILP,
  HYBRID
};

ILOUSERCUTCALLBACK1(Cortes, IloBoolVarArray, x) {

/*	cout << "callback executado!" << endl; */

/*  recupera o ambiente do cplex */
	IloEnv env = getEnv();

/*  pega a solução vigente */
	IloNumArray val(env);
	getValues(val, x);

	IloNumArray y_val(env);
	getValues(y_val, *y_global);

/*	cria os cortes representados pelas desigualdades validas do tipo (18) */
	for(int i = 0; i < v_; i++){
		vector<int> arestas;
		for(int e = 0; e < a_; e++)
			if(val[e] == 1 && (i == origem[e] || i == destino[e])) arestas.push_back(e);

		if(arestas.size() > 2 && y_val[i] < 1.0){
			IloExpr cut(env);
			for(int value : arestas) cut += x[value];
			double lhs = double(arestas.size() - 2);
			cout << cut << " - " << lhs << " * " << *y[i] << " <= " << 2 <<  endl; 
			add(cut - lhs * *y[i] <= 2);
		}
	}

}

ILOLAZYCONSTRAINTCALLBACK1(LazyConstraints, IloBoolVarArray, x) {

//	cout << "lazy constraint executado!" << endl;

/*  recupera o ambiente do cplex */
	IloEnv env = getEnv();

/*  pega a solução vigente */
	IloNumArray val(env);
	getValues(val, x);

	IloNumArray y_val(env);
	getValues(y_val, *y_global);

	vector<vector<int>> components;
	vector<int> inserted(v_, 0);
	vector<int> inserted_(a_, 0);
	int alocado = 0;
	int novo;

	vector<int> row;

/*  avalia se o grafo retornado tem mais de uma componente e, se sim, separa essas componentes */
	do{

		if(row.size() == 0){

/*			cout << "component novo criado" << endl; */

			for(int k = 0; k < v_; k++){

				if(inserted[k] == 0){
					row.push_back(k);
					inserted[k] = 1;
					alocado++;
					break;
				}
			}
		}else{

			novo = 1;

			for(int k : row){
				for(int e = 0; e < a_; e++){

					if(k == origem[e] && val[e] == 1 && inserted[destino[e]] == 0){
/*						cout << "["<< origem[e] << "," << destino[e] << "]" << endl; */
						row.push_back(destino[e]);
						inserted[destino[e]] = 1;
						alocado++;
						novo = 0;
					}

					if(k == destino[e] && val[e] == 1 && inserted[origem[e]] == 0){
/*						cout << "[" << origem[e] << "," << destino[e] << "]" << endl; */
						row.push_back(origem[e]);
						inserted[origem[e]] = 1;
						alocado++;
						novo = 0;
					}
				}
			}

			if(novo){
				components.push_back(row);
				row.clear();
			}
		}

	}while(alocado < 13);

	components.push_back(row);

/*  cria uma restrição tipo mochila (para cada componente) para suprimir subcilclos nas componentes encontradas */
	for(const std::vector<int> &component : components){

		if(component.size() < 13){

/*			for(int value : component) std::cout << value << ' ';
			std::cout << std::endl; */

			IloExpr cut(env);
			int arestas_na_solucao = 0;
			for(int e = 0; e < a_; e++){

				if(std::find(component.begin(), component.end(), origem[e]) != component.end() && std::find(component.begin(), component.end(), destino[e]) != component.end()){
					cut += x[e];
					if(val[e] == 1) arestas_na_solucao++;
				}

			}

			double rhs = double(component.size() - 1.0);
			if(arestas_na_solucao > rhs){
				cout << cut << " <= " << rhs << endl;
				add(cut <= rhs);
			}
		}
	}
}

ILOHEURISTICCALLBACK1(HeuristicaPrimal, IloBoolVarArray, x) {
	cout << "heuristica primal executada!" << endl;
}

int main(int argc, char * argv[]) {

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

/*  inicializa valores de variaveis globais de corte e separacao */
	totcuts = 0;   
	itersep = 0;

	cout << "[Entrar com os dados da instancia: (copiar e colar todos os valores)]" << endl;

/*  le dados de entrada do problema */
	scanf("%d %d", &v_, &a_);
	v.resize(v_);
	origem.resize(a_);
	destino.resize(a_);

	for(int i = 0; i < v_; i++) scanf("%d", &v[i]);
	for(int i = 0; i < a_; i++) scanf("%d %d", &origem[i], &destino[i]);

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

/*  cria objeto do cplex */
	IloCplex cplex(env);

/*  carrega o modelo */
	cplex.extract(model);

/*  atribui valores aos diferentes parametros de controle do CPLEX */
	cplex.setParam(IloCplex::Param::TimeLimit, MAX_CPU_TIME);
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

	cplex.use(HeuristicaPrimal(env, x_temp));
	cplex.use(LazyConstraints(env, x_temp));
	cplex.use(Cortes(env, x_temp));

	cplex.solve();

	cout << endl << "solucao otima: " << cplex.getObjValue() << endl << endl;

	IloNumArray ystar(env);
	cplex.getValues(ystar, y_temp);
	cout << "valores de y:" << endl;
	for(int i = 0; i < v_; i++){
		cout << "\ty[" << i << "]: " << ystar[i] << endl;
	}

	IloNumArray xstar(env);
	cplex.getValues(xstar, x_temp);
	cout << endl << "valores de x:" << endl;
	for(int i = 0; i < a_; i++){
		cout << "\tx[" << origem[i] << "," << destino[i] << "]: " << xstar[i] << endl;
	}

	return 0;
}
