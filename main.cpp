#ifndef HGLOBAIS
#include "globais.hpp"
#endif

#include <map>
#include <string>
#include <chrono>

#include "ilcplex/ilocplex.h"

using namespace std;

int v_; /* armazena o número de vertices */
int a_; /* armazena o número de arestas */
vector<int> v; /* armazena nome dos vertices */
vector<int> origem; /* armazena a origem das arestas */
vector<int> destino; /* armazena o destino das arestas */
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
vector<vector<int>> ciclos_tam_3; /* ciclos de tamanho 3 */
vector<vector<int>> ciclos_arestas; /* ciclos de tamanho 3 */

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

/* ------------------------------------------------------------------------------------- */

enum ModelType {
  ILP,
  HYBRID
};

struct RegAux {
	double valor;
	int indice;
	friend bool operator<(const RegAux& a, const RegAux& b) {
		return (a.valor > b.valor);
	}
};

void caminhos_tam_3(int v, int inicio, int tam, vector<int> &cam) {

    visitado[v] = true;

    // Último vértice do caminho
    if (tam == 0) {

        visitado[v] = false;

        // Verificar se encontrou ciclo
        bool inicio_adj = find(listaAdj[v].begin(), listaAdj[v].end(), inicio) != listaAdj[v].end();

        if (inicio_adj && cam.size() == 2) {

            cam.push_back(inicio);
            sort(cam.begin(), cam.end());

            bool repetido;

            for (vector<int> c : ciclos_tam_3) 
                if (cam == c) repetido = true;
            
            if (!repetido)
                ciclos_tam_3.push_back(cam);

        }

        cam.clear();

    } else {

        // Buscar por profundidade todo caminho de tamanho 3 possível
        for (int i = 0; i < v_; i++) {

            bool eh_adj = find(listaAdj[v].begin(), listaAdj[v].end(), i) != listaAdj[v].end();

            if (!visitado[i] && eh_adj) {
                cam.push_back(i);
                caminhos_tam_3(i, inicio, tam - 1, cam);
            }

        }

        visitado[v] = false;

    }

}

vector<int> ciclo_arestas(vector<int> v_ciclo) {

    int v1 = v_ciclo[0], v2 = v_ciclo[1], v3 = v_ciclo[2];
    vector<int> a_ciclo;

    for (int i = 0; i < 3; i++) {
        for (int e = 0; e < a_; e++) {

            int orig = origem[e], dest = destino[e];
            int next;
            
            (i == 2 ? next = 0 : next = i + 1);

            if ((v_ciclo[i] == orig && v_ciclo[next] == dest) || (v_ciclo[i] == dest && v_ciclo[next] == orig)) {
                a_ciclo.push_back(e);
            }

        }
    }

    return a_ciclo;
}

void encontrar_ciclos() {

    visitado.assign(v_, false);

    for (int i = 0; i < v_ - 2 ; i++) {
        vector<int> cam_vazio;
        caminhos_tam_3(i,i,2,cam_vazio);
        visitado[i] = true;
    }

}

void dfs(int v, vector<vector<int>> lista_adj){

	if (!visitado[v]) {
		visitado[v] = true;
		for (int adj: lista_adj[v]) {
			if (!visitado[adj])
				dfs(adj, lista_adj);
		}
	}
}

bool is_connected(vector<vector<int>> lista_adj){

	visitado.assign(v_, false);
	dfs(0, lista_adj);
	for(int i = 0; i < v_; i++){
		if(i > 0 && !visitado[i]) /* DFS não alcançou nó */
			return false;
	}
	return true;
}

void imprime_solucao(ModelType model_type, int timelimit, int use_primal_heur, string input_path, IloNumArray xstar, double diff_time){

/*	escrita do arquivo .sol */
	ofstream file_sol(input_path + ".sol");

	if(!model_type) file_sol << "s " << timelimit << " " << use_primal_heur << " " << input_path << endl;

	for(int value = 0; value < a_; value++) 
		if(xstar[value] > 0 + EPSILON) file_sol << origem[value] << " " << destino[value] << endl;

	file_sol.close();

/*	escrita do arquivo .est */
	ofstream file(input_path + ".est");

	if(!model_type) file << "s " << timelimit << " " << use_primal_heur << " " << input_path << endl;
	file << contador_sec << endl; /* int */
	file << contador_d18 << endl; /* int */
	file << contador_d19 << endl; /* int */
	file << contador_d34 << endl; /* int */

	file.precision(6);	
	file.setf(ios::fixed);
	file << objval_relax << endl; /* double */
	file << objval_node1 << endl; /* double */

	file.precision(0);	
	file.setf(ios::fixed);
	file << total_no_exp << endl; /* int */
	file << incumbent_node << endl; /* int */
	file << zstar << endl; /* int */

	file.precision(6);	
	file.setf(ios::fixed);
	file << melhor_limitante_dual << endl; /* double */
	file << diff_time << endl; /* double */

	file.close();
}

void imprime_solucao(ModelType model_type, int timelimit, int use_primal_heur, string input_path, double diff_time){

/*	escrita do arquivo .sol */
	ofstream file_sol(input_path + ".sol");

	if(!model_type) file_sol << "s " << timelimit << " " << use_primal_heur << " " << input_path << endl;

	for(int value = 0; value < (v_ - 1); value++) 
		file_sol << -1 << " " << -1 << endl;

	file_sol.close();

/*	escrita do arquivo .est */
	ofstream file(input_path + ".est");

	if(!model_type) file << "s " << timelimit << " " << use_primal_heur << " " << input_path << endl;
	file << contador_sec << endl; /* int */
	file << contador_d18 << endl; /* int */
	file << contador_d19 << endl; /* int */
	file << contador_d34 << endl; /* int */

	file.precision(6);	
	file.setf(ios::fixed);
	file << objval_relax << endl; /* double */
	file << objval_node1 << endl; /* double */

	file.precision(0);	
	file.setf(ios::fixed);
	file << total_no_exp << endl; /* int */
	file << incumbent_node << endl; /* int */
	file << zstar << endl; /* int */

	file.precision(6);	
	file.setf(ios::fixed);
	file << melhor_limitante_dual << endl; /* double */
	file << diff_time << endl; /* double */

	file.close();
}

void dfs_cortes(int v, int pai = -1){

	visitado[v] = true;
	tempo_entrada[v] = tempo_global;
	menor_retorno[v] = tempo_global;
	tempo_global++;

	for(int dest : listaAdj[v]){

		if(dest == pai) continue;

		if(visitado[dest]){
			menor_retorno[v] = min(menor_retorno[v], tempo_entrada[dest]);
		}else{
			dfs_cortes(dest, v);
			menor_retorno[v] = min(menor_retorno[v], menor_retorno[dest]);

			/* Se o vértice dest adjacente à v não alcançar outro cuja entrada na árvore DFS precede ou equivale a de v, a aresta (v, dest) é uma ponte e v um vértice de corte. */
			if(menor_retorno[dest] > tempo_entrada[v]){

//				cout << "Ponte (" << v << ", " << dest << ")" << endl;
				for(int e = 0; e < a_; e++)
					if(origem[e] == v && destino[e] == dest) arestas_ponte.push_back(e);
				if(!(std::find(vertices_corte.begin(), vertices_corte.end(), v) != vertices_corte.end())) vertices_corte.push_back(v);
//              cout << "Vértice de corte: " << v << endl;

			/* Se dest alcançar v, a aresta (v, dest) não é uma ponte mas v é um vértice de corte. */
			}else if(menor_retorno[dest] == tempo_entrada[v]){

//              cout << "Vértice de corte: " << v << endl;
				if(!(std::find(vertices_corte.begin(), vertices_corte.end(), v) != vertices_corte.end())) vertices_corte.push_back(v);

			}
		}
	}
}

void pre_processamento(){

	visitado.assign(v_, false);
    tempo_entrada.assign(v_, -1);
    menor_retorno.assign(v_, -1);

	tempo_global = 0;

/*	pontes e vértices de cortes */
	for(int i = 0; i < v_; i++){
		if(!visitado[i])
			dfs_cortes(i);
	}
}

void makeSet(int v, vector<int>& r, vector<int>& p){ 
	p[v] = v; r[v] = 0; 
}

int findSet(int v, vector<int>& r, vector<int>& p){
	return v == p[v] ? v : p[v] = findSet(p[v], r, p);
} 

int uniteSets(int a, int b, vector<int>& r, vector<int>& p){

	a = findSet(a, r, p);
	b = findSet(b, r, p);

	if (a != b) {
		if(r[a] < r[b])
			swap(a, b);

		p[b] = a; 

		if(r[a] == r[b])
			r[a]++;
 
		return 1;

	}else{

		return 0;
	}
}

ILOUSERCUTCALLBACK1(Cortes, IloBoolVarArray, x) {

/*	cout << "callback executado!" << endl; */

	vector<vector<int>> components;
	vector<int> inserted(v_, 0);
	int alocado = 0;
	int novo;
	int size_comp = 1;

	vector<int> row;

/*  recupera o ambiente do cplex */
	IloEnv env = getEnv();

/*  pega a solução vigente */
	IloNumArray val(env);
	getValues(val, x);

	IloNumArray y_val(env);
	getValues(y_val, *y_global);

	vector<int> r(v_);
	vector<int> p(v_);

	for(int i = 0; i < 1000; i++){

		for (int j = 0; j < v_; j++) makeSet(j, r, p);
		int n_componentes = v_;

/*		Contracao de arestas ate so restar dois vertices */
		while(n_componentes > 2){
			int aresta = rand() % a_;
			n_componentes -= uniteSets(origem[aresta], destino[aresta], r, p);
		}

		double cut = 0;
		for(int j = 0; j < a_; j++){
/* 			somar arestas ligando os dois vertices resultantes */
			if(findSet(origem[j], r, p) != findSet(destino[j], r, p))
				cut += val[j];
		}

/* 		adicionando corte */
		if(cut < 1 - EPSILON){
			IloExpr cut(env);
			for(int j = 0; j < a_; j++)
				if(findSet(origem[j], r, p) != findSet(destino[j], r, p))
					cut += x[j];

			add(cut >= 1);
			contador_sec++;
			break;
		}
	}

/*	cria os cortes representados pelas desigualdades validas do tipo (18) */
	for(int i = 0; i < v_; i++){
		vector<int> arestas;
		for(int e = 0; e < a_; e++)
			if(val[e] > 1.0 - EPSILON && (i == origem[e] || i == destino[e])) arestas.push_back(e);

		if(arestas.size() > 2 && y_val[i] < 1.0 - EPSILON){
			IloExpr cut(env);
			for(int value : arestas) cut += x[value];
			double lhs = double(arestas.size() - 2);
//			cout << cut << " - " << lhs << " * " << *y[i] << " <= " << 2 <<  endl; 
			add(cut - lhs * *y[i] <= 2);
			contador_d18++;
		}
	}
/*  fim do corte pelas desigualdades (18) */

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

					if(k == origem[e] && val[e] > 1 - EPSILON && inserted[destino[e]] == 0){
/*						cout << "["<< origem[e] << "," << destino[e] << "]" << endl; */
						row.push_back(destino[e]);
						inserted[destino[e]] = 1;
						alocado++;
						novo = 0;
					}

					if(k == destino[e] && val[e] > 1 - EPSILON && inserted[origem[e]] == 0){
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
				size_comp++;
			}
		}

	}while(alocado < v_);

	components.push_back(row);

/*	avaliar em qual componente o vertice do corte esta e se ele esta em posicao de separacao [preparacao para (19)] */
	if(size_comp == 2){
		int position = -1;
		for(const std::vector<int> &component : components){
			position++;
			if(component.size() > 2){
				for(int cut_vertex : vertices_corte){
					if(std::find(component.begin(), component.end(), cut_vertex) != component.end() && y_val[cut_vertex] < 1 - EPSILON){

						int index;
						if(position == 0)
							index = 1;
						else if (position == 1)
							index = 0;

						int separado = 0;
						for(int e = 0; e < a_; e++){
							if(origem[e] == cut_vertex){
								if(std::find(components[index].begin(), components[index].end(), destino[e]) != components[index].end()) separado = 1;
							}else if(destino[e] == cut_vertex){
								if(std::find(components[index].begin(), components[index].end(), origem[e]) != components[index].end()) separado = 1;
							}
						}

						if(separado){
							for(int v1 = 0; v1 < component.size() - 1; v1++){
								for(int v2 = v1 + 1; v2 < component.size(); v2++){

									int encontrado = 0;
									IloExpr cut(env);

									for(int e = 0; e < a_; e++){
										if(origem[e] == cut_vertex && destino[e] == component[v1]){ 
											cut += x[e];
											encontrado++;
										}
										if(destino[e] == cut_vertex && origem[e] == component[v1]){
											cut += x[e];
											encontrado++;
										}
										if(origem[e] == cut_vertex && destino[e] == component[v2]){
											cut += x[e];
											encontrado++;
										}
										if(destino[e] == cut_vertex && origem[e] == component[v2]){
											cut += x[e];
											encontrado++;
										}
									}
/*								cria os cortes representados pelas desigualdades validas do tipo (19) */
								if(encontrado == 2 && y_val[cut_vertex] < 1 - EPSILON){
										add(cut - *y[cut_vertex] - 1 <= 0);
										contador_d19++;
//										for(int z : component) cout << z << " ";
//										cout << endl << "vertices: " << cut_vertex << ", (" << component[v1] << ", " << component[v2] << ")" << endl;
									}
								}
							}
						}
					}
				}
			}
		}
	}

/*  cria uma restrição tipo mochila (para cada componente) para suprimir subciclos nas componentes encontradas */
	for(const std::vector<int> &component : components){

		if(component.size() < v_){

/*			for(int value : component) std::cout << value << ' ';
			std::cout << std::endl; */

			IloExpr cut(env);
			int arestas_na_solucao = 0;
			for(int e = 0; e < a_; e++){

				if(std::find(component.begin(), component.end(), origem[e]) != component.end() && std::find(component.begin(), component.end(), destino[e]) != component.end()){
					cut += x[e];
					if(val[e] > EPSILON) arestas_na_solucao++;
				}
			}

			double rhs = double(component.size() - 1.0);
			if(arestas_na_solucao > rhs){
//				cout << cut << " <= " << rhs << endl;
				add(cut <= rhs);
				contador_sec++;
			}
		}
	}

}

ILOLAZYCONSTRAINTCALLBACK1(LazyConstraints, IloBoolVarArray, x) {

/*	recupera ambiente do cplex */
	IloEnv env = getEnv();

/*	pega a solução do LP. */
	IloNumArray val(env);
	getValues(val, x);

	vector<int> r(v_);
	vector<int> p(v_);

	for(int i = 0; i < 1000; i++){

		for (int j = 0; j < v_; j++) makeSet(j, r, p);
		int n_componentes = v_;

/*		Contracao de arestas ate so restar dois vertices */
		while(n_componentes > 2){
			int aresta = rand() % a_;
			n_componentes -= uniteSets(origem[aresta], destino[aresta], r, p);
		}

		double cut = 0;
		for(int j = 0; j < a_; j++){
/* 			somar arestas ligando os dois vertices resultantes */
			if(findSet(origem[j], r, p) != findSet(destino[j], r, p))
				cut += val[j];
		}

/* 		adicionando corte */
		if(cut < 1 - EPSILON){
			IloExpr cut(env);
			for(int j = 0; j < a_; j++)
				if(findSet(origem[j], r, p) != findSet(destino[j], r, p))
					cut += x[j];

			add(cut >= 1);
			contador_sec++;
			break;
		}
	}

	IloNumArray y_val(env);
	getValues(y_val, *y_global);

	vector<vector<int>> components;
	vector<int> inserted(v_, 0);
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

					if(k == origem[e] && val[e] > 1 - EPSILON && inserted[destino[e]] == 0){
/*						cout << "["<< origem[e] << "," << destino[e] << "]" << endl; */
						row.push_back(destino[e]);
						inserted[destino[e]] = 1;
						alocado++;
						novo = 0;
					}

					if(k == destino[e] && val[e] > 1 - EPSILON && inserted[origem[e]] == 0){
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

	}while(alocado < v_);

	components.push_back(row);
/*
	for(vector<int> comp : components){
		for(int k : comp) cout << k << " ";
		cout << endl;
	}
	cout << endl; // */

/*  cria uma restrição tipo mochila (para cada componente) para suprimir subciclos nas componentes encontradas */
	for(const std::vector<int> &component : components){

		if(component.size() < v_){

/*			for(int value : component) std::cout << value << ' ';
			std::cout << std::endl; */

			IloExpr cut(env);
			int arestas_na_solucao = 0;
			for(int e = 0; e < a_; e++){

				if(std::find(component.begin(), component.end(), origem[e]) != component.end() && std::find(component.begin(), component.end(), destino[e]) != component.end()){
					cut += x[e];
					if(val[e] > EPSILON) arestas_na_solucao++;
				}
			}

			double rhs = double(component.size() - 1.0);
			if(arestas_na_solucao > rhs){
//				cout << cut << " <= " << rhs << endl;
				add(cut <= rhs);
				contador_sec++;
			}
		}
	}

}

ILOHEURISTICCALLBACK2(HeuristicaPrimal, IloBoolVarArray, y, IloBoolVarArray, x) {

	IloNumArray val_x(getEnv());
	IloNumArray val_y(getEnv());
	getValues(val_x, x);
	getValues(val_y, y);

//  cout << "testOk() pre heuristica: " << testOk(val_x, val_y, val_z) << endl;
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
		}else{
			val_x[id] = 0.0;
		}
	}
	double custo = 0;
  /* Demarca variaveis y correspondentes a vertices de ramificacoes */
	for (int i = 0; i < v_; i++) {
		if (deg[i] >= 3) {
				val_y[i] = 1.0;
				custo += 1;
			}else{
				val_y[i] = 0.0;
		}
	}

	IloNumVarArray vars(getEnv());
	IloNumArray vals(getEnv());

	for (int i = 0; i < a_; i++) vars.add(x[i]);
	for (int i = 0; i < v_; i++) vars.add(y[i]);

	for (int i = 0; i < a_; i++) vals.add(val_x[i]);
	for (int i = 0; i < v_; i++) vals.add(val_y[i]);

	setSolution(vars, vals, custo);
}

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
	encontrar_ciclos();

	for(vector<int> v_ciclos : ciclos_tam_3){
		vector<int> a_ciclos = ciclo_arestas(v_ciclos);
		ciclos_arestas.push_back(a_ciclos);
		for(int k : a_ciclos) cout << k << " ";
	}
	cout << endl; // */

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
	cplex.setOut(env.getNullStream());

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

	if(use_primal_heur == 1) cplex.use(HeuristicaPrimal(env, y_temp, x_temp));
	cplex.use(LazyConstraints(env, x_temp));
	cplex.use(Cortes(env, x_temp));

	int achou_solucao = 0;

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

	return 0;
}
