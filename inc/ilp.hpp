#ifndef ILP_H_
#define ILP_H_

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
extern vector<vector<int>> ciclos_tam_3; /* ciclos de tamanho 3 */
extern vector<vector<int>> ciclos_arestas; /* ciclos de tamanho 3 */

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

ILOUSERCUTCALLBACK1(CortesILP, IloBoolVarArray, x) {

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

ILOLAZYCONSTRAINTCALLBACK1(LazyConstraintsILP, IloBoolVarArray, x) {

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

ILOHEURISTICCALLBACK2(HeuristicaPrimalILP, IloBoolVarArray, y, IloBoolVarArray, x) {

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

#endif  // ILP_H_