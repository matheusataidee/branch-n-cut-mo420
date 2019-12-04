#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <queue>

#include "globais.hpp"
#include "ilcplex/ilocplex.h"

using namespace std;

// DSU (Disjoint-Set-Union, or Union-Find)  
// creates a new set consist of only the element v 
void makeSet (int v, vector<int>& r, vector<int>& p); 
// returns representative of the set in which v lies 
int findSet (int v, vector<int>& r, vector<int>& p);
// unite the sets in which a and b lie
int uniteSets (int a, int b, vector<int>& r, vector<int>& p);

// Retorna 0 quando solução é coerente, e um codigo quando nao coerente
int testOk(IloNumArray& xstar, IloNumArray& ystar, IloNumArray& zstar);

void dfs(int v, vector<vector<int>> lista_adj);
bool is_connected(vector<vector<int>> lista_adj);
void imprime_solucao(ModelType model_type, int timelimit, int use_primal_heur, string input_path, IloNumArray xstar, double diff_time);
void imprime_solucao(ModelType model_type, int timelimit, int use_primal_heur, string input_path, double diff_time);
void dfs_cortes(int v, int pai);
void pre_processamento();

#endif  // UTILS_H_