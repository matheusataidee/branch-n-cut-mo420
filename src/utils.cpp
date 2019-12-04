#include "utils.hpp"

extern int v_;
extern int a_;
extern vector<int> origem;
extern vector<int> destino;

extern vector<int> arestas_ponte; /* armazena as pontes encontradas no pre processamento */
extern vector<int> vertices_corte; /* armazena os vertices de corte encontrados no pre processamento */

/* declaracao das necessidades do pre processamento */
extern vector<vector<int>> listaAdj;
extern int tempo_global;
extern vector<bool> visitado;
extern vector<int> tempo_entrada; /* Tempo que a árvore DFS chega em v */
extern vector<int> menor_retorno; /* Menor tempo de entrada que v alcança na árvore DFS */

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
		if(xstar[value] == 1) file_sol << origem[value] << " " << destino[value] << endl;

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