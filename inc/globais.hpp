#ifndef GLOBAIS_H_
#define GLOBAIS_H_

#define HGLOBAIS 1

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <numeric>

#include <cassert>

#define EPSILON 0.000001 /* uma tolerância usada em alguns cálculos com doubles*/
#define MAX_CPU_TIME 1800 /* limitando o tempo de CPU */

struct RegAux {
  double valor;
  int indice;
  friend bool operator<(const RegAux &a, const RegAux &b) {
    return (a.valor > b.valor);
  }
};

enum ModelType {
  ILP,
  HYBRID
};

#endif  // GLOBAIS_H_
