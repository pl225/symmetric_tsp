#ifndef FIXAR_VARIAVEIS_H
#define FIXAR_VARIAVEIS_H

#include "christofides/Matching/Graph.h"
#include "arestas_no_um.h"
#include <iostream>
#include <limits>

void fixarVariaveis (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    const std::vector<double> &custosLagrangeanos, 
    const std::vector<double> &u,
    ArestasNoUm arestas,
    ArestasNoUm *fixadasUm
);

#endif