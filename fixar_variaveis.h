#ifndef FIXAR_VARIAVEIS_H
#define FIXAR_VARIAVEIS_H

#include "christofides/Matching/Graph.h"
#include "arestas_no_um.h"
#include <iostream>

void fixarVariaveis (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    Graph &grafoSemUm, 
    std::vector<double> custosLagrangeanos, 
    std::vector<double> u,
    ArestasNoUm arestas
);

#endif