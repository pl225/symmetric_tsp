#ifndef FIXAR_VARIAVEIS_H
#define FIXAR_VARIAVEIS_H

#include "christofides/Matching/Graph.h"
#include "arestas_no_um.h"
#include <unordered_map>
#include <bits/stdc++.h>
#include <iostream>
#include <limits>

typedef std::pair<double, std::pair<int, int>> ArestaCusto;

void fixarVariaveisVerticeUm (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    const std::vector<double> &custosLagrangeanos, 
    const std::vector<double> &u,
    ArestasNoUm arestas,
    ArestasNoUm *fixadasUm
);

void fixarVariaveisOutrosVertices(
    double Z_LB, 
    double Z_UB,
    Graph &grafo,
    Graph &grafoSemUm,
    const std::vector<double> &custosL,
    std::list<int> &mst
);

#endif