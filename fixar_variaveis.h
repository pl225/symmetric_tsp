#ifndef FIXAR_VARIAVEIS_H
#define FIXAR_VARIAVEIS_H

#include "christofides/Matching/Graph.h"
#include "arestas_no_um.h"
#include <unordered_map>
#include <set>
#include <bits/stdc++.h>
#include <queue>
#include <iostream>
#include <limits>
#include "aresta_custo.hpp"

void fixarVariaveisVerticeUm (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    std::vector<double> &custosLagrangeanos, 
    const std::vector<double> &u,
    ArestasNoUm arestas,
    ArestasNoUm *fixadasUm
);

void fixarVariaveisOutrosVertices(
    double Z_LB, 
    double Z_UB,
    Graph &grafo,
    Graph &grafoSemUm,
    std::vector<double> &custosL,
    std::list<int> &mst,
    std::unordered_map<int, std::set<int>> &mapArestasFixadas,
    std::list<std::pair<int, int>> &vecArestasFixadas,
    std::vector<double> &custosD
);

#endif