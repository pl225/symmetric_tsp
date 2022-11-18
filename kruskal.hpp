#ifndef KRUSKAL_H
#define KRUSKAL_H

#include "christofides/Matching/Graph.h"
#include "aresta_custo.hpp"
#include <iostream>

pair< list<int>, double > kruskal_mst(const Graph &G, const std::vector<double> &cost, std::list<std::pair<int, int>> &arestasFixadas);

#endif