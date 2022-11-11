#ifndef LAGRANGEANO_PRINCIPAL_H
#define LAGRANGEANO_PRINCIPAL_H

#include "christofides/Christofides.h"
#include "arestas_no_um.h"
#include "enum_custos.hpp"

#define IT_MAX_PI 100
#define MENOR_PI 0.0005
#define IT_MAX 2000

typedef std::pair<std::pair<std::list<int>, double>, ArestasNoUm> SL;

void ordenarArestasMenores(int *i, int *j, double *cI, double *cJ);

void atualizarArestasMenores(int *i, int *j, int v, double *cI, double *cJ, double cV);

std::vector<double> removerVerticeUm (const Graph &grafo, std::vector<double> &custo, Graph &grafoSemUm);

void atualizarCustosLagrangeanos(
    const Graph &grafoSemUm, 
    const std::vector<double> &custosSemUm, 
    std::vector<double> &custosLagrangeanos, 
    const std::vector<double> &u
);

void atualizarMultiplicadoresLagrangeanos(
    std::vector<double> &u,
    double T,
    const std::vector<double> &G
);

double calcularSubgradiente(
    std::vector<double> &G, 
    const Graph &grafoSemUm, 
    std::list<int> &mst,
    ArestasNoUm arestasNoUm
);

double somaMultiplicadores(const std::vector<double> &u);

bool atualizarMelhoresValores(double Z_LB, double *Z_LB_MAX, double *pi, int *iterSemMelhora);

void melhorarUbCusto(
    const Graph &grafo, 
    const std::vector<double> & custoD, 
    const Graph &grafoSemUm, 
    SL &solucao, 
    double *Z_UB,
    TipoCusto TipoCusto,
    const std::vector<double> &custoLagrangeano
);

bool deveContinuar(double Z_LB, double Z_UB, double pi, int iter, double gQuadrado);

SL resolverSubproblemaLagrangeano(
    const Graph &grafoSemUm, 
    const Graph &grafo, 
    const std::vector<double> &custo,
    const std::vector<double> &custosLagrangeanos,
    const std::vector<double> &u,
    ArestasNoUm *fixadasUm,
    std::list<std::pair<int, int>> &vecArestasFixadas
);

#endif