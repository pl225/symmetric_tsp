#include "fixar_variaveis.h"

void fixarVarsUmVerticeZero (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    const std::vector<double> &custosD,
    const std::vector<double> &u,
    ArestasNoUm arestas
) {
    const int zero = 0;
    list<int>::iterator it = grafo.AdjListWithoutConst(zero).begin();
    const double Z_LB_Cj = Z_LB - arestas.cJ;

    while (it != grafo.AdjList(zero).end()) {
        int v = *it;
        int vSemUm = v - 1;
        
        if (vSemUm != arestas.i && vSemUm != arestas.j) {
            double custoCandidato = custosD[grafo.GetEdgeIndex(zero, v)] - u[vSemUm];
            double Z_LB_Aux = Z_LB_Cj + custoCandidato;

            if (Z_LB_Aux > Z_UB) {
                it = grafo.removeEdge(zero, v);
                continue;
            }
        }

        it++;
    }
}

void fixarVarsUmVerticeUm (
    double Z_LB, 
    double Z_UB, 
    const Graph &grafo, 
    const std::vector<double> &custosD,
    const std::vector<double> &u,
    ArestasNoUm arestas,
    ArestasNoUm *fixadasUm
) {
    if (fixadasUm->i != -1 && fixadasUm->j != -1) return;

    const int zero = 0;
    std::list<int> adjZero = grafo.AdjList(zero);
    if (adjZero.size() > 2) {

        list<int>::iterator it = adjZero.begin();
        int minV = -1;
        double minCost = std::numeric_limits<double>::max();     

        while (it != adjZero.end()) {
            int v = *it;
            int vSemUm = v - 1;

            if (vSemUm != arestas.i && vSemUm != arestas.j) {
                double custoCandidato = custosD[grafo.GetEdgeIndex(zero, v)] - u[vSemUm];

                if (custoCandidato < minCost) {
                    minV = vSemUm;
                    minCost = custoCandidato;
                }
            }
            
            it++;
        }

        if (fixadasUm->j == -1) {
            if (Z_LB - arestas.cJ + minCost > Z_UB) {
                fixadasUm->j = arestas.j;
            }
        }

        if (fixadasUm->i == -1) {
            if (Z_LB - arestas.cI + minCost > Z_UB) {
                fixadasUm->i = arestas.i;
            }
        }
        
    }
}

void fixarVariaveis (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    const std::vector<double> &custosD,
    const std::vector<double> &u,
    ArestasNoUm arestas,
    ArestasNoUm *fixadasUm
) { 
    fixarVarsUmVerticeZero(Z_LB, Z_UB, grafo, custosD, u, arestas);
    fixarVarsUmVerticeUm(Z_LB, Z_UB, grafo, custosD, u, arestas, fixadasUm);
}