#include "fixar_variaveis.h"

void fixarVariaveisZero (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    Graph &grafoSemUm, 
    std::vector<double> custosD,
    std::vector<double> u,
    ArestasNoUm arestas
) {
    const int zero = 0;
    list<int>::iterator it = grafo.AdjListWithoutConst(zero).begin();
    const double Z_LB_Cj = Z_LB - arestas.cJ;

    while (it != grafo.AdjList(zero).end()) {
        int v = *it;
        int vSemUm = v - 1;
        
        if (vSemUm != arestas.i && vSemUm != arestas.j) {
            double custoCandidato = grafo.GetEdgeIndex(zero, v) - u[vSemUm];
            double Z_LB_Aux = Z_LB_Cj + custoCandidato;

            if (Z_LB_Aux > Z_UB) {
                it = grafo.removeEdge(zero, v);
                continue;
            }
        }

        it++;
    }
}

void fixarVariaveis (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    Graph &grafoSemUm, 
    std::vector<double> custosD,
    std::vector<double> u,
    ArestasNoUm arestas
) { 
    fixarVariaveisZero(Z_LB, Z_UB, grafo, grafoSemUm, custosD, u, arestas);
}