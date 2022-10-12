#include "christofides/Matching/Graph.h"
#include "christofides/TSPLIB_parser.h"
#include "christofides/Christofides.h"
#include "arestas_no_um.h"
#include <iostream>

typedef std::pair<std::pair<std::list<int>, double>, ArestasNoUm> SL;

void ordenarArestasMenores(
    int *i,
    int *j,
    double *cI,
    double *cJ
) {
    double auxC = *cI;
    *cI = *cJ;
    *cJ = auxC;
    int auxI = *i;
    *i = *j;
    *j = auxI;
}

void atualizarArestasMenores(
    int *i,
    int *j,
    int v,
    double *cI,
    double *cJ,
    double cV
) {
    if (cV < *cI) {
        *cJ = cV;
        *j = v;
        ordenarArestasMenores(i, j, cI, cJ);
    } else {
        *cJ = cV;
        *j = v;
    }
}

SL resolverSubproblemaLagrangeano(
    Graph &grafoSemUm, 
    std::vector<double> &custoSemUm, 
    const Graph &grafo, 
    std::vector<double> &custo,
    std::vector<double> &custosLagrangeanos
) {
    std::pair< list<int>, double > p = Prim(grafoSemUm, custoSemUm);

    const int zero = 0;
    std::list<int> adjUm = grafo.AdjList(zero);
    list<int>::iterator it = adjUm.begin();

    int i = *it;
    it++;
    int j = *it;
    it++;
    double cI = custo[grafo.GetEdgeIndex(zero, i)] - custosLagrangeanos[i - 1], 
        cJ = custo[grafo.GetEdgeIndex(zero, j)] - custosLagrangeanos[j - 1];
    i--;
    j--;

    if (cI > cJ) {
        ordenarArestasMenores(&i, &j, &cI, &cJ);
    }

    double cV;
    int v, vSemUm;

    while (it != adjUm.end()) {
        v = *it;
        vSemUm = v - 1;
        cV = custo[grafo.GetEdgeIndex(zero, v)] - custosLagrangeanos[vSemUm];
        
        if (cV < cI || cV < cJ) {
            atualizarArestasMenores(&i, &j, vSemUm, &cI, &cJ, cV);
        }

        it++;
    }

    ArestasNoUm no = { i, j, cI, cJ };

    return std::make_pair(p, no);
}

std::vector<double> removerVerticeUm (const Graph &grafo, std::vector<double> &custo, Graph &grafoSemUm) {
    int n = grafoSemUm.GetNumVertices();
    std::vector<double> custoSemUm(n*(n-1)/2);

    for (int i = 1; i < grafo.GetNumVertices(); i++) {
        for (int j: grafo.AdjList(i)) {
            if (j > 0) {
                grafoSemUm.AddEdge(i - 1, j - 1);
                custoSemUm[grafoSemUm.GetEdgeIndex(i - 1, j - 1)] = custo[grafo.GetEdgeIndex(i, j)];
            }
        }
    }

    return custoSemUm;
}

void relaxacaoLagrangeana (const Graph &grafo, std::vector<int> custo) {
    std::vector<double> custoD;
    custoD.assign(custo.begin(), custo.end());
    
    double Z_UB = Christofides(grafo, custoD).second;
    double Z_LB = 0;
    double pi = 2;

    Graph grafoSemUm(grafo.GetNumVertices() - 1);
    std::vector<double> custosSemUm = removerVerticeUm(grafo, custoD, grafoSemUm);
    std::vector<double> u(grafoSemUm.GetNumVertices(), 0);
    std::vector<double> custosLagrangeanos (custosSemUm.size(), 0);

    while (true) {
        SL solucao = resolverSubproblemaLagrangeano(grafoSemUm, custosSemUm, grafo, custoD, custosLagrangeanos);
        std::list mst = solucao.first.first;
        double custoMst = solucao.first.second;
        ArestasNoUm arestasNoUm = solucao.second;

        Z_LB = custoMst + arestasNoUm.cI + arestasNoUm.cJ;

        printf("%lf %d %d\n", Z_LB, arestasNoUm.i, arestasNoUm.j);

        break;
    }

    for (int j: grafo.AdjList(0)) {
        printf("%d %lf\n", j, custoD[grafo.GetEdgeIndex(0, j)]);
    }

}

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        std::cout << "Informe o nome do arquivo" << std::endl;
    }

    string arquivo(argv[1]);

    TSPLIB_parser parser(arquivo);
    Graph g = parser.GetGraph();

    relaxacaoLagrangeana(g, parser.GetCosts());

    return 0;
}
