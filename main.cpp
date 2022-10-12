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
    *cJ = cV;
    *j = v;
    if (cV < *cI) {   
        ordenarArestasMenores(i, j, cI, cJ);
    }
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

void atualizarCustosLagrangeanos(
    const Graph &grafoSemUm, 
    const std::vector<double> &custosSemUm, 
    std::vector<double> &custosLagrangeanos, 
    const std::vector<double> &u
) {
    for (int i = 0; i < grafoSemUm.GetNumVertices(); i++) {
        for (int j: grafoSemUm.AdjList(i)) {
            if (i < j) {
                int index = grafoSemUm.GetEdgeIndex(i, j);
                custosLagrangeanos[index] = custosSemUm[index] - u[i] - u[j];
            }
        }
    }
}

void atualizarMultiplicadoresLagrangeanos(
    std::vector<double> &u,
    double T,
    const std::vector<double> &G
) {
    for (int i = 0; i < u.size(); i++) {
        u[i] = max(0.0, u[i] + T * G[i]);
    }
}

double calcularSubgradiente(
    std::vector<double> &G, 
    const Graph &grafoSemUm, 
    std::list<int> &mst,
    ArestasNoUm arestasNoUm
) {
    G.assign(G.size(), 2);
    G[arestasNoUm.i]--;
    G[arestasNoUm.j]--;
    for(list<int>::iterator it = mst.begin(); it != mst.end(); it++) {
        std::pair<int, int> p = grafoSemUm.GetEdge(*it);
		int u = p.first, v = p.second;
        G[u]--;
        G[v]--;
    }

    double gQuadrado = 0;
    for (double d: G) {
        gQuadrado += (d * d);
    }
    return gQuadrado;
}

double somaMultiplicadores(const std::vector<double> &u) {
    double cU = 0;
    for (double d: u) {
        cU += d;
    }
    return cU * 2;
}

SL resolverSubproblemaLagrangeano(
    const Graph &grafoSemUm, 
    const Graph &grafo, 
    const std::vector<double> &custo,
    const std::vector<double> &custosLagrangeanos,
    const std::vector<double> &u
) {
    std::pair< list<int>, double > p = Prim(grafoSemUm, custosLagrangeanos);

    const int zero = 0;
    std::list<int> adjUm = grafo.AdjList(zero);
    list<int>::iterator it = adjUm.begin();

    int i = *it;
    it++;
    int j = *it;
    it++;
    double cI = custo[grafo.GetEdgeIndex(zero, i)] - u[i - 1], 
        cJ = custo[grafo.GetEdgeIndex(zero, j)] - u[j - 1];
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
        cV = custo[grafo.GetEdgeIndex(zero, v)] - u[vSemUm];
        
        if (cV < cI || cV < cJ) {
            atualizarArestasMenores(&i, &j, vSemUm, &cI, &cJ, cV);
        }

        it++;
    }

    ArestasNoUm no = { i, j, cI, cJ };

    return std::make_pair(p, no);
}

void relaxacaoLagrangeana (const Graph &grafo, std::vector<int> custo) {
    std::vector<double> custoD;
    custoD.assign(custo.begin(), custo.end());
    
    double Z_UB = 1610; //Christofides(grafo, custoD).second;
    double Z_LB = 0;
    double pi = 2;

    Graph grafoSemUm(grafo.GetNumVertices() - 1);
    std::vector<double> custosSemUm = removerVerticeUm(grafo, custoD, grafoSemUm);
    std::vector<double> u(grafoSemUm.GetNumVertices(), 0);
    std::vector<double> custosLagrangeanos (custosSemUm.size(), 0);
    std::vector<double> G (grafoSemUm.GetNumVertices(), 2);

    while (Z_UB - Z_LB > 1) {
        SL solucao = resolverSubproblemaLagrangeano(grafoSemUm, grafo, custoD, custosLagrangeanos, u);
        std::list mst = solucao.first.first;
        double custoMst = solucao.first.second;
        ArestasNoUm arestasNoUm = solucao.second;

        double cU = somaMultiplicadores(u);
        Z_LB = custoMst + arestasNoUm.cI + arestasNoUm.cJ + cU;

        double gQuadrado = calcularSubgradiente(G, grafoSemUm, mst, arestasNoUm);

        double T = (pi * (Z_UB - Z_LB)) / gQuadrado;

        atualizarMultiplicadoresLagrangeanos(u, T, G);

        atualizarCustosLagrangeanos(grafoSemUm, custosSemUm, custosLagrangeanos, u);

        printf("%lf\n", Z_LB);
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
