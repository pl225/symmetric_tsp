#include "christofides/Matching/Graph.h"
#include "christofides/TSPLIB_parser.h"
#include "christofides/Christofides.h"
#include "arestas_no_um.h"
#include <iostream>
#include <chrono>

#define IT_MAX_PI 40
#define MENOR_PI 0.005
#define IT_MAX 2000

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

bool atualizarMelhoresValores(double Z_LB, double *Z_LB_MAX, double *pi, int *iterSemMelhora) {
    if (Z_LB > *Z_LB_MAX) {
        *Z_LB_MAX = Z_LB;
        *iterSemMelhora = 0;
        return true;
    } else {
        *iterSemMelhora++;
        if (*iterSemMelhora == IT_MAX_PI) {
            *pi /= 2;
            *iterSemMelhora = 0;
        }
    }
    return false;
}

void melhorarUbCustoComplementar(
    const Graph &grafo, 
    const std::vector<double> & custoD, 
    const Graph &grafoSemUm, 
    SL &solucao, 
    double *Z_UB
) {
    std::vector<double> custosComplementares(custoD);
    list<int> mst = solucao.first.first;
    ArestasNoUm arestaNoUm = solucao.second;
    int indexArestaGrafo;
    for (int e: mst) {
        std::pair<int, int> edge = grafoSemUm.GetEdge(e);
        indexArestaGrafo = grafo.GetEdgeIndex(edge.first + 1, edge.second + 1);
        custosComplementares[indexArestaGrafo] = 0;
    }

    indexArestaGrafo = grafo.GetEdgeIndex(0, arestaNoUm.i + 1);
    custosComplementares[indexArestaGrafo] = 0;
    indexArestaGrafo = grafo.GetEdgeIndex(0, arestaNoUm.j + 1);
    custosComplementares[indexArestaGrafo] = 0;

    std::pair<std::vector<int>, double> sol = Christofides(grafo, custosComplementares);
    double Z_UB_Novo = 0;
    for (int l: sol.first) {
        Z_UB_Novo += custoD[l];
    }
    if (Z_UB_Novo < *Z_UB) {
        *Z_UB = Z_UB_Novo;
    }
}

bool deveContinuar(double Z_LB, double Z_UB, double pi, int iter, double gQuadrado) {
	return (Z_UB - Z_LB) > 1 && iter < IT_MAX && pi >= MENOR_PI && gQuadrado != 0;
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
    
    double Z_UB = Christofides(grafo, custoD).second,
        Z_LB = 0,
        Z_LB_MAX = 0,
        pi = 2,
        gQuadrado = -1;
    int iterSemMelhora = 0, iter = 0;

    printf("%lf, ", Z_UB);

    Graph grafoSemUm(grafo.GetNumVertices() - 1);
    std::vector<double> custosSemUm = removerVerticeUm(grafo, custoD, grafoSemUm);
    std::vector<double> u(grafoSemUm.GetNumVertices(), 0);
    std::vector<double> custosLagrangeanos (custosSemUm);
    std::vector<double> G (grafoSemUm.GetNumVertices(), 2);

    while (deveContinuar(Z_LB, Z_UB, pi, iter, gQuadrado)) {
        SL solucao = resolverSubproblemaLagrangeano(grafoSemUm, grafo, custoD, custosLagrangeanos, u);
        std::list<int> mst = solucao.first.first;
        double custoMst = solucao.first.second;
        ArestasNoUm arestasNoUm = solucao.second;

        double cU = somaMultiplicadores(u);
        Z_LB = custoMst + arestasNoUm.cI + arestasNoUm.cJ + cU;

        gQuadrado = calcularSubgradiente(G, grafoSemUm, mst, arestasNoUm);

        double T = (pi * (1.01*Z_UB - Z_LB)) / gQuadrado;

        atualizarMultiplicadoresLagrangeanos(u, T, G);

        atualizarCustosLagrangeanos(grafoSemUm, custosSemUm, custosLagrangeanos, u);

        bool deveMelhorarUb = atualizarMelhoresValores(Z_LB, &Z_LB_MAX, &pi, &iterSemMelhora);
        if (deveMelhorarUb) {
            melhorarUbCustoComplementar(grafo, custoD, grafoSemUm, solucao, &Z_UB);
        }

        iter++;

        //printf("%lf\n", Z_LB);
    }
    
    printf("%lf, %lf, %lf", Z_LB_MAX, Z_LB, Z_UB);
}

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        std::cout << "Informe o nome do arquivo" << std::endl;
    }

    string arquivo(argv[1]);

    std::cout << arquivo.substr(11) << ", ";
    
    TSPLIB_parser parser(arquivo);
    Graph g = parser.GetGraph();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    relaxacaoLagrangeana(g, parser.GetCosts());
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << ", " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    return 0;
}
