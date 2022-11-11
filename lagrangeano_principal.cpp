#include "lagrangeano_principal.h"

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
        *iterSemMelhora = *iterSemMelhora + 1;
        if (*iterSemMelhora == IT_MAX_PI) {
            *pi = *pi / 2;
            *iterSemMelhora = 0;
        }
    }
    return false;
}

std::vector<double> construirCustosIniciais(
    const Graph &grafo, 
    const std::vector<double> & custoD, 
    const Graph &grafoSemUm, 
    SL &solucao
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

    return custosComplementares;
}

void melhorarUbCusto(
    const Graph &grafo, 
    const std::vector<double> & custoD, 
    const Graph &grafoSemUm, 
    SL &solucao, 
    double *Z_UB,
    TipoCusto tipoCusto,
    const std::vector<double> &custoLagrangeano
) {
    std::vector<double> custos;
    if (tipoCusto == TipoCusto::COMPLEMENTAR) {
        custos = construirCustosIniciais(grafo, custoD, grafoSemUm, solucao);
    }

    std::pair<std::vector<int>, double> sol = Christofides(grafo, custos);
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
    const std::vector<double> &u,
    ArestasNoUm *fixadasUm,
    std::list<std::pair<int, int>> &vecArestasFixadas
) {
    std::pair< list<int>, double > p = Prim(grafoSemUm, custosLagrangeanos, vecArestasFixadas);

    const int zero = 0;
    std::list<int> adjUm = grafo.AdjList(zero);
    list<int>::iterator it = adjUm.begin();

    int i = *it;
    if (fixadasUm->i != -1) {
        i = fixadasUm->i + 1;
    } else {
        it++;
    }
    
    int j = *it;
    if (fixadasUm->j != -1) {
        j = fixadasUm->j + 1;
    } else {
        it++;
    }

    double cI = custo[grafo.GetEdgeIndex(zero, i)] - u[i - 1], 
        cJ = custo[grafo.GetEdgeIndex(zero, j)] - u[j - 1];
    i--;
    j--;

    if (fixadasUm->i != -1 && fixadasUm->j != -1) {
        ArestasNoUm no = { i, j, cI, cJ };
        return std::make_pair(p, no);
    }

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