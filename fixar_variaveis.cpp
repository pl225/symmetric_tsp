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

void fixarVariaveisVerticeUm (
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

std::vector<ArestaCusto> construirArvore(
    Graph &grafoSemUm,
    const std::vector<double> &custosL,
    std::list<int> &mst,
    std::vector<std::vector<bool>> &mapaArestas
) {
    std::vector<ArestaCusto> arestasOrdenadas;
    for(list<int>::iterator it = mst.begin(); it != mst.end(); it++) {
		pair<int, int> p = grafoSemUm.GetEdge(*it);
		int u = p.first, v = p.second;
        
        arestasOrdenadas.push_back(std::make_pair(custosL[*it], std::make_pair(u, v)));

        mapaArestas[u][v] = true;
        mapaArestas[v][u] = true;
	}

    sort(arestasOrdenadas.begin(), arestasOrdenadas.end());
    return arestasOrdenadas;
}

void fixarVariaveisOutrosVertices(
    double Z_LB, 
    double Z_UB,
    Graph &grafo,
    Graph &grafoSemUm,
    const std::vector<double> &custosL,
    std::list<int> &mst,
    std::unordered_map<int, std::set<int>> &arestasFixadasUm
) {
    std::vector componentes(grafoSemUm.GetNumVertices(), 0);
    std::vector<std::vector<bool>> mapaArestas (grafoSemUm.GetNumVertices(), std::vector<bool>(grafoSemUm.GetNumVertices(), false));

    std::vector<ArestaCusto> arestasOrdenadas = construirArvore(grafoSemUm, custosL, mst, mapaArestas);
    int nComponentes = 0;

    while (!arestasOrdenadas.empty()) {
        ArestaCusto e = arestasOrdenadas.back();
        int u = e.second.first, w = e.second.second;
        arestasOrdenadas.pop_back();
        mapaArestas[w][u] = false;
        mapaArestas[u][w] = false;
        std::queue<int> fila;
        std::vector<bool> visitados (grafoSemUm.GetNumVertices(), false);
        std::vector<int> componenteU, componenteW;
        
        fila.push(u);
        visitados[u] = true;
        componenteU.push_back(u);
        componenteW.push_back(w);
        nComponentes++;

        while (!fila.empty()) {
            int v = fila.back();
            fila.pop();

            for (int a: grafoSemUm.AdjList(v)) {
                if (!visitados[a] && mapaArestas[v][a]) {
                    fila.push(a);
                    visitados[a] = true;
                    componentes[a] = nComponentes;
                    componenteU.push_back(a);
                }
            }
        }
        componentes[u] = nComponentes;

        nComponentes++;
        for (int i = 0; i < grafoSemUm.GetNumVertices(); i++) {
            if (i != w && componentes[i] == componentes[w]) {
                componenteW.push_back(i);
                componentes[i] = nComponentes;
            }
        }
        componentes[w] = nComponentes;

        double minCost = std::numeric_limits<double>::max();
        int y = -1, z = -1;
        for (int i: componenteU) {
            list<int>::iterator it = grafoSemUm.AdjListWithoutConst(i).begin();

            while (it != grafoSemUm.AdjListWithoutConst(i).end()) {
                int j = *it;
                if ((i != e.second.first || j != e.second.second) && componentes[j] == componentes[w]) {
                    int index = grafoSemUm.GetEdgeIndex(i, j);
                    double novoCusto = Z_LB - e.first + custosL[index];

                    if (custosL[index] < minCost) {
                        minCost = custosL[index];
                        y = i; z = j;
                    }

                    if (novoCusto > Z_UB) {
                        it = grafoSemUm.removeEdge(i, j);
                        mapaArestas[i][j] = mapaArestas[j][i] = false;
                        continue;
                    }
                }
                it++;
            }
        }

        if (y != - 1 && Z_LB - e.first + minCost > Z_UB) {
            arestasFixadasUm.insert(std::make_pair(y, std::set<int>()));
            arestasFixadasUm[y].insert(z);
            arestasFixadasUm.insert(std::make_pair(z, std::set<int>()));
            arestasFixadasUm[z].insert(y);
        }
    }
}