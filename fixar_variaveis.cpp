#include "fixar_variaveis.h"

void fixarVarsUmVerticeZero (
    double Z_LB, 
    double Z_UB, 
    Graph &grafo, 
    std::vector<double> &custosD,
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
                grafo.changeCostEdge(zero, v, custosD);
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

        if (fixadasUm->j == -1 && arestas.j != fixadasUm->i) {
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
    std::vector<double> &custosD,
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
    std::vector<double> &custosL,
    std::list<int> &mst,
    std::unordered_map<int, std::set<int>> &mapArestasFixadas,
    std::list<std::pair<int, int>> &vecArestasFixadas,
    std::vector<double> &custosD
) {
    std::vector<int> componentes(grafoSemUm.GetNumVertices(), 0);
    std::vector<std::vector<bool>> mapaArestas (grafoSemUm.GetNumVertices(), std::vector<bool>(grafoSemUm.GetNumVertices(), false));

    std::vector<ArestaCusto> arestasOrdenadas = construirArvore(grafoSemUm, custosL, mst, mapaArestas);
    int nComponentes = 0;

    while (!arestasOrdenadas.empty()) {
        ArestaCusto e = arestasOrdenadas.back();
        const int u = e.second.first, w = e.second.second;
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
            if (i != w && componentes[i] == componentes[w] && mapaArestas[i][w]) {
                componenteW.push_back(i);
                componentes[i] = nComponentes;
            }
        }
        componentes[w] = nComponentes;

        double minCost = std::numeric_limits<double>::max();
        bool escolheuMenor = false;
        for (int i: componenteU) {
            list<int>::iterator it = grafoSemUm.AdjListWithoutConst(i).begin();

            while (it != grafoSemUm.AdjListWithoutConst(i).end()) {
                int j = *it;
                if ((i != e.second.first || j != e.second.second) && componentes[j] == componentes[w]) {
                    int index = grafoSemUm.GetEdgeIndex(i, j);
                    double novoCusto = Z_LB - e.first + custosL[index];

                    if (custosL[index] < minCost) {
                        escolheuMenor = true;
                        minCost = custosL[index];
                    }

                    if (novoCusto > Z_UB) {
                        grafoSemUm.changeCostEdge(i, j, custosL);
                        grafo.changeCostEdge(i + 1, j + 1, custosD);
                        mapaArestas[i][j] = mapaArestas[j][i] = false;
                    }
                }
                it++;
            }
        }

        if (escolheuMenor && Z_LB - e.first + minCost > Z_UB) {
            bool inserir = true;
            if (mapArestasFixadas.find(u) != mapArestasFixadas.end()) {
                if (mapArestasFixadas[u].find(w) != mapArestasFixadas[u].end()) {
                    inserir = false;        
                }
            }

            if (inserir) {
                mapArestasFixadas.insert(std::make_pair(u, std::set<int>()));
                mapArestasFixadas[u].insert(w);
                mapArestasFixadas.insert(std::make_pair(w, std::set<int>()));
                mapArestasFixadas[w].insert(u);
                vecArestasFixadas.push_back(std::make_pair(w, u));
            }
        }
    }
}