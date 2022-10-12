#include "christofides/Matching/Graph.h"
#include "christofides/TSPLIB_parser.h"
#include "christofides/Christofides.h"
#include <iostream>

std::vector<double> removerVerticeUm (const Graph &grafo, std::vector<double> custo, Graph &grafoSemUm) {
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
    double pi = 2;

    Graph grafoSemUm(grafo.GetNumVertices() - 1);
    std::vector<double> custosSemUm = removerVerticeUm(grafo, custoD, grafoSemUm);
    std::vector<double> u(grafoSemUm.GetNumVertices(), 0);
    std::vector<double> custosLagrangeanos (custosSemUm.size(), 0);
    

    printf("%lf\n", Z_UB);
    printf("%lf %lf\n", custoD[grafo.GetEdgeIndex(5, 12)], custosSemUm[grafoSemUm.GetEdgeIndex(4, 11)]);


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
