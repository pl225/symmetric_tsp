#include "christofides/Matching/Graph.h"
#include "christofides/TSPLIB_parser.h"
#include "christofides/Christofides.h"
#include <iostream>

void relaxacaoLagrangeana (const Graph &grafo, std::vector<int> custo) {
    std::vector<double> custoD;
    custoD.assign(custo.begin(), custo.end());
    
    double Z_UB = Christofides(grafo, custoD).second;

    printf("%lf\n", Z_UB);
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
