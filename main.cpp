#include "lagrangeano_principal.h"
#include "christofides/TSPLIB_parser.h"
#include "fixar_variaveis.h"
#include <iostream>
#include <chrono>

void relaxacaoLagrangeana (Graph &grafo, std::vector<int> custo) {
    std::vector<double> custoD;
    custoD.assign(custo.begin(), custo.end());
    
    double Z_UB = Christofides(grafo, custoD).second,
        Z_LB = 0,
        Z_LB_MAX = 0,
        pi = 2,
        gQuadrado = -1;
    int iterSemMelhora = 0, iter = 0;
    ArestasNoUm fixadasUm = { -1, -1, -1, -1 };

    printf("%lf, ", Z_UB);

    Graph grafoSemUm(grafo.GetNumVertices() - 1);
    std::vector<double> custosSemUm = removerVerticeUm(grafo, custoD, grafoSemUm);
    std::vector<double> u(grafoSemUm.GetNumVertices(), 0);
    std::vector<double> custosLagrangeanos (custosSemUm);
    std::vector<double> G (grafoSemUm.GetNumVertices(), 2);

    while (deveContinuar(Z_LB, Z_UB, pi, iter, gQuadrado)) {
        SL solucao = resolverSubproblemaLagrangeano(grafoSemUm, grafo, custoD, custosLagrangeanos, u, &fixadasUm);
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
            fixarVariaveisVerticeUm(Z_LB, Z_UB, grafo, custoD, u, arestasNoUm, &fixadasUm);
            fixarVariaveisOutrosVertices(Z_LB, Z_UB, grafo, grafoSemUm, custosLagrangeanos, mst);
        }

        iter++;

        //printf("Z_LB=%lf, Z_LB_MAX=%lf, Z_UB=%lf, T=%lf, PI=%lf, iter=%d\n", Z_LB, Z_LB_MAX, Z_UB, T, pi, iter);
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
