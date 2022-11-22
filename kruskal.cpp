// C++ program for the above approach

#include <bits/stdc++.h>
#include "kruskal.hpp"

// DSU data structure
// path compression + rank by union

class DSU
{
    int *parent;
    int *rank;

public:
    DSU(int n)
    {
        parent = new int[n];
        rank = new int[n];

        for (int i = 0; i < n; i++)
        {
            parent[i] = -1;
            rank[i] = 1;
        }
    }

    ~DSU()
    {
        delete parent;
        delete rank;
    }

    // Find function
    int find(int i)
    {
        if (parent[i] == -1)
            return i;

        return parent[i] = find(parent[i]);
    }

    // Union function
    void unite(int x, int y)
    {
        int s1 = find(x);
        int s2 = find(y);

        if (s1 != s2)
        {
            if (rank[s1] < rank[s2])
            {
                parent[s1] = s2;
                rank[s2] += rank[s1];
            }
            else
            {
                parent[s2] = s1;
                rank[s1] += rank[s2];
            }
        }
    }
};

pair< list<int>, double > kruskal_mst(const Graph &G, const std::vector<double> &cost, std::list<std::pair<int, int>> &arestasFixadas)
{

    // Initialize the DSU
	DSU s(G.GetNumVertices());
    double obj = 0;
    std::vector<ArestaCusto> edgelist;
    std::list<int> edges;

    for (const std::pair<int, int> &ac : arestasFixadas)
    {
        int x = ac.first;
        int y = ac.second;
        double w = cost[G.GetEdgeIndex(x, y)];
        edges.push_back(G.GetEdgeIndex(x, y));
        s.unite(x, y);

        obj += w;
        //cout << x << " -- " << y << " == " << w << endl;
    }

    for (int i = 0; i < G.GetNumVertices(); i++) {
        for (int j: G.AdjList(i)) {
            if (i < j) {
                edgelist.push_back(std::make_pair(cost[G.GetEdgeIndex(i, j)], std::make_pair(i, j) ));
            }
        }
    }
    // 1. Sort all edges
    sort(edgelist.begin(), edgelist.end());

    /*cout << "Following are the edges in the "
            "constructed MST"
         << endl;*/
    for (auto &edge : edgelist)
    {
        double w = edge.first;
        int x = edge.second.first;
        int y = edge.second.second;

        // Take this edge in MST if it does
        // not forms a cycle
        if (s.find(x) != s.find(y))
        {
            s.unite(x, y);
            edges.push_back(G.GetEdgeIndex(x, y));
            obj += w;
            /*cout << x << " -- " << y << " == " << w
                 << endl;*/
        }
    }

    //cout << "Minimum Cost Spanning Tree Kruskal: " << obj << ' ' << edges.size() << endl;

    return pair< list<int>, double >(edges, obj);
}
