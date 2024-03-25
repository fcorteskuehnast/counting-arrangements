#include "prt_ord.h"

PartialOrder::PartialOrder(size_t numOfV) : 
    numOfV(numOfV),
    adj(std::vector<std::vector<Vertex>>(numOfV, std::vector<Vertex>())),
    indeg(std::vector<size_t>(numOfV, 0))
{}

void PartialOrder::add_edge(Vertex a, Vertex b)
{
    adj[a].push_back(b);
    indeg[b]++;
}

void PartialOrder::edges(){
    for (Vertex i = 0; i < numOfV; i++) {
        for (Vertex j : adj[i]) {
            std::cout << "(" << i << "," << j << ")";
        }
    }
    std::cout << std::endl;
}

std::vector<std::vector<Vertex>> PartialOrder::all_linear_extensions()
{
    std::vector<bool> visited(numOfV, false);
    std::vector<Vertex> res;

    all_linear_extensions_recursion(res, visited);

    return results;
}

void PartialOrder::all_linear_extensions_recursion(std::vector<Vertex>& res, std::vector<bool>& visited)
{
    bool newVertexViseted = false;

    for (Vertex i = 0; i < numOfV; i++) {
        if (indeg[i] == 0 && !visited[i]) {
            for (Vertex neighbor : adj[i]) {
                indeg[neighbor]--;
            }

            res.push_back(i);
            visited[i] = true;

            all_linear_extensions_recursion(res, visited);

            visited[i] = false;
            res.pop_back();
            for (Vertex neighbor : adj[i]) {
                indeg[neighbor]++;
            }

            newVertexViseted = true;
        }
    }

    if (!newVertexViseted) {
        results.push_back(res);

    }
}
