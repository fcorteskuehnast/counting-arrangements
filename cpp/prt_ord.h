#pragma once

#include <vector>
#include <iostream>

#include "utils.h"

using Vertex = Line;

//#include <stddef.h>

class PartialOrder {
public:
    PartialOrder(size_t numOfV);
    void add_edge(Vertex a, Vertex b);
    std::vector<std::vector<Vertex>> all_linear_extensions();

    std::vector<Vertex> lineLabels;

    void edges();

private:
    size_t numOfV;
    std::vector<std::vector<Vertex>>	adj;
    std::vector<size_t>				    indeg;
    std::vector<std::vector<Vertex>>	results;

    void all_linear_extensions_recursion(std::vector<Vertex>& res, std::vector<bool>& visited);
};

