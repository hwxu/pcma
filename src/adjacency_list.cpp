#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>

#include "adjacency_list.hpp"

// Constructor
AdjacencyList::AdjacencyList(unsigned n):
size(n), num_edges(0), degree(new unsigned[n]), neighbors(new std::unique_ptr<unsigned[]>[n]{})
{ }

// Move Constructor
AdjacencyList::AdjacencyList(AdjacencyList&& other) = default;

// Destructor
AdjacencyList::~AdjacencyList() { }

std::unique_ptr<AdjacencyList>
AdjacencyList::LoadFromBinaryFile(const char* filename, bool load_id_mapping)
{
    std::ifstream f (filename, std::ios::binary);
    if (f.is_open())
    {
        unsigned n, d, tmp;
        // read the network size {n}
        f.read((char*)&n, sizeof(unsigned));
        auto al = std::unique_ptr<AdjacencyList> (new AdjacencyList(n));
        
        // load the adjacency list
        for (unsigned i = 0; i < n; ++i)
        {
            // read the degree and the following {d} neighbors' id for Vertex i
            f.read((char*)&d, sizeof(unsigned));
            std::unique_ptr<unsigned[]> list (new unsigned[d]);
            for (unsigned j = 0; j < d; ++j)
            {
                f.read((char*)&tmp, sizeof(unsigned));
                list[j] = tmp;
                if (i <= tmp) al->num_edges++;
            }
            al->degree[i] = d;
            al->neighbors[i] = std::move(list);
        }
        // load id mapping if applicable
        if (load_id_mapping && !f.eof()) {
            for (unsigned i = 0; i < n; ++i) {
                f.read((char*)&tmp, sizeof(unsigned));
                al->mapping[i] = tmp;
            }
        }
        f.close();
        return al;
    }
    else {
        std::cout << "Can't open file \"" << filename << "\"." << std::endl;
        return std::unique_ptr<AdjacencyList> (nullptr);
    }
}

void AdjacencyList::SaveToBinaryFile(const char* filename) const
{
    std::ofstream f (filename, std::ofstream::binary);
    
    // size of the network
    f.write((char*)&size, sizeof(unsigned));
    for (unsigned i = 0; i < size; i++) {
        // degree of Vertex i
        f.write((char*)&degree[i], sizeof(unsigned));
        // neighbors of Vertex i
        for (unsigned j = 0; j < degree[i]; j++)
            f.write((char*)&neighbors[i][j], sizeof(unsigned));
    }
    // mapping (optional)
    if (mapping.size() == size) {
        for (unsigned i = 0; i < size; i++)
            f.write((char*)&mapping.at(i), sizeof(unsigned));
    }
    f.close();
}

std::unique_ptr<AdjacencyList>
AdjacencyList::CreateSubAdjacencyList(std::vector<unsigned>& vertices) const
{
    unsigned i, j;
    std::unordered_map<unsigned, unsigned> mapping_old_new;
    
    std::unique_ptr<AdjacencyList> sub (new AdjacencyList(vertices.size()));
    
    for (i = 0; i < vertices.size(); i++) {
        sub->mapping[i] = vertices[i];
        mapping_old_new[vertices[i]] = i;
    }

    std::vector<unsigned> tmp;
    for (i = 0; i < vertices.size(); i++) {
        unsigned degree_i = degree[vertices[i]];
        auto& neighbors_i = neighbors[vertices[i]];
        
        for (j = 0; j < degree_i; j++) {
            auto iter = mapping_old_new.find(neighbors_i[j]);
            if (iter != mapping_old_new.end()) {
                tmp.push_back(iter->second);
                if (i <= iter->second)
                    sub->num_edges++;
            }
        }
        sub->degree[i] = tmp.size();
        sub->neighbors[i] =
            std::unique_ptr<unsigned[]> (new unsigned[sub->degree[i]]);
        std::copy(tmp.begin(), tmp.end(), sub->neighbors[i].get());
        tmp.clear();
    }

    return sub;
}

std::unique_ptr<AdjacencyList>
AdjacencyList::FirstNeighborhoodOfVertex(unsigned id, bool exclude_ego) const
{
    std::vector<unsigned> tmp;
    tmp.reserve(degree[id]+1);
    tmp.insert(tmp.end(), neighbors[id].get(), neighbors[id].get() + degree[id]);
    if (!exclude_ego) tmp.push_back(id);
    return CreateSubAdjacencyList(tmp);
}
