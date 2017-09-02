#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <unistd.h>

int main(int argc, char* argv[])
{
    auto path_el = "edge_list.txt";
    auto path_al = "adjacency_list.dat";
    auto preserve_id = false;
    
    int ch;
    while ((ch = getopt(argc, argv, "e:a:p")) != -1) {
        switch (ch) {
            case 'e':
                path_el = optarg;
                break;
            case 'a':
                path_al = optarg;
                break;
            case 'p':
                preserve_id = true;
                break;
            default:
                std::cout << "unknown option: " << ch << std::endl;
                break;
        }
    }
    
    // load edge list
    std::ifstream f_el (path_el);
    if (!f_el.is_open()) {
        std::cout << "can't open file \"" << path_el << "\"" << std::endl;
        return 1;
    }

    std::unordered_map<unsigned, std::vector<unsigned>> neighbors;
    unsigned max_id = 0;
    unsigned v1, v2;
    
    std::string line;
    while (std::getline(f_el, line)) {
        if (line.find('#') != -1) continue;
        std::stringstream linestream(line);
        linestream >> v1;
        linestream >> v2;
        if (v1 == v2)
            std::cout << "warning: self edge (" << v1 << ")" << std::endl;
        max_id = std::max(max_id, std::max(v1, v2));
        neighbors[v1].push_back(v2);
        neighbors[v2].push_back(v1);
    }
    f_el.close();
    
    // mapping
    unsigned id = 0;
    std::unordered_map<unsigned, unsigned> old_to_new, new_to_old;
    for (unsigned i = 0; i < max_id + 1; i++)
        if (preserve_id || neighbors[i].size() > 0) {
            old_to_new[i] = id;
            new_to_old[id] = i;
            id++;
        }

    // save adjacency list
    std::ofstream f_al (path_al, std::ios::binary);
    if (!f_al.is_open()) {
        std::cout << "can't write file \"" << path_al << "\"" << std::endl;
        return 2;
    }
    
    unsigned network_size = new_to_old.size();
    std::cout << "number of vertices: " << network_size << std::endl;
    f_al.write((char*)&network_size, sizeof(unsigned));
    
    unsigned num_of_edges = 0;
    for (unsigned i = 0; i < network_size; i++) {
        unsigned old_id = new_to_old[i];
        
        std::set<unsigned> tmp_set (neighbors[old_id].begin(), neighbors[old_id].end());
        neighbors[old_id] = std::vector<unsigned> (tmp_set.begin(), tmp_set.end());
        unsigned degree = neighbors[old_id].size();
        
        // vertex degree
        f_al.write((char*)&degree, sizeof(unsigned));
        // neighbors
        for (auto j : neighbors[old_id]) {
            if (old_to_new.find(j) == old_to_new.end())
                std::cout << "error: vertex " << j << " (original id) is not in mapping table." << std::endl;
            else
                f_al.write((char*)&old_to_new[j], sizeof(unsigned));
        }
        num_of_edges += degree;
        if (degree == 0) {
            std::cout << "warning: vertex " << i << " has no edges." << std::endl;
        }
    }
    
    // save mapping
    for (unsigned i = 0; i < network_size; i++)
        f_al.write((char*)&new_to_old[i], sizeof(unsigned));
    f_al.close();
    
    std::cout << "number of edges (reciprocal edges counted twice): " << num_of_edges << std::endl;
}
