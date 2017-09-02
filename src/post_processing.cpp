#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <unistd.h>

void save_communities(const char* filename,
                      std::vector<std::vector<unsigned>>& communities)
{
    std::ofstream f (filename, std::ofstream::out);
    for (auto c: communities) {
        if (c.size() == 0) continue;
        for (unsigned i = 0; i < c.size() - 1; i++)
            f << c[i] << ' ';
        f << c[c.size()-1] << std::endl;
    }
    f.close();
}

int main(int argc, char* argv[])
{
    auto path_input = "merged_communities.dat";
    auto path_output = "detected_communities.txt";
    unsigned threshold_l = 10;
    unsigned threshold_s_min = 2;
    double threshold_f_core = 0.1;
    
    int ch;
    while ((ch = getopt(argc, argv, "i:o:l:f:s:")) != -1)
        switch (ch) {
            case 'i':
                path_input = optarg;
                break;
            case 'o':
                path_output = optarg;
                break;
            case 'l':
                threshold_l = std::atoi(optarg);
                break;
            case 'f':
                threshold_f_core = std::atof(optarg);
                break;
            case 's':
                threshold_s_min = std::atoi(optarg);
                break;
            default:
                return 1;
        }
    
    std::ifstream f(path_input, std::ifstream::binary);
    if (!f.is_open()) {
        std::cout << "Can't open " << path_input << "." << std::endl;
        return 1;
    }
    
    std::vector<std::vector<unsigned>> communities;
    while (!f.eof()) {
        unsigned id, l, size, member, score, ego;
        double g;
        std::vector<unsigned> members;
        std::vector<unsigned> egos;
        
        f.read((char*)&id, 4);
        f.read((char*)&l, 4);
        f.read((char*)&g, 8);
        f.read((char*)&size, 4);
        
        egos.reserve(l);
        for (unsigned i = 0; i < l; i++) {
            f.read((char*)&ego, 4);
            egos.push_back(ego);
        }
        
        members.reserve(size);
        for (unsigned i = 0; i < size; i++) {
            f.read((char*)&member, 4);
            f.read((char*)&score, 4);
            // exclude the member's own contribution to its score
            score -= std::count(egos.begin(), egos.end(), member);
            if (score / double(l) > threshold_f_core && score >= threshold_s_min) {
                members.push_back(member);
            }
        }
        
        if (l >= threshold_l && members.size() > 0) {
            communities.push_back(std::move(members));
        }
    }
    f.close();
    save_communities(path_output, communities);
}
