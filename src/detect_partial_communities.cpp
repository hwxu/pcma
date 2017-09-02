#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <queue>
#include <sys/time.h>
#include <unistd.h>

#include "adjacency_list.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif


double count_overlap(std::vector<unsigned>& a, std::vector<unsigned>& b)
{
    unsigned count = 0, i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) i++;
        else if (a[i] > b[j]) j++;
        else {
            count++;
            i++;
            j++;
        }
    }
    return (double)count / std::min(a.size(), b.size());
}

std::vector<std::vector<unsigned>>
merge_similar_partial_communities(std::vector<std::vector<unsigned>> communities, double threshold)
{
    std::priority_queue<std::pair<double, std::pair<unsigned, unsigned>>> queue;
    for (unsigned z = 0; z < communities.size(); z++) {
        for (unsigned z2 = z + 1; z2 < communities.size(); z2++) {
            auto overlap = count_overlap(communities[z], communities[z2]);
            if (overlap > threshold)
                queue.push({overlap, {z, z2}});
        }
    }
    
    while (queue.size()) {
        auto pair = queue.top().second;
        queue.pop();
        unsigned i = pair.first;
        unsigned j = pair.second;
        if (communities[i].size() == 0 || communities[j].size() == 0)
            continue;
        
        // merge
        std::set<unsigned> tmp (communities[i].begin(), communities[i].end());
        tmp.insert(communities[j].begin(), communities[j].end());
        std::vector<unsigned> tmp2 (tmp.begin(), tmp.end());
        communities.push_back(std::move(tmp2)); // append
        communities[i].clear();
        communities[j].clear();
        
        // update queue
        for (unsigned z = 0; z < communities.size() - 1; z++) {
            if (communities[z].size() == 0) continue;
            auto overlap = count_overlap(communities.back(), communities[z]);
            if (overlap > threshold)
                queue.push({overlap, {communities.size() - 1, z}});
        }
    }
    
    std::vector<std::vector<unsigned>> merged;
    for (auto& c : communities)
        if (c.size())
            merged.push_back(std::move(c));

    return merged;
}

// An implementation of the method proposed by Brian Ball et al.
// doi:10.1103/PhysRevE.84.036103
std::unique_ptr<double[]>
find_communities(std::unique_ptr<AdjacencyList>& al,
                unsigned num_communities, unsigned loop)
{
    const double RAND_F = 10.0 / RAND_MAX;
    const double EPSILON = 0.001;
    
    unsigned i, j, z, jj, degree_i, total_size, total_bytes, kappa_bytes, iter;
    double log_p, d;
    double max_log_p = -DBL_MAX;
    double* kappa, *q_ij;
    double** k, **k_new, **k_swap, **k_best;
    
    kappa = new double[num_communities];
    kappa_bytes = sizeof(double) * num_communities;
    q_ij = new double[num_communities];
    // use a 1D array (k[0]) to simulate a 2D array (k[z][i])
    total_size = al->size * num_communities;
    total_bytes = sizeof(double) * total_size;
    k = new double*[num_communities];
    k_new = new double*[num_communities];
    k_best = new double*[num_communities];
    k[0] = new double[total_size];
    k_new[0] = new double[total_size];
    k_best[0] = new double[total_size];
    
    std::memset(k_best[0], 0, total_bytes);
    for (z = 1, j = al->size; z < num_communities; z++) {
        k[z] = k[0] + j;
        k_new[z] = k_new[0] + j;
        k_best[z] = k_best[0] + j;
        j += al->size;
    }
    
    while (loop-- > 0) {
        
        iter = 0;
        
        // random initial condition
        for (i = 0; i < total_size; i++) {
            k[0][i] = rand() * RAND_F + 1;
        }
        
        while (1) {
            double* p, *p_end;
            
            // calculate kappa[z] and take the reciprocal
            std::memset(kappa, 0, kappa_bytes);
            for (z = 0; z < num_communities; z++) {
                p = k[z];
                p_end = p + al->size;
                while (p < p_end)
                    kappa[z] += *p++;
                kappa[z] = 1.0 / kappa[z];
            }
            
            // go through all edges with colors
            std::memset(k_new[0], 0, total_bytes);
            for (i = 0; i < al->size; i++) {
                degree_i = al->degree[i];
                for (jj = 0; jj < degree_i; jj++) {
                    j = al->neighbors[i][jj];
                    if (i <= j) {
                        d = 0.;
                        for (z = 0; z < num_communities; z++)
                            if (k[z][i] && k[z][j])
                                d += q_ij[z] = k[z][i] * k[z][j] * kappa[z];
                        d = 1.0 / d;
                        for (z = 0; z < num_communities; z++)
                            if (k[z][i] && k[z][j]) {
                                q_ij[z] *= d;
                                k_new[z][i] += q_ij[z];
                                k_new[z][j] += q_ij[z];
                            }
                    }
                }
            }
            // set to 0 for those below EPSILON
            for (i = 0; i < total_size; i++)
                if (k_new[0][i] && k_new[0][i] < EPSILON)
                    k_new[0][i] = 0;
            
            // compare k and k_new
            for (i = 0; i < total_size; i++)
                if (fabs(k[0][i] - k_new[0][i]) > EPSILON)
                    break;
            // break loop if they are close enough
            if (i == total_size) {
                break;
            }
            else {
                k_swap = k;
                k = k_new;
                k_new = k_swap;
            }
            if (++iter > 100000) {
                std::cout << "Warning: iteration limit reached (vertex id: "
                    << al->mapping[al->size - 1] << ")" << std::endl;
                break;
            }
        }
        
        // calculate the Log Likelihood
        log_p = 0.;
        for (i = 0; i < al->size; i++) {
            degree_i = al->degree[i];
            for (jj = 0; jj < degree_i; jj++) {
                j = al->neighbors[i][jj];
                if (i <= j) {
                    d = 0.;
                    for (z = 0; z < num_communities; z++) {
                        if (k[z][i] && k[z][j])
                            d += q_ij[z] = k[z][i] * k[z][j] * kappa[z];
                    }
                    if (i < j) log_p += log(d) * 2;
                    else log_p += log(d);
                }
            }
        }
        for (i = 0; i < total_size; i++)
            log_p -= k[0][i];
        
        // check if the new result is better
        if (log_p > max_log_p) {
            max_log_p = log_p;
            std::copy(&k[0][0], &k[0][total_size], k_best[0]);
        }
    }
    
    delete[] k[0];
    delete[] k_new[0];
    delete[] k;
    delete[] k_new;
    delete[] kappa;
    delete[] q_ij;
    
    // normalize k
    for (i = 0; i < al->size; i++) {
        auto inverse_degree_i = 0.;
        for (z = 0; z < num_communities; z++)
            inverse_degree_i += k_best[z][i];
        inverse_degree_i = 1.0 / inverse_degree_i;
        for (z = 0; z < num_communities; z++)
            k_best[z][i] *= inverse_degree_i;
    }
    
    std::unique_ptr<double[]> result (k_best[0]);
    delete[] k_best;
    return result;
}

int main(int argc, char* argv[])
{
    auto path_input = "network.dat";
    auto path_output = "partial_communities.dat";
    unsigned num_communities_min = 5;
    unsigned num_communities_max = -1;
    unsigned num_communities_every = 0;
    unsigned repeat = 1;
    double threshold_belonging = 0.2;
    double threshold_merging = 0.3;
    double threshold_clustering = 0.;
    unsigned threshold_min_degree = 0;
    
    int ch;
    while ((ch = getopt(argc, argv, "i:o:n:u:e:r:b:m:c:d:")) != -1)
        switch (ch) {
            case 'i':
                path_input = optarg;
                break;
            case 'o':
                path_output = optarg;
                break;
            case 'n':
                num_communities_min = std::atoi(optarg);
                break;
            case 'u':
                num_communities_max = std::atoi(optarg);
                break;
            case 'e':
                num_communities_every = std::atoi(optarg);
                break;
            case 'r':
                repeat = std::atoi(optarg);
                break;
            case 'b':
                threshold_belonging = std::atof(optarg);
                break;
            case 'm':
                threshold_merging = std::atof(optarg);
                break;
            case 'c':
                threshold_clustering = std::atof(optarg);
                break;
            case 'd':
                threshold_min_degree = std::atoi(optarg);
                break;
            default:
                return 1;
        }
    
    auto global = AdjacencyList::LoadFromBinaryFile(path_input);
    if (!global) {
        std::cout << "Failed to load the network. Aborted." << std::endl;
        return 1;
    }
    // id mapping is not used, save memory
    global->mapping.clear();
    
    std::ofstream f (path_output, std::ofstream::binary);
    
    timeval ms1, ms2;
    gettimeofday(&ms1, NULL);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (unsigned id = 0; id < global->size; id++) {
        // screening
        // skip low degree vertices
        if (threshold_min_degree != 0 && global->degree[id] < threshold_min_degree)
            continue;
        
        // create the subnetwork of the first neighborhood of vertex id, exclude the vertex itself
        auto sub = global->FirstNeighborhoodOfVertex(id, true);
        // remove isolated vertices (degree 0 and degree 1 pairs)
        std::vector<unsigned> sub_vertices;
        sub_vertices.reserve(sub->size);
        for (unsigned i = 0; i < sub->size; i++)
            if (sub->degree[i] > 1 || (sub->degree[i] == 1 && sub->degree[sub->neighbors[i][0]] > 1))
                sub_vertices.push_back(sub->mapping[i]);
        // add the ego itself
        sub_vertices.push_back(id);
        sub = global->CreateSubAdjacencyList(sub_vertices);
        
        // calculate the number of communities to be found
        unsigned num_communities = num_communities_every > 0 ? std::min(num_communities_max, sub->size / num_communities_every + num_communities_min) : num_communities_min;
        
        // calculate k
        auto k = find_communities(sub, num_communities, repeat);
        
        // extract information in k to construct partial communities
        std::vector<std::vector<unsigned>> raw_communities (num_communities);
        for (unsigned offset = 0, z = 0; z < num_communities; offset += sub->size, z++) {
            for (unsigned i = 0; i < sub->size - 1; i++)
                if (k[offset+i] > threshold_belonging)
                    raw_communities[z].push_back(i);
        }
        
        // merge similar partial communities
        auto communities = merge_similar_partial_communities(raw_communities, threshold_merging);
        
        for (auto& c : communities) {
            // the ego vertex is a default member of all its partial communities
            c.push_back(sub->size - 1);
            // screening - remove false communities with clustering coefficient lower than a threshold
            if (threshold_clustering > 0) {
                auto sub_c = sub->CreateSubAdjacencyList(c);
                if (2*(sub_c->num_edges - sub_c->size + 1.0)/((sub_c->size - 1)*(sub_c->size - 2)) < threshold_clustering)
                    c.clear();
            }
            // id mapping
            for (auto& t : c) t = sub->mapping[t];
            // sort
            std::sort(c.begin(), c.end());
        }
        
        // save communities
        unsigned buffer_size = 0;
        for (auto& c: communities)
            buffer_size += c.size() + 2;
        std::vector<unsigned> buffer;
        buffer.reserve(buffer_size);
        for (auto& c: communities)
            if (c.size() > 2) {
                buffer.push_back(id);
                buffer.push_back(c.size());
                buffer.insert(buffer.end(), c.begin(), c.end());
            }
        if (buffer.size() > 0) {
            #pragma omp critical
            f.write((char*)buffer.data(), sizeof(unsigned) * buffer.size());
        }
    }
    f.close();
    
    gettimeofday(&ms2, NULL);
    std::cout << (ms2.tv_sec - ms1.tv_sec) + (ms2.tv_usec - ms1.tv_usec) * 1e-6 << " seconds." << std:: endl;
}
