#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <unordered_map>
#include <map>
#include <vector>
#include <sys/time.h>
#include <unistd.h>

class Community {
public:
    unsigned id;
    unsigned size;
    unsigned l;
    unsigned w;
    double g;
    
    unsigned best_candidate;
    float best_fs;
    bool visited;
    
    std::vector<unsigned> members;
    std::vector<unsigned> score;
    std::vector<unsigned> origins;
    std::vector<unsigned> egos;
    std::unordered_map<unsigned, double> candidates;
    
    Community(unsigned, unsigned, unsigned);
    double CountWeightedOverlap(Community&);
    void Merge(Community&, double);
    void FindBestCandidate();
};

Community::Community(unsigned cid, unsigned ego, unsigned s): id(cid), size(s), l(1), w(0),
g(1.), best_candidate(cid), best_fs(0.), visited(false), members(s), score(s)
{
    origins.push_back(cid);
    egos.push_back(ego);
}

double Community::CountWeightedOverlap(Community& c)
{
    unsigned count{0}, i{0}, j{0};
    
    while (i < size && j < c.size) {
        if (members[i] < c.members[j])
            i++;
        else if (members[i] > c.members[j])
            j++;
        else
            count += score[i++] * c.score[j++];
    }
    return (double)(2*count) / (w*c.l + c.w*l);
}

void Community::Merge(Community& c, double fs)
{
    std::vector<unsigned> new_members, new_score;
    unsigned i{0}, j{0};
    while (i < size && j < c.size) {
        if (members[i] == c.members[j]) {
            new_members.push_back(members[i]);
            new_score.push_back(score[i++] + c.score[j++]);
        }
        else if (members[i] < c.members[j]) {
            new_members.push_back(members[i]);
            new_score.push_back(score[i++]);
        }
        else {
            new_members.push_back(c.members[j]);
            new_score.push_back(c.score[j++]);
        }
    }
    while (i < size) {
        new_members.push_back(members[i]);
        new_score.push_back(score[i++]);
    }
    while (j < c.size) {
        new_members.push_back(c.members[j]);
        new_score.push_back(c.score[j++]);
    }
    size = new_members.size();
    members = std::move(new_members);
    score = std::move(new_score);
    origins.insert(origins.end(), c.origins.begin(), c.origins.end());
    egos.insert(egos.end(), c.egos.begin(), c.egos.end());
    
    g = (l-1)*w*g + (c.l-1)*c.w*c.g + (w*c.l + c.w*l)*fs;
    w += c.w;
    l += c.l;
    g = g / (w*(l-1));
}

void Community::FindBestCandidate()
{
    if (candidates.size() == 0) {
        best_candidate = id;
        best_fs = 0.;
        return;
    }
    auto it = std::max_element(candidates.begin(), candidates.end(),
                               [](std::pair<unsigned, double> a,
                                  std::pair<unsigned, double> b)
                               { return a.second < b.second; });
    best_candidate = it->first;
    best_fs = it->second;
}

std::unordered_map<unsigned, std::vector<unsigned>> node_communities;
std::vector<Community*> communities;

void load(const char* filename)
{
    std::ifstream f (filename, std::ios::binary);
    if (!f.is_open()) {
        std::cout << "Can't open " << filename << ".";
        return;
    }
    unsigned cid{0}, id, size;
    while (!f.eof()) {
        f.read((char*)&id, sizeof(unsigned));
        f.read((char*)&size, sizeof(unsigned));
        Community* c = new Community(cid, id, size);
        for (unsigned i = 0; i < size; i++) {
            // members must be in asc order
            f.read((char*)&c->members[i], sizeof(unsigned));
            c->score[i] = 1;
            c->w++;
        }
        node_communities[id].push_back(cid);
        communities.push_back(c);
        cid++;
    }
    f.close();
    std::cout << "Number of partial communities loaded: " << cid << std::endl;
}

void save(const char* filename)
{
    std::ofstream f;
    f.open(filename, std::ios::binary);
    for (unsigned i = 0; i < communities.size(); ++i) {
        // check whether the community has been merged
        if (communities[i]->id != i) continue;
        Community& c = *communities[i];
        f.write((char*)&c.id, 4);
        f.write((char*)&c.l, 4);
        f.write((char*)&c.g, 8);
        f.write((char*)&c.size, 4);
        f.write((char*)c.egos.data(), 4 * c.l);
        for (unsigned j = 0; j < c.size; ++j) {
            f.write((char*)&c.members[j], 4);
            f.write((char*)&c.score[j], 4);
        }
    }
    f.close();
}

void update_candidates(Community& c, double threshold_fs, double threshold_f0x2, bool refresh)
{
    // possible candidates' ids
    std::set<unsigned> cids;
    for (auto i : c.members)
        for (auto j : node_communities[i])
            cids.insert(communities[j]->id);
    if (refresh) {
        for (auto& it : c.candidates) cids.insert(it.first);
        c.best_candidate = c.id;
        c.best_fs = 0.;
    }
    // remove ids that are already in candidates
    else for (auto& it : c.candidates) cids.erase(it.first);
    // remove self if in cids
    cids.erase(c.id);
    
    // calculate fs, update candidates
    for (auto cid : cids) {
        auto& c2 = *communities[cid];
        auto fs = c.CountWeightedOverlap(c2);
        
        // add candidate or update fs
        if (fs > threshold_fs &&
            fs*(c.w*c2.l+c2.w*c.l)/std::max(c.l, c2.l) >= threshold_f0x2) {
            c.candidates[c2.id] = fs;
            if (fs > c.best_fs) {
                c.best_fs = fs;
                c.best_candidate = c2.id;
            }
            c2.candidates[c.id] = fs;
            if (fs > c2.best_fs) {
                c2.best_fs = fs;
                c2.best_candidate = c.id;
            }
            // if c was c2's best candidate, find the new best candidate for c2
            else if (refresh && c2.best_candidate == c.id)
                c2.FindBestCandidate();
        }
        // remove candidate
        else if (refresh && c.candidates.count(c2.id) > 0) {
            c.candidates.erase(c2.id);
            c2.candidates.erase(c.id);
            if (c2.best_candidate == c.id)
                c2.FindBestCandidate();
        }
    }
}

void merge(Community& c, Community& c2, double fs, double threshold_fs, double threshold_f0x2)
{
    c.Merge(c2, fs);
    
    for (auto& it : c2.candidates) {
        auto& c3 = *communities[it.first];
        c3.candidates.erase(c2.id);
        if (c3.best_candidate == c2.id)
            c3.FindBestCandidate();
    }
    for (auto cid : c2.origins)
        communities[cid] = &c;
    
    delete &c2;
    
    update_candidates(c, threshold_fs, threshold_f0x2, true);
}

int main(int argc, char* argv[])
{
    auto path_input = "partial_communities.dat";
    auto path_output = "merged_communities.dat";
    double threshold_fs = 0.1;
    double threshold_f0x2 = 0.;
    
    int ch;
    while ((ch = getopt(argc, argv, "i:o:t:z:")) != -1)
        switch (ch) {
            case 'i':
                path_input = optarg;
                break;
            case 'o':
                path_output = optarg;
                break;
            case 't':
                threshold_fs = std::atof(optarg);
                break;
            case 'z':
                threshold_f0x2 = std::atof(optarg) * 2;
                break;
            default:
                return 1;
        }

    
    timeval ms1, ms2;
    gettimeofday(&ms1, NULL);
    
    load(path_input);
    
    bool flag = true;
    while (flag) {
        flag = false;
        Community* c, *c2;
        for (unsigned i = 0; i < communities.size(); i++) {
            c = communities[i];
            if (c->id != i) continue;
            if (!c->visited) {
                update_candidates(*c, threshold_fs, threshold_f0x2, false);
                c->visited = true;
            }
            while (c->best_candidate != c->id) {
                c2 = communities[c->best_candidate];
                if (!c2->visited) {
                    update_candidates(*c2, threshold_fs, threshold_f0x2, false);
                    c2->visited = true;
                }
                if (c2->best_candidate == c->id ||
                    c2->best_fs == c->best_fs) {
                    merge(*c, *c2, c->best_fs, threshold_fs, threshold_f0x2);
                    flag = true;
                }
                else c = c2;
            }
        }
    }

    save(path_output);
    
    gettimeofday(&ms2, NULL);
    std::cout << (ms2.tv_sec - ms1.tv_sec) + (ms2.tv_usec - ms1.tv_usec) * 1e-6 << " seconds." << std:: endl;
}
