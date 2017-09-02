#include <memory>
#include <unordered_map>
#include <vector>

class AdjacencyList {

public:
    unsigned size;
    unsigned num_edges;
    std::unique_ptr<unsigned[]> degree;
    std::unique_ptr<std::unique_ptr<unsigned[]>[]> neighbors;
    std::unordered_map<unsigned, unsigned> mapping;
    
    // Constructor
    AdjacencyList(unsigned);
    // Move Constructor
    AdjacencyList(AdjacencyList&&);
    // Destructor
    virtual ~AdjacencyList();
    
    static std::unique_ptr<AdjacencyList> LoadFromBinaryFile(const char*, bool = true);
    void SaveToBinaryFile(const char*) const;
    
    std::unique_ptr<AdjacencyList> CreateSubAdjacencyList(std::vector<unsigned>&) const;
    std::unique_ptr<AdjacencyList> FirstNeighborhoodOfVertex(unsigned, bool = false) const;
};
