CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall
VPATH = ./src

all : detect merge post-processing

detect : detect_partial_communities.cpp adjacency_list.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $^ -o $@

merge : merge.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

post-processing : post_processing.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY : clean
clean :
	-rm -f detect merge post-processing
