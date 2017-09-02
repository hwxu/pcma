import networkx as nx
import numpy as np
import random
import itertools
import struct
import argparse

def write(f, i, size=4):
    buffer = bytearray(size)
    struct.pack_into('<I', buffer, 0, i)
    f.write(buffer)

parser = argparse.ArgumentParser()
parser.add_argument('-n', type=int, required=True, help='number of vertices')
parser.add_argument('-k', type=float, required=True, help='mean degree of background noise')
parser.add_argument('-p', type=float, required=True, help='expected intra-community edge density')
parser.add_argument('-s', type=float, required=True, help='average size of communities')
parser.add_argument('-c', type=float, required=True, help='average number of communities a vertex has')
parser.add_argument('-d', type=str, required=False, help='output folder')
args = parser.parse_args()

n = args.n
k_mean = args.k
p0 = k_mean / float(n-1)
p1 = args.p
poisson_lambda = args.s
num_communities = n / poisson_lambda * args.c

output_dir = args.d
if output_dir is None:
    output_dir = './'
elif not output_dir.endswith('/'):
    output_dir += '/'

nodes = list(np.arange(0, n, 1))
G = nx.fast_gnp_random_graph(n, p0)

community_sizes = np.random.poisson(poisson_lambda, num_communities)
communities = list()
for s in community_sizes:
    community = random.sample(nodes, s)
    communities.append(community)
    for (u,v) in itertools.combinations(community, 2):
        if G.has_edge(u,v) or np.random.random_sample() > p1:
            continue
        G.add_edge(u,v)

with open(output_dir+'generated_communities.txt', 'w') as f:
    for community in communities:
        f.write(' '.join([str(m) for m in community]))
        f.write('\n')

with open(output_dir+'network.dat','wb') as f:
    write(f, n)
    for node in G:
        write(f, G.degree(node))
        for neighbor in G.neighbors(node):
            write(f, neighbor)
