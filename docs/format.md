# Network Format used by PCMA

The input network is encoded into an array of 4-byte unsigned integers in binary format:

<code>
n d<sub>0</sub> [Vertex 0's d<sub>0</sub> neighbors] d<sub>1</sub> [Vertex 1's d<sub>1</sub> neighbours] ... d<sub>n-1</sub> [Vertex n-1's d<sub>n-1</sub> neighbors]
</code>
<br />

`n` is the number of vertices in the network, <code>d<sub>i</sub></code> is the degree of Vertex i. Note that vertices are numbered by consecutive integers starting from 0.

## Example

Consider a network represented by the following adjacency matrix:

```
0 1 1 0
1 0 1 1
1 1 0 0
0 1 0 0
```

The corresponding encoded array is:

```
4 2 1 2 3 0 2 3 2 0 1 1 1
```

## Converting from Edge List
A program [el_to_al_.cpp](../tools/el_to_al.cpp) is provided for converting edge list to the format used by PCMA.

```
$ g++ -std=c++11 el_to_al.cpp -o el_to_al
$ ./el_to_al -e edge_list.txt -a adjacency_list.dat -p
```

By default, vertices are relabelled by consecutive integers starting from 0 after conversion. The old IDs are appended to the converted network `adjacency_list.dat` in ascending order of the corresponding new ID. To keep the vertex ID unchanged, add option `-p`.
