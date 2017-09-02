# Options

## STEP 1 of PCMA: Detecting partial communities

This step detects partial communities in all (or sampled) vertices’ ego networks. For more information, read Appendix A of the paper.

#### Basic Options

**`-n integer`** The number of partial communities to be found in the ego network of each vertex. It is recommended to overestimate the number as isolated vertices also need to be accommodated by some communities. The overestimation will be handled properly by merging similar partial communities. (default: 5)

**`-i filename`** Input network file. (default: "network.dat")

**`-o filename`** Output file. (default: "partial_communities.dat")

#### Advanced Options

**`-e integer`** If this option is used, the number of partial communities to be found in a vertex's ego network is the number specified by `-n` plus a number proportional to the size of the ego network. This option sets the inverse of the proportional constant.

**`-u integer`** Similar to `-n`, but specify the upper bound. This option must be used together with `-e`.

**`-b threshold`** A vertex is considered a member of a partial community if its belonging coefficient to this community is above the threshold. (default: 0.2)

**`-m threshold`** Merge partial communities of the same vertex if the overlap between two partial communities is above the threshold. (default: 0.3)

**`-c threshold`** Remove partial communities of which the clustering coefficient (take community as subnetwork) is below the threshold. This option is to remove those obviously false communities. By default, this feature is disabled. (default: 0)

**`-d min_degree`** Skip vertices with a degree lower than the specified value as it’s less likely to find partial communities in these vertices’ ego networks. By default, this feature is disabled. (default: 0)

**`-r integer`** The number of random initializations. The larger the `r`, the more accurate the result. As PCMA does not require STEP 1 to be precise, a small value is recommended. By default, the minimum value `r = 1` is used for fast analysis. (default: 1)



## STEP 2 of PCMA: Merging similar partial communities

This step repeatedly merges the most similar pair of communities to reconstruct complete communities. For more information, read Secion II B of the paper.

#### Options

**`-t threshold`** Repeatedly merge the most similar pair of partial / merged communities until the similarity is below the threshold. (default: 0.1)

**`-z threshold`** Prevent merger of two communities if the number of common members is below the threshold. This option is intended to suppress unwanted mergers of small partial communities. (default: 0)

**`-i filename`** Input file. (default: "partial_communities.dat")

**`-o filename`** Output file. (default: "merged_communities.dat")



## STEP 3 of PCMA: Post-processing

This step cleans the noise accumulated in the merged communities to sift out real communities. The provided [post_processing.cpp](../src/post_processing.cpp) is a minimal implementation and it outputs only the IDs of each community's members. Customize it to output more useful information for in-depth analysis. For details, see Section II C of the paper.

#### Options

**`-l threshold`** Set the threshold <code>t<sub>l</sub></code> as described in the paper. (default: 10)

**`-f threshold`** Set the threshold <code>t<sub>S/l</sub></code> as described in the paper. (default: 0.1)

**`-s threshold`** Set the threshold <code>t<sub>S'</sub></code> as described in the paper. (default: 2)

**`-i filename`** Input file. (default: "merged_communities.dat")

**`-o filename`** Output file. (default: "detected_communities.txt")
