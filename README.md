# PCMA - Partial Community Merger Algorithm

PCMA is an efficient algorithm for detecting community structure in large-scale networks with millions / billions of vertices. Especially suitable for analyzing huge social networks.

For more information, read the paper: 

E. H. W. Xu and P. M. Hui, Network Science (2017). [doi:10.1017/nws.2017.32](https://doi.org/10.1017/nws.2017.32), [arXiv:1509.00556](https://arxiv.org/abs/1509.00556)

## Installation

```
$ git clone https://github.com/hwxu/pcma.git
$ cd pcma
$ make
```

Three executables are generated: `detect`, `merge`, and `post-processing`.

Remove the option `-fopenmp` in [Makefile](Makefile) if your compiler does not support OpenMP.

## Quick Start

Generate a synthetic network with community structure using [benchmark.py](tools/benchmark.py):

```
$ mkdir output
$ python3 ./tools/benchmark.py -n 10000 -k 20 -p 0.3 -s 40 -c 2 -d ./output
```

The generated `network.dat` is in binary format. Read [format.md](docs/format.md) for detail information.

PCMA consists of three steps:
1. Detect the partial communities in the ego network of each vertex.
2. Merge partial communities that are parts of the same community to reconstruct complete communities.
3. Clean the merged communities to sift out real communities.

```
$ ./detect -n 5 -i ./output/network.dat -o ./output/partial_communities.dat
$ ./merge -t 0.1 -i ./output/partial_communities.dat -o ./output/merged_communities.dat
$ ./post-processing -l 10 -f 0.1 -i ./output/merged_communities.dat -o ./output/detected_communities.txt
```

For advanced usage, see [options.md](docs/options.md).

The detected communities are saved in `detected_communities.txt`, with each line showing the member IDs of a community. To evaluate the accuracy, compare the result to the ground-truth `generated_communities.txt` with an evaluation measure, such as the Normalized Mutual Information or Omega Index.

## License

PCMA is licensed under the GNU General Public License v2.0 or any later version. See [LICENSE](LICENSE) for the full license text.
