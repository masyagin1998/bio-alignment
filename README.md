# needleman-wunsch
Implementation of Needleman-Wunsch bioinformatics algorithm for comparing biological sequences.

### Tech

Algorithm is coded in pure `C89` without any dependencies.

### Installation

`needleman-wunsch` requires only `C89`-compatible compiler and `make` utility.

```sh
$ cd needleman-wunsch
$ make
$ ./bin/nw --help
$ ./bin/nw -i data/in.fasta -o out.fasta -s blosum62 -g -5
```

### Short description

`needleman-wunsch` algorithm is an algorithm used in bioinformatics to align protein or nucleotide sequences, using `dynamic programming`.
