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

### Tests

The `in.fasta` and` in1.fasta`, `in2.fasta` files, used for testing the utility, are located in the` data` folder.

The `in.fasta` test can be run as follows:
```sh
$ ./bin/nw -i data/in.fasta
```

The `in1.fasta`, `in2.fasta` test can be run as follows:
```sh
$ ./bin/nw -i data/in1.fasta,data/in2.fasta
```
