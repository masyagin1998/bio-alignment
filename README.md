# bio-alignment
 Implementation of Needleman-Wunsch, Smith-Waterman, Hirschberg and affine bioinformatics algorithms for alighning biological sequences.

### Tech

Algorithm is coded in pure `C89` without any dependencies.

### Installation

`bio-alignment` requires only `C89`-compatible compiler and `make` utility.

```sh
$ cd bio-alignment
$ make
$ ./bin/bio-alignment --help
$ ./bin/bio-alignment -i data/in.fasta -o out.fasta -s blosum62 -g -5 -a nw
```

### Tests

The `in1.fasta`, `in2.fasta`, `in1_1.fasta`, `in1_2.fasta`, `in2_1.fasta`, `in2_1.fasta` files, used for testing the utility, are located in the` data` folder.

The `in1.fasta` test can be run as follows:
```sh
$ ./bin/bio-alignment -i data/in1.fasta
```

The `in1_1.fasta`, `in1_2.fasta` test can be run as follows:
```sh
$ ./bin/bio-alignment -i data/in1_1.fasta,data/in1_2.fasta
```
