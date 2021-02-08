# Tarantula

[![Latest GitHub release](https://img.shields.io/github/release/lbcb-sci/tarantula.svg)](https://github.com/lbcb-sci/tarantula/releases/latest)
[![Build status for c++/clang++](https://travis-ci.com/lbcb-sci/tarantula.svg?branch=master)](https://travis-ci.com/lbcb-sci/tarantula)

Tarantula is a Hi-C scaffolder.

## Usage
To build tarantula run the following commands:

```bash
git clone https://github.com/lbcb-sci/tarantula && cd tarantula && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

which will create tarantula executable and unit tests (running `make install` will install the executable to your system). Running the executable will display the following usage:

```bash
usage: tarantula [options ...] <target> <sequences> [<sequences> ...]

  # default output is to stdout
  <target>/<sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage;
```

#### Build options
- `tarantula_build_tests`: build unit tests

#### Dependencies
- gcc 4.8+ | clang 3.5+
- cmake 3.11+
- zlib 1.2.8+

###### Hidden
- lbcb-sci/ram 2.0.0
- rvaser/bioparser 3.0.13
- (tarantula_test) google/googletest 1.10.0

## Acknowledgment
This work has been supported in part by the Genome Institute of Singapore (A\*STAR).
