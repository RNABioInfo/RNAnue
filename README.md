# RNAnue - 1.0.0

[![docker-release](https://github.com/Ibvt/RNAnue/actions/workflows/docker.yml/badge.svg)](https://github.com/Ibvt/RNAnue/actions/workflows/docker.yml)

## About

RNAnue is a comprehensive tool for detecting RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.

## Installation

We recommend using the provided installer to set up RNAnue.
If your operating system is not supported, you can use RNAnue via Docker or build it from source.

### Installer

_Not available yet_

Download OS-specific installers for macOS and Debian here: [Latest Release](https://github.com/RNABioInfo/RNAnue/releases/latest).

### Docker

_Not available yet_

We provide a ready-to-use [Docker container](https://hub.docker.com/repository/docker/cobirna/rnanue) with RNAnue pre-configured.

```bash
docker pull cobirna/rnanue:latest
docker run -ti  cobirna/rnanue
```

### Singularity

_Not available yet_

The Docker container can also be used with Singularity.

```bash
singularity pull docker://cobirna/rnanue:latest
singularity exec --bind /path/to/data:/data rnanue_latest.sif RNAnue <subcall> --config /data/params.cfg
```

### Building from source

#### Prerequisites

To build RNAnue, you need `cmake (>=v3.24.0)` and GCC/G++ 14 or newer.
RNAnue uses C++23 standard-library APIs that are not fully available in GCC 13,
including `<print>`, `std::forward_like`, `std::ranges::to`, and
`std::ranges::zip_view`.
For fast local builds, install Boost.Program_options through your system package
manager before configuring RNAnue. On Ubuntu, this package is
`libboost-program-options-dev`.
If you need to compile ViennaRNA, you also need `autoconf`, `automake`, and `libtool` (see [Dependencies](#dependencies)).

#### Downloading

Clone the repository and its submodules:

```bash
git clone --recurse-submodules <git-repo-here>
cd RNAnue
git submodule update --init --recursive
```

> **Information**
> If you'd like to clone a specific branch, use the `--branch <branch-name>` flag.

#### Building and Installing

Build and install RNAnue with the following commands:

```bash
cmake --preset release
cmake --build --preset release --parallel <num-threads-here>
cmake --install build/release
```

The presets test the default C++ compiler before the main CMake language
configuration. If the default compiler is not compatible but a compatible
`g++-14` or newer compiler is on `PATH`, RNAnue selects that compiler and the
matching `gcc` automatically.

For local debugging, use:

```bash
cmake --preset debug
cmake --build --preset debug --parallel <num-threads-here>
```

To build and run the tests, use:

```bash
cmake --preset test
cmake --build --preset test --parallel <num-threads-here>
ctest --preset test
```

The presets do not require Ninja; CMake will use the default generator for your system unless you override it.

> **IMPORTANT**
> RNAnue can only be compiled with [gcc](https://gcc.gnu.org) 14 or newer (tested with v14.2.0; needs to support **C++23**)

> **IMPORTANT – MacOS**
> RNAnue will try to select Homebrew or MacPorts GCC automatically to avoid AppleClang.
> If auto-detection cannot find GCC, specify the compiler explicitly:
> `cmake --preset release -DCMAKE_CXX_COMPILER=<path-to-g++> -DCMAKE_C_COMPILER=<path-to-gcc>`
> System Boost packages on macOS are commonly built with AppleClang/libc++ and are
> not ABI-compatible with RNAnue's GCC/libstdc++ build. The `release` preset will
> use the bundled Boost fallback on macOS. If you use `debug` or `test` on macOS,
> add `-DRNANUE_BOOST_PROVIDER=BUNDLED`.

#### Dependencies

RNAnue includes or uses the following dependencies:

- [Boost.Program_options](https://github.com/boostorg/boost) (system package preferred; bundled fallback v1.86.0)
- [Segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (v0.3.4)
- [SeqAn](https://github.com/seqan/seqan3) (v3.3.0)

The `debug` and `test` presets require system Boost.Program_options and fail early if it is missing.
The `release` preset tries system Boost.Program_options first and falls back to the bundled Boost build.

The following dependencies will be used if present on the system, otherwise they will be fetched (internet connection required):

- [htslib](https://github.com/samtools/htslib.git) (v1.20)
- [Vienna Package](https://www.tbi.univie.ac.at/RNA/#binary_packages) (v2.6.4)

On Linux, [oneTBB](https://github.com/oneapi-src/oneTBB) is also fetched and built automatically if no system CMake package is found.

## Overview

![Principle](principle.png)

## Usage

### Positional Arguments

RNAnue provides different sub-calls for individual pipeline steps. These include `preprocess`,
`align`, `detect`, `analyze`. In addition, `complete` applies the whole workflow.

## Input

> **IMPORTANT** RNAnue requires the sequencing files to be in a specific folder structure. Files can also be compressed with gzip.

> **IMPORTANT** In order to process paired-end files, files must end with "\_forward.fastq" and "\_reverse.fastq", "\_R1.fastq" and "\_R2.fastq" or "\_1.fastq" and "\_2.fastq".

The root folders of the treatments (`--trtms`; required) and controls (`--ctrls`; optional) are specified accordingly. These folders contain sub-folders
with arbitrary samples that in turn contain the read files.

### Example folder structure

```text
./trtms/
    sample1
        *.fastq
    sample2
        *_R1.fastq
        *_R2.fastq
./ctrls/
    sample3
        *forward.fastq
        *reverse.fastq
    sample4
        *.fastq
```

## Parameters

RNAnue accepts parameter settings both from the command-line and through a configuration file.
For the latter, we provide a template configuration file ([params.cfg](./example/params.cfg)) that
allows to set the parameters in a more convenient fashion. This means that the call of RNAnue
is reduced to the following call.

```bash
RNAnue <sub-call-here> --config <params.cfg-here>
```

In any case, the specifying parameters over the command lines has precedence over the config file.
Boolean parameters accept explicit values on the command line, for example
`--deduplicate=false`. Options that are enabled by default also provide clearer inverse flags,
such as `--no-deduplicate`, `--no-preprocess`, `--no-trimpolyg`, `--no-maskmulticopy`,
`--no-multimap`, and `--keep-altsplice`.

## Results

In principle, the results of the analysis are stored in the specified output folder and its sub-folders
(e.g., ./preprocess, ./align, ./detect, ./analyze). RNAnue reports the split reads in SAM format, the clusters
and the RNA-RNA interactions. RNAnue reports the split reads in SAM format. Additionally, the complementarity
scores and hybridization energies are stored in the tags FC and FE, respectively. We report the clusters in a
custom format that includes the IDs of the clusters, its length, size and genomic coordinates.

### Split Reads (.SAM)

RNAnue reports the detected splits in .SAM format (RNAnue `detect`). In this file, pairs of rows represent the
split reads, consisting of the individual segments, e.g

```text
A00551:29:H73LYDSXX:1:1101:7274:10645 16 gi|170079663|ref|NC_010473.1| 3520484 22 1X51= * 0 0 AGGGGTCTTTCCGTCTTGCCGCGGGTACACTGCATCTTCACAGCGAGTTCAA * XA:Z:TTTCTGG XC:f:0.714286 XE:f:-15.6 XL:i:7 XM:i:5 XN:i:0 XR:f:0.0735294 XS:i:5 XX:i:1 XY:i:52
A00551:29:H73LYDSXX:1:1101:7274:10645 16 gi|170079663|ref|NC_010473.1| 3520662 22 11=5S * 0 0 TTCGATCAAGAAGAAC * XA:Z:GAAGAAC XC:f:0.714286 XE:f:-15.6 XL:i:7 XM:i:5 XN:i:0 XR:f:0.0735294 XS:i:5 XX:i:53 XY:i:68

```

In the following the tags are listed that are reported in the detected split reads. Please note that in the upper
segment the alignment is in reverse as done in the calculation of the complementarity to represent the 3'-5' and 5'-3'
duplex.

| tag  | description           |
| ---- | --------------------- |
| XC:f | complementarity       |
| XL:f | length of alignment   |
| XR:f | site length ratio     |
| XM:i | matches in alignment  |
| XA:Z | alignment of sequence |
| XE:f | hybridization energy  |

### Clustering results

### Interaction table

The interaction table reports `no_splits` as the weighted split-read contribution assigned to an
interaction cluster; this can be fractional when a read group is shared across competing hit groups.
RNAnue also reports coverage-shaped support density metrics. `support_per_total_bp` is the simple
weighted support divided by the merged arm span. `effective_coverage_span_bp` is derived from
weighted per-base arm coverage as `(sum coverage)^2 / sum(coverage^2)`, and
`support_per_effective_bp` divides `no_splits` by this effective span. This keeps the intuitive
support-per-base idea while making broad, permissively merged intervals visible through
`coverage_concentration`, `coverage_components`, `arm_balance`, and `coverage_profile`. These
metrics are reported by default; filtering by them is opt-in through the analyze options.
The optional coverage filters are `mineffdens`, `maxcovcomp`, and `minarmbal`.
Each sample also gets an aggregate weighted bedGraph of retained interaction-arm coverage for
visual inspection.
The `p_value` column is an analytic within-sample enrichment statistic: among the candidate
interaction clusters detected in the sample, RNAnue tests whether a cluster has more weighted
split-read support than expected from independent feature abundance estimated from contiguous
background reads. `padj_value` is the Benjamini-Hochberg adjusted value across the evaluated
candidate clusters. These values are intended for prioritization and filtering; they do not by
themselves establish biological mechanism, physical validation, or treatment-control differential
interaction.

### Testing

## Troubleshooting

contact <cobi@ibvt.uni-stuttgart.de> or create an issue
