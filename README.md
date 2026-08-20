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
For fast local builds, install available dependencies through your system package
manager before configuring RNAnue. On Ubuntu, useful packages include
`libboost-program-options-dev`, `libtbb-dev`, `zlib1g-dev`, `libpng-dev`,
`libbz2-dev`, and `liblzma-dev`.
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
> not ABI-compatible with RNAnue's GCC/libstdc++ build. The presets use AUTO-first
> dependency discovery and fall back to bundled Boost.Program_options on macOS
> when the installed Boost package is not safe for the selected compiler.

#### Dependencies

RNAnue includes or uses the following dependencies:

- [Boost.Program_options](https://github.com/boostorg/boost) (installed package preferred; bundled fallback v1.86.0)
- [Segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (v0.3.4)
- [SeqAn](https://github.com/seqan/seqan3) (v3.3.0)
- [Matplot++](https://github.com/alandefreitas/matplotplusplus) (installed package preferred; pinned bundled fallback)

RNAnue uses an AUTO-first dependency policy by default. CMake first tries
installed packages from system paths, Conda, Homebrew, MacPorts, and paths passed
through `CMAKE_PREFIX_PATH`; missing dependencies then use pinned bundled
fallbacks where available.

The following dependencies will be used if present, otherwise they will be fetched
in AUTO mode (internet connection required):

- [htslib](https://github.com/samtools/htslib.git) (v1.20)
- [Vienna Package](https://www.tbi.univie.ac.at/RNA/#binary_packages) (v2.6.4)
- [oneTBB](https://github.com/oneapi-src/oneTBB) (v2022.0.0)
- [zlib](https://github.com/madler/zlib) (v1.3.1, only when needed by bundled htslib)

Dependency provider controls:

```bash
cmake --preset release -DRNANUE_DEPENDENCY_PROVIDER=AUTO
cmake --preset release -DRNANUE_DEPENDENCY_PROVIDER=SYSTEM
cmake --preset release -DRNANUE_DEPENDENCY_PROVIDER=BUNDLED
```

`AUTO` is the default. `SYSTEM` never downloads fallbacks and fails with install
guidance if a required dependency is missing. `BUNDLED` forces pinned bundled
fallbacks where RNAnue provides one. Individual dependencies can override the
global default with `RNANUE_BOOST_PROVIDER`, `RNANUE_VIENNARNA_PROVIDER`,
`RNANUE_HTSLIB_PROVIDER`, `RNANUE_TBB_PROVIDER`, `RNANUE_MATPLOT_PROVIDER`, and
`RNANUE_ZLIB_PROVIDER`.

When building inside Conda or another non-system prefix, inspect the final
runtime bindings if CMake warns about hidden implicit libraries:

```bash
ldd build/release/RNAnue | egrep 'libz|libgomp|libomp|libtbb|libRNA|libhts'
otool -L build/release/RNAnue
```

Use one dependency prefix consistently when possible, for example by passing
`-DCMAKE_PREFIX_PATH=/path/to/env -DRNANUE_DEPENDENCY_PREFIX=/path/to/env` for
Conda-style Linux builds. RNAnue uses that active prefix to prefer common runtime
libraries such as zlib, bzip2, libpng, and oneTBB from the same prefix instead of
mixing them with `/usr/lib` libraries in an unsafe runtime search path.

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

#### Coverage-shaped support metrics

The interaction table reports `no_splits` as the weighted split-read contribution assigned to an
interaction cluster; this can be fractional when a read group is shared across competing hit groups.
RNAnue also reports coverage-shaped support metrics by default. These are structural read-support
metrics: they describe how split-read evidence is distributed across the two interaction arms, but
do not by themselves establish a physical RNA interaction mechanism.

Input coverage is built from the aligned bases of each split-read fragment:

- CIGAR `M`, `=`, and `X` blocks add covered intervals.
- CIGAR `D` and `N` blocks move along the reference but add no covered bases.
- Each covered interval is weighted by the fragment's `XB` split-read contribution.
- The weighted intervals are reduced to piecewise-constant per-base coverage runs for each arm.

The reported metrics are calculated from those coverage runs:

- `total_span_bp`: length of the first merged arm plus length of the second merged arm.
- `support_per_total_bp`: `no_splits / total_span_bp`.
- `effective_coverage_span_bp`: `(sum coverage)^2 / sum(coverage^2)`, using integrated coverage
  across both arms. Evenly spread support gives a larger effective span; sharply peaked support
  gives a smaller effective span.
- `support_per_effective_bp`: `no_splits / effective_coverage_span_bp`.
- `coverage_concentration`: `1 - effective_coverage_span_bp / total_span_bp`, clamped to `[0, 1]`.
  Values near 0 indicate coverage spread across most of the merged span; values near 1 indicate
  coverage concentrated in a small part of a larger span.
- `coverage_components`: number of contiguous coverage components across both arms. A run counts as
  part of a component only when its coverage is at least `max(1.0, 0.05 * maxCoverage)`, which avoids
  counting very weak shoulders around stronger peaks.
- `arm_balance`: `min(first_arm_integrated_coverage, second_arm_integrated_coverage) /
  max(first_arm_integrated_coverage, second_arm_integrated_coverage)`. Values near 1 indicate
  balanced support on both arms; values near 0 indicate mostly one-sided support.

Coverage-shaped filtering is opt-in through `mineffdens`, `maxcovcomp`, and `minarmbal`. Each sample
also gets an aggregate weighted bedGraph of retained interaction-arm coverage for visual inspection.

#### Coverage profile assignment

`coverage_profile` is assigned from the metrics above using ordered labels:

- `low_support_density`: zero effective span or `support_per_effective_bp < 0.02`.
- `imbalanced_support`: `arm_balance < 0.25`.
- `multi_peak_refine`: more than two coverage components.
- `broad_diffuse`: span above 300 bp with `coverage_concentration <= 0.25`.
- `compact_dense`: all other retained coverage shapes.

Together, these labels flag suspicious cluster shapes after statistical scoring: broad intervals with
localized support, multi-peak merges, and interactions supported mainly by one arm.

#### Significance columns

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
