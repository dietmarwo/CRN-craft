# Third-party notices

## Exponential Ziggurat sampler

The exponential tables and algorithm in `crates/crn-craft/src/ssa/ziggurat*`
are adapted from Julia 1.13.3, `stdlib/Random/src/normal.jl`. The adaptation
uses CRN-craft random streams, explicit work/cancellation checks, iterative
rejection and strict positive-wait checks; it does not promise Julia RNG replay.

Upstream: https://github.com/JuliaLang/julia/blob/v1.13.3/stdlib/Random/src/normal.jl

MIT License

Copyright (c) 2009-2025: Jeff Bezanson, Stefan Karpinski, Viral B. Shah, and other contributors: https://github.com/JuliaLang/julia/contributors

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

end of terms and conditions

[Julia's LICENSE.md ends with a pointer to its own THIRDPARTY.md, which
covers other software in the Julia repository; none of it is used here.]

## Optional parameter-optimization references

The isolated `validation/parameter-optimization/worker` workspace uses rebop
0.10.0 as an optional benchmark reference. Its generated wrappers invoke the
upstream macro API; it is not a dependency or implementation of the default
CRN-craft engine. The upstream MIT notice is preserved in
`validation/parameter-optimization/REBOP_LICENSE.txt`. Cargo.lock pins its
reference dependencies. Python fcmaes-rust/rebop dependencies are separately
listed in that directory's requirements file.

## Rust crates in the Python extension

The extension module of the Python package (`crn_craft._native`, built from
`bindings/python` with its own `Cargo.lock`) is built from the Rust crates
that lockfile pins; procedural macros and build scripts among them run
only at build time. Besides the crates of the core engine they are PyO3 (the
Python binding) and faer with its dependencies (gemm, nano-gemm, pulp,
dyn-stack, reborrow and others), which supply the sparse LU factorization
of the implicit projection integrator and of stationary laws. Every crate
and its declared license is listed in
`validation/dependency-policy-python.json`, and
`validation/check_dependencies.py` fails when the lockfile holds a crate
or a license outside that list. The licenses are MIT, Apache-2.0
(including Apache-2.0 WITH LLVM-exception), BSD-2-Clause, Zlib, Unicode-3.0
and the Unlicense, most crates offering a choice among them; none is a
copyleft license. faer is built without its default features (no rayon
thread pool) and needs no system BLAS or LAPACK library. The upstream
license texts are in each crate's published source.

## Reference data of the guide case studies

The MIT license of this repository covers the code and text written for
CRN-craft. The files below are numerical tables derived from other authors'
published data. They are not covered by the MIT license and remain under
the terms of their sources. They are included so that the checks and plots
of tutorials S4 to S7 and O3 and the figures of optimization lessons C6 to C9
run without network access. No code from these sources is included. Where
an example has a `fetch.py`, it rebuilds the tables from the stated
revision; the lessons download their sources by commit with
`tutorials/optimization/publications/sources.py`.

| Files | Derived from | License of the source |
| --- | --- | --- |
| `examples/mtdna_spread/reference-1c.csv`, `reference-2e.csv` | Simulated figure data of Insalata, Hoitzing, Aryaman and Jones (2022), doi:10.1073/pnas.2122073119: https://github.com/ferdinando17/Survival-of-the-Densest at `45462390b4c27cf2ffe4860c4f526e67d49ad6c8` | No license file at that revision |
| `examples/mtdna_network/conditions.csv`, `reference-variance.csv` | Parameter tables and simulation summaries of Aryaman, Bowles, Jones and Johnston (2019), doi:10.1534/genetics.119.302423: https://github.com/ImperialCollegeLondon/MitoNetworksGenetics at `8f158fb7a61eaaa622d6e715d654c8bdfb44db89` | GPL-3.0 |
| `examples/mtdna_bottleneck/reference-table1.csv`, `reference-table2.csv`, `reference-hb.csv`, `reference-figure3.csv` | Tables 1 and 2, the Figure 4 source data and readings of Figure 3 of Johnston et al. (2015), doi:10.7554/eLife.07464 | CC BY 4.0 |
| `examples/two_stage_expression/reference-figure1d.csv` | Marker heights read from Figure 1D of Shahrezaei and Swain (2008), doi:10.1073/pnas.0803850105 | Readings of a published figure |
| `examples/yeast_rna_fit/source-counts.csv`, `source-sites.csv`, `source-audit.csv`, `source-parameters.json`, `reference-site-activity.csv`, `reference-site-distribution.csv` | Cell counts, parameter vectors and figure series of Munsky, Li, Fox, Shepherd and Neuert (2018), doi:10.1073/pnas.1804060115: https://github.com/MunskyGroup/Munsky_et_al_2018 at `1aa0b5ce25d21484c7c4e8a894477b1719209ce7` | GPL-3.0 |
| `tutorials/optimization/results/C6/distributions.csv` (column `observed`) | Cell counts of Gómez-Schiavon, Chen, West and Buchler (2017), doi:10.1186/s13059-017-1297-9: https://github.com/mgschiavon/BayFish at `70d1b685da4b4da8aa8a64934e200a5e3593791d` | GPL-3.0 |
| `tutorials/optimization/results/C7/archive-groups.csv`, `archive-front.csv` | Benchmark results of Felton, Rittig and Lapkin (2021), doi:10.1002/cmtd.202000051: https://github.com/sustainable-processes/summit at `1de682d05e97adcfb96cd8376e876cef2d6160d3` | MIT |
| `tutorials/optimization/results/C8/burst-sweep.csv` (column `archive_eopt`), `native-pair.csv` and `domain-expansion.csv` (columns `a_archive`, `b_archive`) | Result files of Fox and Munsky (2019), doi:10.1371/journal.pcbi.1006365: https://github.com/Munsky/FSPFIM_2018 at `2e2f76fdeeb0d9290b48aa981210422d254b07c2` | No license file at that revision |
| `tutorials/optimization/results/C9/measurements.csv`, `nll.csv` (row `archive`) | Measurements of Boehm et al. (2014), doi:10.1021/pr5006923, and archived simulations: https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab at `49900951a5594a757ad5e651684c7f0a59272881`, problem `Boehm_JProteomeRes2014` | BSD-3-Clause for the collection |
