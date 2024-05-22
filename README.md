# Introduction

The new version of Lagrange, a tool to infer ancestral ranges and splits using the DEC model. For a detailed
explanation of the please see our publication in Systematic Biology.

# Building 

First, this project relies on submodule to track its depenancies. This means that in order to clone the repository, the
following command should be used

```
git clone https://github.com/computations/lagrange-ng --recursive --depth=1
```

There is a helpful makefile which will generate and build the cmake directory. So, you can just run `make` to build.
Lagrange-NG uses CMake to configure the build, and so options for building (including options such as `Debug` or
`Release` builds) can be configured using tools which interact with cmake configurations, such as `ccmake`. There are
two important build options for the purposes of building Lagrange-NG:

- `ENABLE_MKL`: if set to `ON`, will attempt to find an Intel Math Kernel Library install and use that for linear
  algebra computation. If it is set to `OFF`, OpenBLAS will be used instead.
- `MKL_PREFIX_DIR`: Implies `ENABLE_MKL=ON`, and can be used to point to the MKL install directory.

## Troubleshooting

Some users have reported issues with CMake not finding NLOpt properly. In this case, I would recommend installing NLOpt
as a library via `apt` or `brew`. If this fails to fix the issue, please let me know.

# Running

In order to run Lagrange-NG, three files are needed:

1. A phylogeny in newick format with branch lengths specified in time,
2. An absence presence matrix for the ranges of each of the taxa in the phylogeny as a phylip file and,
3. A control file, which contains the path to both the phylogeny and the matrix file.

An example of a control file is

```
datafile = matrix.phy
treefile = phylogeny.nwk
areanames = Region1 Region2 Region3 Region4
states
workers = 4
threads-per-worker = 1
```

The options for the control file are detailed below. To run Lagrange-NG, the following command should be used

```
lagrange-ng <CONTROL FILE>
```

After computation, there will be 3 results files: `<treefile>.nodes.tre`, `<treefile>.results.json`, and
`<treefile>.scaled.tre`. The `<treefile>.nodes.tre` contains the input tree with nodes labeled according to the
`<treefile.results.json` indices. This is so that the entries in the results file can be matched to the inner nodes of
the tree. The `<treefile>.results.json` file contains the

## Example

Some example files have been provided in the `example` dir. Specifically, there is an example configuration file
`example/example.conf` with the contents

```
treefile = example/example.nwk
datafile = example/example.phy
areanames = RA RB RC RD RE
states
workers = 1
```

As well as the example tree and the example alignment in `example/example.nwk` and `example/example.phy`:

```
# example/example.nwk
((g:2.33641,f:2.33641):0.128694,((e:0.507286,d:0.507286):1.75297,((c:1.93387,(b:1.03499,a:1.03499):0.898877):0.13294,(j:0.210444,(i:0.203881,h:0.203881):0.0065634):1.85637):0.193449):0.20484);

# example/example.phy
10 5
g 00110
f 11100
e 01110
d 10010
c 01011
b 10110
a 10101
j 11011
i 00100
h 10100
```

To run an analysis with these files, simply run in the terminal the command

```
./bin/lagrange-ng example/example.conf
```

If Lagrange-NG has been built correctly, we should see the following output (or something similar):

```
[  0.00s] Reading tree...
[  0.00s] Tree number 1 has 10 leaves
[  0.00s] Reading data...
[  0.00s] Checking data...
[  0.00s] Running analysis...
[  0.00s] Using Pade's method for EXPM computation
[  0.00s] Starting Workers
[  0.00s] Making Worker#1
[  0.00s] Waiting for Workers to finish
[  0.00s] Initial LLH: -66.235818
[  0.01s] Current LLH: -66.235818
[  0.01s] Current LLH: -37.599879
[  0.02s] Current LLH: -31.542749
[  0.02s] Current LLH: -31.447068
[  0.03s] Current LLH: -31.428544
[  0.04s] Current LLH: -31.425485
[  0.04s] Current LLH: -31.424314
[  0.04s] Current LLH: -31.424298
[  0.05s] Current LLH: -31.424297
[  0.05s] Final LLH: -31.424296
[  0.05s] Computing reverse operations
[  0.05s] Computing ancestral states
[  0.05s] Computing ancestral splits
[  0.05s] Writing results to example/example.nwk.results.json
[  0.06s] Writing node annotated tree to example/example.nwk.nodes.tre
[  0.06s] Writing node annotated tree to example/example.nwk.scaled.tre
[  0.06s] Analysis took 0.056379s
```

After Lagrange-NG is finished, there should be 3 new files in `example`:

- `example.nwk.nodes.tre` the tree with internal node labels as node ids for identification with nodes in the results
  file;
- `example.nwk.results.json`, which is the file containing all the results; and
- `example.nwk.scaled.tre`, which contains the scaled tree that Lagrange-NG used to infer model parameters.

## Control File Options

- `treefile`: Path to the newick file containing the phylogeny. Required.
- `datafile`: Path to the matrix in either PHYLIP or FASTA format. In both
  cases, leading zeros are significant, I.E. `0011` is not the same as `11`.
- `areanames`: A space separated list of region names, which will be used for reference in the output. Required.
- `states`: Compute the ancestral states after fitting model parameters.
- `splits`: Compute ancestral splits after fitting model parameters.
- `workers`: Number of workers to use for execution. The optimal number of workers depends on both the size of the
  phylogeny, and the number of regions in the matrix. However, except for small problems (low number of regions and a
  small phylogeny), setting the number of workers to the number of available physical cores will give good results.
- `threads-per-worker`: Number of threads allocated to each worker. Almost always should be left unspecified. If you
  have a large number of cores (over 16), marginally faster runtimes can be obtained by setting this to `2`.
- `dispersal`/`extinction`: These two should be set together. The function depends on the execution mode. In `optimize`,
  this is the starting value for the optimization routine. In `evaluate`, computes the likelihood and (optionally) 
  ancestral states/splits for the given parameter values.
- `lh-epsilon`: Change the stopping criterion for the optimization step. Default is `1e-9`.
- `maxareas`: Only consider ranges with a number of areas equal to or less than the specified number of `maxareas`. For
  example, if the `maxareas` is set to 4 then only ranges with 4 regions or less are considered, even if the total
  number of ranges is greater.
- `expm-mode`: May be set to: `krylov`, `pade`, or `adaptive`. Set to `adaptive` by default. Generally, this shouldn't
  set. See the following section for more.
- `mode`: May be set to `optimize` or `evaulate`. The `optimize` mode will optimize the model parameters and then
  compute the ancestral states/splits. The `evaluate` mode will _not_ optimize the model parameters before computing
  ancestral state/splits. The `dispersal`/`extinction` options should be set if a particular set of parameter values
  should be used.
- `mrca`: Specify an interior node by a list of tips. The node specified is the MRCA of those tips. For example, when
  using the tree `(a,(b, c)1)2`, `mrca foo = b c` will refer to node `1`, and `mrca bar = a c` will refer to node `2`.
  Use this option to specify interior nodes for other options.
- `fossil`: Specify constraints on inner nodes. There are currently two different modes supported:
  - `node`: Constrain a node to contain at least the areas specified. For example, `fossil = node foo 011` will constrain
    the node `foo` (specified with `mrca`) to the distributions `111` and `011`.
  - `fixed`: Constrain a node to be exactly the distribution specified. For example, `fossil = fixed foo 011` will
    constrain the node `foo` to the distribution `011`.
- `logfile`: Specify a file to log the output of the program (the log messages).

## Expm Modes

Lagrange-NG implements 2 methods of computing the matrix exponential, which can have fairly significant implications for
results. The first mode is `pade`, which is the much slower but (in my experience) more accurate method of computing the
exponential. The second method is `krylov`, which uses Krylov subspaces, and is often 2 orders of magnitude faster than
the `pade` mode. However, it can be numerically unstable for this particular application.

Lagrange-NG also implements a third mode, `adaptive`, which attempts to combine the strengths of the `krylov` and `pade`
modes. By default, `adpative` uses `krylov` mode for matrix exponential computation, but watches for easily detectable
errors. If such an error is found, then the computation for that particular matrix is reverted to the `pade` mode. In my
experience, this mode fixes all the numerical errors that I could generate. I would recommend to use `adaptive` in all
cases, except when the results seem very suspect. In this case, the `pade` mode should be used. Please note that the
`pade` mode is very slow compared to `adaptive` or `krylov`, so only use this when a genuine numerical error is
suspected.

## Results JSON file

The results file (`<treefile>.results.json`) contains the results from the optimization and ancestral state/split
computation. There are 3 top level keys to the file:

- `attributes`: This contains an object with information about the dataset, such as the number of regions and taxa.
  Additionally, it contains the node tree  in `node-tree`.
- `node-results`: Contains a list of the results for each node. Discussed in detail later.
- `params`: Contains the inferred values for `dispersion` and `extinction`.

As mentioned, `node-results` contains a list of results for each node. Each entry of the list has the following keys and
values:

- `number`: Contains the node number, as specified in the
  `<treefile>.nodes.tre`. When an `mrca` is specified in the config file, the
  name of the MRCA is used instead.
- `states` (optionally): Contains the results as a list from the ancestral state computation with the following keys:
    - `distribution`: The current distribution, as a binary number, output in decimal.
    - `distribution-string`: The current distribution, but expressed using the region names from the control file.
    - `regions`: A list of the regions present in the current distribution. Each region is given as a string, with names
      from the control file.
    - `llh`: The log-likelihood for the current distribution
    - `ratio`: The likelihood weight ratio (LWR) for the current distribution
- `splits` (optionally): Contains the result as a list of the ancestral split computation with the following keys:
  - `anc-dist`, `left-dist`, `right-dist`:
    - `distribution`: The internal representation of the distribution as an integer.
    - `distribution-string`: Distribution as a string using the region names from the control file.
    - `regions`: List of regions in the distribution where the region names are taken from the control file.
  - `llh`: The log likelihood of this split.
  - `ratio`: the LWR for this split.
