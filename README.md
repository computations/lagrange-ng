[![Build and run tests](https://github.com/computations/lagrange-ng/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/computations/lagrange-ng/actions/workflows/cmake-multi-platform.yml)
[![Release](https://github.com/computations/lagrange-ng/actions/workflows/create-release.yml/badge.svg)](https://github.com/computations/lagrange-ng/actions/workflows/create-release.yml)

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

In order to run Lagrange-NG, three files are required:

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
`<treefile>.states.tre`. The `<treefile>.nodes.tre` contains the input tree with nodes labeled according to the
`<treefile.results.json` indices. This is so that the entries in the results file can be matched to the inner nodes of
the tree. `<treefile>.states.tre` contains the best states for each node which was specified in the config file (to
specify individual nodes, see the section "Control File Options").

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
[  0.06s] Writing scaled tree to example/example.nwk.scaled.tre
[  0.06s] Analysis took 0.056379s
```

After Lagrange-NG is finished, there should be 3 new files in `example`:

- `example.nwk.nodes.tre` the tree with internal node labels as node ids for identification with nodes in the results
  file;
- `example.nwk.results.json`, which is the file containing all the results; and
- `example.nwk.scaled.tre`, which contains the scaled tree that Lagrange-NG used to infer model parameters.

Don't forget to clean up the generated files to keep the git repository clean.

## Control File Options

- `treefile`: Path to the newick file containing the phylogeny. Required.
- `datafile`: Path to the matrix in either PHYLIP or FASTA format. In both
  cases, leading zeros are significant, I.E. `0011` is not the same as `11`.
- `areanames`: A space separated list of region names, which will be used for reference in the output. Required.
- `states`: Compute the ancestral states after fitting model parameters. When this is set, the best ancestral node
  states will be written to the file `<prefix>.states.tre` as inner node labels.
- `splits`: Compute ancestral splits after fitting model parameters. When this option is set, the best ancestral split
  _conditioned on the best ancestral state_ will be written to the file `<prefix>.splits.tre`. Splits are given as node
  labels on inner nodes, in the form `<left range> | <right range>`, and the range is given as a list of area names with
  underscores (`_`) between each region.
- `workers`: Number of workers to use for execution. The optimal number of workers depends on both the size of the
  phylogeny, and the number of regions in the matrix. However, except for small problems (low number of regions and a
  small phylogeny), setting the number of workers to the number of available physical cores will give good results.
- `threads-per-worker`: Number of threads allocated to each worker. Almost always should be left unspecified. If you
  have a large number of cores (over 16), marginally faster run times can be obtained by setting this to `2`.
- `dispersion`/`extinction`: These two should be set together. The function depends on the execution mode. In `optimize`,
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
- `period`: Specify and modify periods. None of the below options are required,
  except that periods cannot be overlapping. As the default start and end times
  for a period are `0` and `inf`, this means that if the control file specifies
  more than one period, the boundry between them must be defined with start and
  end times. Please see the section on periods for an example.
  - `<VALUE>`: declare a period with the name `<VALUE>`. Periods must be
  declared before they are modified.
  - `start`: The start time for the period. If not specified, defaults to `0.0`.
  - `end`: The end time for the period. If not specified, defaults to `inf`.
  - `matrix`: Adjustment matrix filename for the period. Used to specify
    region connectivities. If specified, the parameter `distance-penalty` will
    be inferred, which estimates the "strength" of the connectivity. Specified
    as a CSV file. Please see the section periods for details
  - `include`: A range mask to define which areas _must_ be included for this
    period. Given as a binary string.
  - `exclude`: A range mask to define which areas will be excluded for this
    period. Given as a binary string.
- `fossil`: Specify constraints on inner nodes. There are currently four different modes supported:
  - `include`: Constrain a node to contain at least the areas specified. For 
    example, `fossil include foo = 011` will constrain the node `foo` 
    (specified with `mrca`) to the distributions `111` and `011`.
  - `exclude`: Constrain a node to contain _none_ of the areas specified. For 
    example, `fossil exclude foo = 001` will exclude the distributions `001`,
    `011`, `101`, and `111` from node `foo`.
  - `fixed`: Constrain a node to be exactly the distribution specified. For example, `fossil fixed foo = 011` will
    constrain the node `foo` to the distribution `011`.
  - `node`: This option is an alias for `include`.
- `logfile`: Specify a file to log the output of the program (the log messages).
- `output-type`: One of either `json` or `csv`:
  - `json`: Compile all the results into one json file, who's contents are 
    detailed below. If your tree is large, this has a chance of running out of
    memory, particularly when computing splits.
  - `csv`: Output results as a set of CSV files.
- `opt-method`: Specify the optimization algorithm used to optimize the rates.
  Available methods are:
  - `nelder-mead`: Default in the old `lagrange-cpp`. Use this for parity with
  the old version, but in general it is slow. 
  - `bfgs`: A quasi-newton method that is in practice quite fast. This is my
  recommended option when using a large number of periods (more than 4).
  - `bobyqa`: An alternative to the `nelder-mead` option above that is much
  faster for a large number of periods. If there are issues with `bfgs` try this
  instead.
- `allow-ambiguous`: Allow coercion to an ambiguous state in the case that the
  tip range is too large for the specified max areas. By default this is
  _enabled_. Set this to false to enable it.
- `lwr-threshold`: Limit output of states and splits to results which have a
  LWR above this value. For example `lwr-threshold = 0.5` will only report
  results with a LWR above `0.5`.

If any of the arguments to the options has spaces, then it needs
to be quoted. For example, if the tree file is `My Super Awesome Tree.nwk`, then
control file needs to have 
```
treefile = 'My Super Awesome Tree.nwk`
```
to be able to parse the tree. This holds true for taxa names in the `mrca`
option, `mrca` names, user specified region names, and filenames.

## Periods

At times, researchers might want to model the changing geography of a region
over the course of the phylogenetic tree. For this purpose, there are
**periods**. Periods break up the time along the tree, and give each section of
time:

- independent extinction and dispersion parameters;
- independent connectivities and distance penalties; and
- period wide area inclusion and exclusion masks.

To specify a period and all of the possible options, Lagrange-NG uses a "declare
and use" model. That means that periods must be declared before modification.
For example, the following period statements are valid.

```ini
period early
period early exclude = 010
```

But the following are not, because the period is "used" before it is "declared".


```ini
period late matrix = late_distance.csv
period late
```

As mentioned above, each period has several parameters that can be modified.
These are

- `start` and `end` times;
- `include` and `exclude` masks; and
- a `matrix` csv file.example/periods.3.periods.csv

Below is an example showing all of these being set.

```ini
period early
period early end = 1.0
period early matrix = early_matrix.csv

period middle
period middle start = 1.0
period middle end = 2.0
period middle include = 100

period late
period late start = 2.0
period exclude = 010
```

Here, we declare three periods `early`, `middle`, and `late`. Period `early`
will be used for range evolution between `0.0` and `1.0` because the default
start time is `0.0`. In addition `early` has a connectivity matrix which is
specified in the file `early_matrix.csv`. Period `middle` runs from times `1.0`
to `2.0`, and requires that all ranges considered will include the first area.
Finally, the `late` period starts at time `2.0`, and has an unbounded end time.
Additionally, `late` has an exclude mask which will exclude any range which
includes the second area.

Important: Period times are measured _backwards_ through time, where `0.0`
indicates the present.

### Example

```ini
period early
period early exclude = 010
period early end = 0.2

period middle
period middle start = 0.2
period middle end = 0.4

period late
period late start = 0.4
period late exclude = 010
```

### Adjustment Matrices

Adjustment matrices are given as CSV files, with each row specifying an
adjustment between two areas. For two areas `A` and `B`, the adjustment value
between them represents the ease at which a population can spread to `B`
when starting from `A`. Adjustments can take on any value (including negative
values), but I suggest that only positive values are used. Adjustment values
greater than `1.0` indicate that the link between the two areas is "easy", and
values less than `1.0` indicate the link is "difficult".

While performing inference, Lagrange-NG will modify the per region dispersion
rate based on the adjustment values in the matrix. The contribution of the link
between `A` and `B` to the dispersion rate is `Adj(A,B)^x` where `x` is the
distance penalty. The distance penalty is an inferred value, and can be positive
or negative. If the distance penalty is above `1.0`, adjustment values which are
"easy" become more easy, and adjustment values which are "difficult" become more
difficult. When the distance penalty is below `-1.0`, the opposite is true. As
this is the case, I would suggest normalizing the values such that they are all
either greater than or equal to `1.0` or between `0.0` and `1.0`.

Adjustment matrices are given as CSV files with a header. Below is an example of
a (symmetric) adjustment matrix. 

```csv
from, to, dist
A, B, 0.2
A, C, 0.4
B, C, 0.1
```

Adjustment matrices can be specified as either symmetric, or non-symmetric. In
the symmetric case, like the matrix above, the adjustment between `A` and `B` is
assumed to be symmetric. When specifying a symmetric matrix, the form provided
should be in "upper triangular" form. That is, for `n` areas, the matrix should
have `n-1` rows where the first area is in the `from` column, and it should have
`n-2` rows where the second area is in the `from` column, and so on. The rows do
not have to be in any order however.

In the non-symmetric case, the rows can be in any order, and they are not
subject to any constraints. For example, below is a non-symmetric adjustment
matrix where the rows have no particular order.

```csv
from, to, dist
A, B, 0.2
A, C, 0.4
C, A, 0.3
C, B, 0.3
B, C, 0.1
B, A, 0.3
```

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
