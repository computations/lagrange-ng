# Introduction

The new version of Lagrange, a tool to infer ancestral ranges and splits using the DEC model. For a detailed
explanation of the please see our publication in Systematic Biology.

# Building 

There is a helpful makefile which will generate and build the cmake directory. So, you can just run `make` to build.

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
lagrange <CONTROL FILE>
```

After computation, there will be 3 results files: `<treefile>.nodes.tre`, `<treefile>.results.json`, and
`<treefile>.scaled.tre`. The `<treefile>.nodes.tre` contains the input tree with nodes labeled according to the
`<treefile.results.json` indices. This is so that the entries in the results file can be matched to the inner nodes of
the tree. The `<treefile>.results.json` file contains the

## Control File Options

- `treefile`: Path to the newick file containing the phylogeny. Required.
- `datafile`: Path to the matrix in a phylip format. At the moment, only phylip is supported. Required.
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

## Expm Modes

Lagrange-ng implements 2 methods of computing the matrix exponential, which can have fairly significant implications for
results. The first mode is `pade`, which is the much slower but (in my experience) more accurate method of computing the
exponential. The second method is `krylov`, which uses Krylov subspaces, and is often 2 orders of magnitude faster than
the `pade` mode. However, it can be numerically unstable for this particular application.

Lagrange-ng implements a third mode, `adaptive`, which attempts to combine the strengths of the `krylov` and `pade`
modes. By default, `adpative` uses `krylov` mode for matrix exponential computation, but watches for easily detectable
errors. If such an error is found, then the computation for that particular matrix is reverted to the `pade` mode.
In my experience, this mode fixes all the numerical errors that I could generate. I would recommend to use `adaptive` 
in all cases, except when the results seem very suspect. In this case, the `pade` mode should be used. Please note that
the `pade` mode is very slow compared to `adaptive` or `krylov`, so only use this when a genuine numerical error is
suspected.

## Results JSON file

The results file (`<treefile>.results.json`) contains the results from the optimization and ancestral state/split
computation. There are 3 top level keys to the file:

- `attributes`: This contains an object with information about the dataset, such as the number of regions and taxa.
- `node-results`: Contains a list of the results for each node. Discussed in detail later.
- `params`: Contains the inferred values for `dispersion` and `extinction`.

As mentioned, `node-results` contains a list of results for each node. Each entry of the list has the following keys and
values:

- `number`: Contains the node number, as specified in the `<treefile>.nodes.tre`.
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
